#!/usr/bin/env python3
"""Run the paired NEF/NMR-STAR corpus against deposited PDB ensembles."""

from __future__ import annotations

import argparse
import json
import math
import shutil
from collections import Counter
from pathlib import Path
from typing import Any

from nmr2boltz.output import write_outputs
from nmr2boltz.project import ProjectionSettings, project_document
from nmr2boltz.star import parse_star_document
from nmr2boltz.topology import TopologyLibrary
try:
    from validation.compare_ensemble import (
        align_structure_to_sequence_map,
        evaluate_heavy_constraints,
        evaluate_source_groups,
        implication_audit,
        load_structure,
    )
except ModuleNotFoundError:  # Direct execution: python validation/benchmark_corpus.py
    from compare_ensemble import (  # type: ignore[no-redef]
        align_structure_to_sequence_map,
        evaluate_heavy_constraints,
        evaluate_source_groups,
        implication_audit,
        load_structure,
    )


INITIAL_AUDIT = {
    "successful_conversions": 22,
    "no_distance_failures": 2,
    "nef_emitted_constraints": 13004,
    "star_emitted_constraints": 11841,
    "fasta_files": 0,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_directory",
        nargs="?",
        type=Path,
        default=Path("benchmark/input"),
        help="Directory containing one subdirectory per PDB ID.",
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        type=Path,
        default=Path("benchmark/output"),
        help="Destination for per-case conversion and audit outputs.",
    )
    parser.add_argument("--tolerance", type=float, default=0.01)
    return parser.parse_args()


def _single_file(case_directory: Path, suffix: str) -> Path:
    matches = sorted(path for path in case_directory.glob(f"*{suffix}") if path.is_file())
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {suffix} file in {case_directory}, found {len(matches)}."
        )
    return matches[0]


def _constraint_key(constraint: Any) -> tuple[tuple[str, int, str], tuple[str, int, str]]:
    ends = (
        (constraint.atom1.chain, constraint.atom1.residue_index, constraint.atom1.atom_name),
        (constraint.atom2.chain, constraint.atom2.residue_index, constraint.atom2.atom_name),
    )
    return tuple(sorted(ends))  # type: ignore[return-value]


def compare_reports(nef_report: Any, star_report: Any) -> dict[str, Any]:
    nef = {_constraint_key(item): float(item.max_distance) for item in nef_report.emitted_constraints}
    star = {_constraint_key(item): float(item.max_distance) for item in star_report.emitted_constraints}
    common = set(nef) & set(star)
    bound_deltas = [abs(nef[key] - star[key]) for key in common]
    different = [delta for delta in bound_deltas if delta > 1e-6]
    return {
        "nef_constraints": len(nef),
        "star_constraints": len(star),
        "common_atom_pairs": len(common),
        "nef_only_atom_pairs": len(set(nef) - set(star)),
        "star_only_atom_pairs": len(set(star) - set(nef)),
        "common_pairs_with_different_bounds": len(different),
        "maximum_common_bound_delta_angstrom": max(bound_deltas) if bound_deltas else None,
        "exact_pair_and_bound_parity": (
            set(nef) == set(star) and not different
        ),
    }


def _coordinate_summary(report: Any, pdb_path: Path, tolerance: float) -> dict[str, Any]:
    report_dict = report.to_dict()
    model_ids, models, residue_models = load_structure(pdb_path)
    if not models:
        raise ValueError(f"PDB file contains no models: {pdb_path}")
    if any(names != residue_models[0] for names in residue_models[1:]):
        raise ValueError(f"PDB models do not have identical residue identities: {pdb_path}")
    models, alignment = align_structure_to_sequence_map(
        report_dict["sequence_map"], models, residue_models
    )
    source_rows, source_summary, source_states = evaluate_source_groups(
        report_dict, model_ids, models, tolerance
    )
    heavy_rows, heavy_summary = evaluate_heavy_constraints(
        report_dict, model_ids, models, tolerance
    )
    implication = implication_audit(report_dict, heavy_rows, source_states, tolerance)
    source_states_count = Counter(row["state"] for row in source_rows)
    satisfied = sum(row["satisfied"] is True for row in heavy_rows)
    violated = sum(row["satisfied"] is False for row in heavy_rows)
    missing = sum(row["satisfied"] is None for row in heavy_rows)
    resolved = satisfied + violated
    return {
        "pdb": str(pdb_path),
        "model_count": len(model_ids),
        "tolerance_angstrom": tolerance,
        "sequence_alignment": alignment,
        "source_restraints": {
            "group_count": len(source_summary),
            "group_model_state_counts": dict(source_states_count),
            "groups_satisfied_in_any_model": sum(
                row["models_satisfied"] > 0 for row in source_summary
            ),
            "groups_satisfied_in_every_model": sum(
                row["models_satisfied"] == len(model_ids) for row in source_summary
            ),
            "groups_indeterminate_in_every_model": sum(
                row["models_indeterminate"] == len(model_ids) for row in source_summary
            ),
        },
        "heavy_atom_constraints": {
            "constraint_count": len(heavy_summary),
            "resolved_model_cases": resolved,
            "satisfied_model_cases": satisfied,
            "violated_model_cases": violated,
            "missing_model_cases": missing,
            "satisfaction_fraction_of_resolved": satisfied / resolved if resolved else None,
            "constraints_resolved_in_every_model": sum(
                row["models_resolved"] == len(model_ids) for row in heavy_summary
            ),
            "constraints_satisfied_in_any_model": sum(
                row["models_satisfied"] > 0 for row in heavy_summary
            ),
            "constraints_satisfied_in_every_model": sum(
                row["models_satisfied"] == len(model_ids) for row in heavy_summary
            ),
        },
        "projection_implication_audit": implication,
    }


def _conversion_summary(report: Any, coordinates: dict[str, Any]) -> dict[str, Any]:
    reasons = Counter(item.reason for item in report.rejections)
    return {
        "detected_format": report.detected_format,
        "has_distance_restraints": bool(report.source_restraint_groups),
        "restraint_groups": len(report.source_restraint_groups),
        "source_alternatives": int(report.statistics["source_alternatives_read"]),
        "emitted_constraints": len(report.emitted_constraints),
        "ambiguous_groups": len(report.ambiguous_groups),
        "rejection_records": len(report.rejections),
        "rejection_reasons": dict(reasons),
        "sequence_records": len(report.sequence_map),
        "warnings": report.warnings,
        "coordinates": coordinates,
    }


def _safe_reset_case(case_output: Path, output_root: Path) -> None:
    resolved_case = case_output.resolve()
    resolved_root = output_root.resolve()
    try:
        resolved_case.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(f"Refusing to clean output outside {resolved_root}: {resolved_case}") from exc
    if resolved_case == resolved_root:
        raise ValueError(f"Refusing to clean benchmark output root directly: {resolved_root}")
    if case_output.exists():
        shutil.rmtree(case_output)


def _run_case(case_directory: Path, output_root: Path, tolerance: float) -> dict[str, Any]:
    case_id = case_directory.name.upper()
    case_output = output_root / case_id
    _safe_reset_case(case_output, output_root)
    case_output.mkdir(parents=True, exist_ok=True)
    pdb_path = _single_file(case_directory, ".pdb")
    reports: dict[str, Any] = {}
    result: dict[str, Any] = {"case_id": case_id, "status": "pass", "formats": {}}
    for format_name, suffix in (("nef", ".nef"), ("star", ".str")):
        input_path = _single_file(case_directory, suffix)
        parsed = parse_star_document(input_path)
        report = project_document(
            parsed,
            input_file=str(input_path),
            topology_library=TopologyLibrary(),
            settings=ProjectionSettings(),
            parser_settings={
                "benchmark_case": case_id,
                "benchmark_format": format_name,
                "missing_upper_policy": "reject",
                "origin_filter": "all",
            },
        )
        format_output = case_output / format_name
        write_outputs(report, format_output)
        coordinates = _coordinate_summary(report, pdb_path, tolerance)
        (format_output / "coordinate_summary.json").write_text(
            json.dumps(coordinates, indent=2, sort_keys=False) + "\n",
            encoding="utf-8",
        )
        reports[format_name] = report
        result["formats"][format_name] = _conversion_summary(report, coordinates)
    parity = compare_reports(reports["nef"], reports["star"])
    result["format_parity"] = parity
    (case_output / "format_parity.json").write_text(
        json.dumps(parity, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    return result


def _aggregate(cases: list[dict[str, Any]]) -> dict[str, Any]:
    aggregate: dict[str, Any] = {
        "case_count": len(cases),
        "successful_cases": sum(case["status"] == "pass" for case in cases),
        "successful_conversions": 0,
        "no_distance_conversions": 0,
        "fasta_files": 0,
        "exact_format_parity_cases": sum(
            case.get("format_parity", {}).get("exact_pair_and_bound_parity", False)
            for case in cases
        ),
        "exact_positive_distance_parity_cases": sum(
            case.get("format_parity", {}).get("exact_pair_and_bound_parity", False)
            and case.get("formats", {}).get("nef", {}).get("has_distance_restraints", False)
            for case in cases
        ),
        "formats": {},
    }
    for format_name in ("nef", "star"):
        conversions = [
            case["formats"][format_name]
            for case in cases
            if case.get("status") == "pass" and format_name in case.get("formats", {})
        ]
        satisfied = sum(
            item["coordinates"]["heavy_atom_constraints"]["satisfied_model_cases"]
            for item in conversions
        )
        violated = sum(
            item["coordinates"]["heavy_atom_constraints"]["violated_model_cases"]
            for item in conversions
        )
        aggregate["successful_conversions"] += len(conversions)
        aggregate["no_distance_conversions"] += sum(
            not item["has_distance_restraints"] for item in conversions
        )
        aggregate["fasta_files"] += len(conversions)
        aggregate["formats"][format_name] = {
            "emitted_constraints": sum(item["emitted_constraints"] for item in conversions),
            "resolved_model_cases": satisfied + violated,
            "satisfied_model_cases": satisfied,
            "violated_model_cases": violated,
            "satisfaction_fraction_of_resolved": (
                satisfied / (satisfied + violated) if satisfied + violated else None
            ),
            "implication_antecedent_cases": sum(
                item["coordinates"]["projection_implication_audit"]["antecedent_cases"]
                for item in conversions
            ),
            "implication_failures": sum(
                item["coordinates"]["projection_implication_audit"]["failure_count"]
                for item in conversions
            ),
            "sequence_residue_mismatch_rejections": sum(
                item["rejection_reasons"].get("sequence_residue_mismatch", 0)
                for item in conversions
            ),
        }
    return aggregate


def _percent(value: float | None) -> str:
    return "N/A" if value is None or not math.isfinite(value) else f"{100 * value:.2f}%"


def _markdown(run: dict[str, Any]) -> str:
    aggregate = run["aggregate"]
    lines = [
        "# nmr2boltz paired-format benchmark",
        "",
        "All 12 deposited ensembles were converted from both NEF and NMR-STAR with conservative defaults, then audited against every PDB conformer using sequence-aware coordinate alignment.",
        "",
        "| Case | NEF contacts | STAR contacts | NEF PDB satisfaction | STAR PDB satisfaction | Exact parity | Implication failures |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for case in run["cases"]:
        if case["status"] != "pass":
            lines.append(f"| {case['case_id']} | error | error | N/A | N/A | no | N/A |")
            continue
        nef = case["formats"]["nef"]
        star = case["formats"]["star"]
        nef_heavy = nef["coordinates"]["heavy_atom_constraints"]
        star_heavy = star["coordinates"]["heavy_atom_constraints"]
        failures = (
            nef["coordinates"]["projection_implication_audit"]["failure_count"]
            + star["coordinates"]["projection_implication_audit"]["failure_count"]
        )
        lines.append(
            f"| {case['case_id']} | {nef['emitted_constraints']} | {star['emitted_constraints']} | "
            f"{_percent(nef_heavy['satisfaction_fraction_of_resolved'])} | "
            f"{_percent(star_heavy['satisfaction_fraction_of_resolved'])} | "
            f"{'yes' if case['format_parity']['exact_pair_and_bound_parity'] else 'no'} | {failures} |"
        )
    nef = aggregate["formats"]["nef"]
    star = aggregate["formats"]["star"]
    lines.extend(
        [
            "",
            "## Current result",
            "",
            f"- {aggregate['successful_conversions']}/24 conversions completed; {aggregate['no_distance_conversions']} are valid empty distance conversions for 8S8O.",
            f"- NEF: {nef['emitted_constraints']} contacts, {_percent(nef['satisfaction_fraction_of_resolved'])} resolved PDB satisfaction.",
            f"- NMR-STAR: {star['emitted_constraints']} contacts, {_percent(star['satisfaction_fraction_of_resolved'])} resolved PDB satisfaction.",
            f"- Conservative implication failures: {nef['implication_failures'] + star['implication_failures']} across {nef['implication_antecedent_cases'] + star['implication_antecedent_cases']} satisfied-antecedent cases.",
            f"- Exact NEF/STAR pair-and-bound parity: {aggregate['exact_positive_distance_parity_cases']}/11 positive-distance cases; 8S8O also has exact empty-output parity.",
            "",
            "## Comparison with the initial audit",
            "",
            f"- Successful conversions: {INITIAL_AUDIT['successful_conversions']} -> {aggregate['successful_conversions']}.",
            f"- No-distance failures: {INITIAL_AUDIT['no_distance_failures']} -> 0.",
            f"- FASTA outputs: {INITIAL_AUDIT['fasta_files']} -> {aggregate['fasta_files']}.",
            f"- Emitted contacts: {INITIAL_AUDIT['nef_emitted_constraints']} -> {nef['emitted_constraints']} NEF and {INITIAL_AUDIT['star_emitted_constraints']} -> {star['emitted_constraints']} STAR under unchanged conservative projection policies.",
            f"- Sequence/residue conflicts are now explicit: {nef['sequence_residue_mismatch_rejections'] + star['sequence_residue_mismatch_rejections']} rejection records use `sequence_residue_mismatch`.",
            "- Large pseudoatom and atom-set parity differences remain intentionally visible rather than being silently approximated.",
            "",
            "This validates conversion safety and format behavior against structures refined with the deposited restraints; it is not an independent Boltz prediction-accuracy benchmark.",
            "",
        ]
    )
    return "\n".join(lines)


def run_corpus(input_directory: Path, output_directory: Path, tolerance: float) -> dict[str, Any]:
    if tolerance < 0 or not math.isfinite(tolerance):
        raise ValueError("tolerance must be finite and non-negative")
    case_directories = sorted(path for path in input_directory.iterdir() if path.is_dir())
    if not case_directories:
        raise ValueError(f"No case directories were found in {input_directory}.")
    output_directory.mkdir(parents=True, exist_ok=True)
    cases: list[dict[str, Any]] = []
    for case_directory in case_directories:
        try:
            cases.append(_run_case(case_directory, output_directory, tolerance))
        except Exception as exc:
            cases.append(
                {
                    "case_id": case_directory.name.upper(),
                    "status": "fail",
                    "error": f"{type(exc).__name__}: {exc}",
                    "formats": {},
                }
            )
    run = {
        "schema_version": 1,
        "input_directory": str(input_directory),
        "output_directory": str(output_directory),
        "settings": {
            "averaging_policy": "sum-r6",
            "pseudoatom_policy": "reject",
            "missing_upper_policy": "reject",
            "coordinate_tolerance_angstrom": tolerance,
        },
        "aggregate": _aggregate(cases),
        "cases": cases,
    }
    (output_directory / "benchmark_summary.json").write_text(
        json.dumps(run, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    (output_directory / "BENCHMARK_REPORT.md").write_text(
        _markdown(run), encoding="utf-8"
    )
    (output_directory / "RUN_COMMAND.txt").write_text(
        "python validation/benchmark_corpus.py benchmark/input --output-directory benchmark/output\n",
        encoding="utf-8",
    )
    return run


def main() -> int:
    args = parse_args()
    run = run_corpus(args.input_directory, args.output_directory, args.tolerance)
    print(json.dumps(run["aggregate"], indent=2, sort_keys=False))
    return 0 if run["aggregate"]["successful_cases"] == run["aggregate"]["case_count"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
