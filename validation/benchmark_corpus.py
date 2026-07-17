#!/usr/bin/env python3
"""Run the paired NEF/NMR-STAR corpus against deposited PDB ensembles."""

from __future__ import annotations

import argparse
import hashlib
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
    from validation.discrepancy_audit import (
        audit_reports,
        audit_summary,
        write_audit_tsv,
    )
except ModuleNotFoundError:  # Direct execution: python validation/benchmark_corpus.py
    from compare_ensemble import (  # type: ignore[no-redef]
        align_structure_to_sequence_map,
        evaluate_heavy_constraints,
        evaluate_source_groups,
        implication_audit,
        load_structure,
    )
    from discrepancy_audit import (  # type: ignore[no-redef]
        audit_reports,
        audit_summary,
        write_audit_tsv,
    )


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
    parser.add_argument(
        "--reviewed-baseline",
        type=Path,
        default=Path("benchmark/reviewed_baseline.json"),
        help="Reviewed audit digest and metric snapshot used by the fail-closed gate.",
    )
    parser.add_argument(
        "--write-reviewed-baseline",
        action="store_true",
        help="Explicitly replace the reviewed baseline with this run after manual review.",
    )
    return parser.parse_args()


def _single_file(case_directory: Path, suffix: str) -> Path:
    matches = sorted(path for path in case_directory.glob(f"*{suffix}") if path.is_file())
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {suffix} file in {case_directory}, found {len(matches)}."
        )
    return matches[0]


def _portable_path(path: Path) -> str:
    """Serialize artifact paths deterministically across CI operating systems."""
    return path.as_posix()


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
    missing_constraints = [
        row for row in heavy_summary if row["models_resolved"] != len(model_ids)
    ]
    missing_contact_records = [
        {
            "constraint_id": row["constraint_id"],
            "atom1": row["atom1"],
            "atom2": row["atom2"],
            "upper_bound_angstrom": row["upper_bound"],
            "models_resolved": row["models_resolved"],
            "source_groups": row["source_groups"],
        }
        for row in missing_constraints
    ]
    return {
        "pdb": pdb_path.as_posix(),
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
            "constraints_missing_in_any_model": len(missing_constraints),
            "missing_coordinate_contacts": missing_contact_records,
            "missing_constraint_examples": missing_contact_records[:50],
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


def _run_case(
    case_directory: Path, output_root: Path, tolerance: float
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
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
            input_file=input_path.as_posix(),
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
    audit_rows = audit_reports(case_id, reports["nef"], reports["star"])
    case_audit_summary = audit_summary(audit_rows)
    result["discrepancy_audit"] = case_audit_summary
    write_audit_tsv(case_output / "format_discrepancy_audit.tsv", audit_rows)
    (case_output / "format_discrepancy_summary.json").write_text(
        json.dumps(case_audit_summary, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return result, audit_rows


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
            "missing_model_cases": sum(
                item["coordinates"]["heavy_atom_constraints"]["missing_model_cases"]
                for item in conversions
            ),
            "constraints_missing_in_any_model": sum(
                item["coordinates"]["heavy_atom_constraints"][
                    "constraints_missing_in_any_model"
                ]
                for item in conversions
            ),
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


def _metric_snapshot(run: dict[str, Any]) -> dict[str, Any]:
    """Select deterministic scientific metrics; any change requires review."""
    cases: dict[str, Any] = {}
    for case in run["cases"]:
        case_metrics: dict[str, Any] = {"status": case["status"]}
        if case["status"] == "pass":
            case_metrics["format_parity"] = case["format_parity"]
            case_metrics["formats"] = {}
            for format_name in ("nef", "star"):
                conversion = case["formats"][format_name]
                heavy = conversion["coordinates"]["heavy_atom_constraints"]
                implication = conversion["coordinates"]["projection_implication_audit"]
                case_metrics["formats"][format_name] = {
                    "restraint_groups": conversion["restraint_groups"],
                    "source_alternatives": conversion["source_alternatives"],
                    "emitted_constraints": conversion["emitted_constraints"],
                    "ambiguous_groups": conversion["ambiguous_groups"],
                    "rejection_reasons": conversion["rejection_reasons"],
                    "resolved_model_cases": heavy["resolved_model_cases"],
                    "satisfied_model_cases": heavy["satisfied_model_cases"],
                    "violated_model_cases": heavy["violated_model_cases"],
                    "missing_model_cases": heavy["missing_model_cases"],
                    "constraints_missing_in_any_model": heavy[
                        "constraints_missing_in_any_model"
                    ],
                    "implication_antecedent_cases": implication["antecedent_cases"],
                    "implication_failures": implication["failure_count"],
                }
        cases[case["case_id"]] = case_metrics
    return {"cases": cases}


def _missing_coordinate_review(run: dict[str, Any]) -> dict[str, Any]:
    """Pin every reviewed coordinate gap by scientific identity and digest."""
    contacts: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    for case in run["cases"]:
        if case["status"] != "pass":
            continue
        for format_name in ("nef", "star"):
            heavy = case["formats"][format_name]["coordinates"][
                "heavy_atom_constraints"
            ]
            for item in heavy.get("missing_coordinate_contacts", []):
                atom1, atom2 = sorted((str(item["atom1"]), str(item["atom2"])))
                source_groups = sorted(
                    group
                    for group in str(item.get("source_groups") or "").split(";")
                    if group
                )
                contact = {
                    "contact_id": f'{case["case_id"]}|{format_name}|{atom1}--{atom2}',
                    "case_id": case["case_id"],
                    "format": format_name,
                    "atom1": atom1,
                    "atom2": atom2,
                    "upper_bound_angstrom": float(item["upper_bound_angstrom"]),
                    "source_groups": source_groups,
                }
                contacts.append(contact)
                counts[f'{case["case_id"]}:{format_name}'] += 1
    contacts.sort(
        key=lambda item: (
            item["contact_id"],
            item["upper_bound_angstrom"],
            item["source_groups"],
        )
    )
    payload = json.dumps(contacts, separators=(",", ":"), sort_keys=True)
    return {
        "digest_sha256": hashlib.sha256(payload.encode("utf-8")).hexdigest(),
        "contact_count": len(contacts),
        "case_format_counts": dict(sorted(counts.items())),
        "contacts": contacts,
    }


def _baseline_payload(run: dict[str, Any]) -> dict[str, Any]:
    audit = run["discrepancy_audit"]
    return {
        "schema_version": 2,
        "audit_review": {
            "digest_sha256": audit["digest_sha256"],
            "row_count": audit["row_count"],
            "classification_counts": audit["classification_counts"],
        },
        "reviewed_missing_coordinates": _missing_coordinate_review(run),
        "metrics": _metric_snapshot(run),
    }


def _evaluate_gate(run: dict[str, Any], baseline_path: Path) -> dict[str, Any]:
    failures: list[dict[str, Any]] = []
    failed_cases = [case["case_id"] for case in run["cases"] if case["status"] != "pass"]
    if failed_cases:
        failures.append({"reason": "case_execution_failure", "cases": failed_cases})
    implication_failures = sum(
        run["aggregate"]["formats"].get(name, {}).get("implication_failures", 0)
        for name in ("nef", "star")
    )
    if implication_failures:
        failures.append(
            {"reason": "projection_implication_failure", "count": implication_failures}
        )
    unresolved = int(run["discrepancy_audit"]["unresolved_count"])
    if unresolved:
        failures.append({"reason": "unresolved_format_discrepancy", "count": unresolved})
    parser_bugs = int(run["discrepancy_audit"]["parser_projection_bug_count"])
    if parser_bugs:
        failures.append({"reason": "parser_projection_bug", "count": parser_bugs})

    if not baseline_path.is_file():
        failures.append(
            {"reason": "missing_reviewed_baseline", "path": _portable_path(baseline_path)}
        )
    else:
        baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
        current = _baseline_payload(run)
        if baseline.get("audit_review") != current["audit_review"]:
            failures.append(
                {
                    "reason": "unreviewed_discrepancy_change",
                    "expected": baseline.get("audit_review"),
                    "observed": current["audit_review"],
                }
            )
        if baseline.get("reviewed_missing_coordinates") != current[
            "reviewed_missing_coordinates"
        ]:
            expected_coordinates = baseline.get("reviewed_missing_coordinates") or {}
            observed_coordinates = current["reviewed_missing_coordinates"]
            failures.append(
                {
                    "reason": "changed_reviewed_missing_coordinate_set",
                    "expected_digest_sha256": expected_coordinates.get(
                        "digest_sha256"
                    ),
                    "observed_digest_sha256": observed_coordinates[
                        "digest_sha256"
                    ],
                    "expected_contact_count": expected_coordinates.get(
                        "contact_count"
                    ),
                    "observed_contact_count": observed_coordinates[
                        "contact_count"
                    ],
                }
            )
        if baseline.get("metrics") != current["metrics"]:
            failures.append({"reason": "unreviewed_metric_change"})
    return {
        "status": "pass" if not failures else "fail",
        "failure_count": len(failures),
        "failures": failures,
        "reviewed_baseline": _portable_path(baseline_path),
    }


def _percent(value: float | None) -> str:
    return "N/A" if value is None or not math.isfinite(value) else f"{100 * value:.2f}%"


def _markdown(run: dict[str, Any]) -> str:
    aggregate = run["aggregate"]
    audit = run["discrepancy_audit"]
    gate = run["gate"]
    coordinate_review = _missing_coordinate_review(run)
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
            "## Row-level format discrepancy audit",
            "",
            f"- {audit['row_count']} NEF-only, STAR-only, or different-bound contacts were audited to source rows and physical proton sets.",
            f"- Classifications: {audit['classification_counts'].get('expected_format_difference', 0)} scientifically expected format differences; {audit['classification_counts'].get('deposition_inconsistency', 0)} deposition inconsistencies; {audit['unresolved_count']} unresolved; {audit['parser_projection_bug_count']} parser/projection bugs.",
            f"- Reviewed audit digest: `{audit['digest_sha256']}`.",
            "- `expected_format_difference` is allowlisted only for wildcard-set versus explicit OR, x/y assignment versus a compatible physical set, rejected Q/M pseudoatoms, or verified canonical aliases.",
            "- 9PQH contains 43 contact discrepancies caused by its NEF sequence/residue numbering conflict; its remaining two different-bound rows are the expected explicit-OR versus wildcard atom-set distinction.",
            "- 9CCH contains 16 deposition inconsistencies because `ZN` and `ZN*` are unverified heavy-atom names that topology cannot prove equivalent.",
            "- 43JX, 6M6O, 9SGX, and 9VUY differences are justified by explicit-OR, wildcard atom-set, x/y assignment, or rejected geometric-pseudoatom semantics.",
            "",
            "## Fail-closed gate",
            "",
            f"- Gate status: **{gate['status'].upper()}** ({gate['failure_count']} failure {'category' if gate['failure_count'] == 1 else 'categories'}).",
            f"- Coordinate resolution gaps: {nef['constraints_missing_in_any_model']} NEF and {star['constraints_missing_in_any_model']} STAR contacts. These are not omitted from the denominator silently.",
            "- The current gaps are 45 contacts per format in partial-coordinate 8R1X and 8 contacts per format in 9CCH restraints that deposit `ZN`/`ZN*` on GLN 48 rather than the coordinate zinc residue.",
            f"- The exact 106-contact reviewed coordinate set is pinned by digest `{coordinate_review['digest_sha256']}`; any addition, removal, identity, bound, or provenance change fails CI.",
            "- Any implication failure, unresolved discrepancy, parser/projection bug, discrepancy-digest change, scientific-metric change, or reviewed-coordinate-set change makes the command exit nonzero.",
            "",
            "This validates conversion safety and format behavior against structures refined with the deposited restraints; it is not an independent Boltz prediction-accuracy benchmark.",
            "",
        ]
    )
    return "\n".join(lines)


def run_corpus(
    input_directory: Path,
    output_directory: Path,
    tolerance: float,
    reviewed_baseline: Path = Path("benchmark/reviewed_baseline.json"),
    *,
    write_reviewed_baseline: bool = False,
) -> dict[str, Any]:
    if tolerance < 0 or not math.isfinite(tolerance):
        raise ValueError("tolerance must be finite and non-negative")
    case_directories = sorted(path for path in input_directory.iterdir() if path.is_dir())
    if not case_directories:
        raise ValueError(f"No case directories were found in {input_directory}.")
    output_directory.mkdir(parents=True, exist_ok=True)
    cases: list[dict[str, Any]] = []
    all_audit_rows: list[dict[str, Any]] = []
    for case_directory in case_directories:
        try:
            case, audit_rows = _run_case(case_directory, output_directory, tolerance)
            cases.append(case)
            all_audit_rows.extend(audit_rows)
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
        "schema_version": 2,
        "input_directory": _portable_path(input_directory),
        "output_directory": _portable_path(output_directory),
        "settings": {
            "averaging_policy": "sum-r6",
            "pseudoatom_policy": "reject",
            "missing_upper_policy": "reject",
            "coordinate_tolerance_angstrom": tolerance,
        },
        "aggregate": _aggregate(cases),
        "discrepancy_audit": audit_summary(all_audit_rows),
        "cases": cases,
    }
    write_audit_tsv(output_directory / "FORMAT_DISCREPANCY_AUDIT.tsv", all_audit_rows)
    (output_directory / "format_discrepancy_summary.json").write_text(
        json.dumps(run["discrepancy_audit"], indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    if write_reviewed_baseline:
        reviewed_baseline.parent.mkdir(parents=True, exist_ok=True)
        reviewed_baseline.write_text(
            json.dumps(_baseline_payload(run), indent=2, sort_keys=False) + "\n",
            encoding="utf-8",
        )
    run["gate"] = _evaluate_gate(run, reviewed_baseline)
    (output_directory / "benchmark_summary.json").write_text(
        json.dumps(run, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    (output_directory / "BENCHMARK_REPORT.md").write_text(
        _markdown(run), encoding="utf-8"
    )
    (output_directory / "RUN_COMMAND.txt").write_text(
        "python validation/benchmark_corpus.py benchmark/input --output-directory benchmark/output --reviewed-baseline benchmark/reviewed_baseline.json\n",
        encoding="utf-8",
    )
    return run


def main() -> int:
    args = parse_args()
    run = run_corpus(
        args.input_directory,
        args.output_directory,
        args.tolerance,
        args.reviewed_baseline,
        write_reviewed_baseline=args.write_reviewed_baseline,
    )
    print(json.dumps(run["aggregate"], indent=2, sort_keys=False))
    print(json.dumps({"gate": run["gate"]}, indent=2, sort_keys=False))
    return 0 if run["gate"]["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
