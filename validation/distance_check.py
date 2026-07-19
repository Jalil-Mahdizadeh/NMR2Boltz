#!/usr/bin/env python3
"""Write per-entry NEF/STAR exact-contact bounds beside PDB-model distances."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import yaml

from nmr2boltz.output import _outward_decimal

try:
    from validation.compare_ensemble import (
        align_structure_to_sequence_map,
        distance,
        load_structure,
    )
except ModuleNotFoundError:  # Direct execution: python validation/distance_check.py
    from compare_ensemble import (  # type: ignore[no-redef]
        align_structure_to_sequence_map,
        distance,
        load_structure,
    )


AtomKey = tuple[str, int, str]
PairKey = tuple[AtomKey, AtomKey]
Point = tuple[float, float, float]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_directory",
        nargs="?",
        type=Path,
        default=Path("benchmark/input"),
        help="Directory containing one subdirectory per PDB ID.",
    )
    parser.add_argument(
        "--conversion-output",
        type=Path,
        default=Path("benchmark/output"),
        help="Existing corpus output containing NEF and STAR conversion artifacts.",
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        type=Path,
        default=Path("benchmark/distance_check"),
        help="Destination for the per-entry CSV files and summary.",
    )
    parser.add_argument(
        "--expected-case-count",
        type=int,
        default=12,
        help="Fail unless exactly this many benchmark entry directories are present.",
    )
    parser.add_argument(
        "--mapping-tolerance",
        type=float,
        default=1e-6,
        help="Maximum NEF/STAR coordinate-distance disagreement for a common pair (A).",
    )
    parser.add_argument(
        "--satisfaction-tolerance",
        type=float,
        default=0.01,
        help="Tolerance used only for the summary satisfaction comparison (A).",
    )
    return parser.parse_args(argv)


def _single_file(directory: Path, suffix: str) -> Path:
    matches = sorted(path for path in directory.glob(f"*{suffix}") if path.is_file())
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {suffix} file in {directory}, found {len(matches)}."
        )
    return matches[0]


def _atom_key(payload: dict[str, Any]) -> AtomKey:
    return (
        str(payload["chain"]),
        int(payload["residue_index"]),
        str(payload["atom_name"]),
    )


def _pair_key(atom1: dict[str, Any], atom2: dict[str, Any]) -> PairKey:
    return tuple(sorted((_atom_key(atom1), _atom_key(atom2))))  # type: ignore[return-value]


def _pair_text(pair: PairKey) -> str:
    return "--".join(":".join(map(str, atom)) for atom in pair)


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}.")
    return payload


def _report_bounds(report: dict[str, Any], path: Path) -> dict[PairKey, float]:
    bounds: dict[PairKey, float] = {}
    for constraint in report.get("emitted_constraints", []):
        pair = _pair_key(constraint["atom1"], constraint["atom2"])
        bound = float(constraint["max_distance"])
        if not math.isfinite(bound) or bound < 0:
            raise ValueError(f"Invalid exact-contact bound for {_pair_text(pair)} in {path}.")
        if pair in bounds:
            raise ValueError(f"Duplicate exact-contact pair {_pair_text(pair)} in {path}.")
        bounds[pair] = bound
    return bounds


def _yaml_bounds(path: Path) -> dict[PairKey, float]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or set(payload) != {"constraints"}:
        raise ValueError(f"Unexpected exact-contact YAML root in {path}.")
    constraints = payload["constraints"]
    if not isinstance(constraints, list):
        raise ValueError(f"Expected a constraints list in {path}.")

    bounds: dict[PairKey, float] = {}
    for wrapper in constraints:
        if not isinstance(wrapper, dict) or set(wrapper) != {"atom_contact"}:
            raise ValueError(f"Unexpected exact-contact record in {path}.")
        contact = wrapper["atom_contact"]
        if not isinstance(contact, dict) or set(contact) != {
            "atom1",
            "atom2",
            "max_distance",
            "force",
        }:
            raise ValueError(f"Unexpected exact-contact schema in {path}.")
        if contact["force"] is not True:
            raise ValueError(f"Exact contact is not forced in {path}.")
        atom1 = {
            "chain": contact["atom1"][0],
            "residue_index": contact["atom1"][1],
            "atom_name": contact["atom1"][2],
        }
        atom2 = {
            "chain": contact["atom2"][0],
            "residue_index": contact["atom2"][1],
            "atom_name": contact["atom2"][2],
        }
        pair = _pair_key(atom1, atom2)
        bound = float(contact["max_distance"])
        if not math.isfinite(bound) or bound < 0:
            raise ValueError(f"Invalid exact-contact bound for {_pair_text(pair)} in {path}.")
        if pair in bounds:
            raise ValueError(f"Duplicate exact-contact pair {_pair_text(pair)} in {path}.")
        bounds[pair] = bound
    return bounds


def load_exact_bounds(format_directory: Path) -> tuple[dict[str, Any], dict[PairKey, float]]:
    report_path = format_directory / "conversion_report.json"
    yaml_path = format_directory / "atom_constraints_exact.yaml"
    report = _load_json(report_path)
    report_bounds = _report_bounds(report, report_path)
    yaml_bounds = _yaml_bounds(yaml_path)
    if set(report_bounds) != set(yaml_bounds):
        missing_yaml = sorted(set(report_bounds) - set(yaml_bounds))
        extra_yaml = sorted(set(yaml_bounds) - set(report_bounds))
        raise ValueError(
            f"Exact YAML/report pair mismatch in {format_directory}: "
            f"missing_yaml={len(missing_yaml)}, extra_yaml={len(extra_yaml)}."
        )
    for pair, raw_bound in report_bounds.items():
        serialized = float(_outward_decimal(raw_bound))
        if yaml_bounds[pair] != serialized:
            raise ValueError(
                f"Exact YAML/report bound mismatch for {_pair_text(pair)} in "
                f"{format_directory}: YAML={yaml_bounds[pair]:.6f}, "
                f"expected={serialized:.6f}."
            )
    return report, yaml_bounds


def _aligned_models(
    report: dict[str, Any],
    raw_models: list[dict[AtomKey, Point]],
    residue_models: list[dict[tuple[str, int], str]],
) -> list[dict[AtomKey, Point]]:
    aligned, _alignment = align_structure_to_sequence_map(
        report["sequence_map"], raw_models, residue_models
    )
    return aligned


def _pair_distances(
    pair: PairKey,
    models: list[dict[AtomKey, Point]],
) -> list[float | None]:
    atom1, atom2 = pair
    return [
        None
        if atom1 not in atoms or atom2 not in atoms
        else distance(atoms[atom1], atoms[atom2])
        for atoms in models
    ]


def _merge_distances(
    pair: PairKey,
    nef_values: list[float | None] | None,
    star_values: list[float | None] | None,
    tolerance: float,
) -> list[float | None]:
    if nef_values is None:
        return list(star_values or [])
    if star_values is None:
        return list(nef_values)
    if len(nef_values) != len(star_values):
        raise ValueError(f"Model-count mismatch for common pair {_pair_text(pair)}.")

    merged: list[float | None] = []
    for model_index, (nef_value, star_value) in enumerate(
        zip(nef_values, star_values), start=1
    ):
        if (nef_value is None) != (star_value is None):
            raise ValueError(
                f"NEF/STAR coordinate-resolution mismatch for {_pair_text(pair)}, "
                f"model {model_index}."
            )
        if nef_value is None:
            merged.append(None)
            continue
        assert star_value is not None
        if abs(nef_value - star_value) > tolerance:
            raise ValueError(
                f"NEF/STAR aligned distance mismatch for {_pair_text(pair)}, "
                f"model {model_index}: NEF={nef_value:.9f}, STAR={star_value:.9f}."
            )
        merged.append(nef_value)
    return merged


def _format_bound(value: float | None) -> str:
    return "" if value is None else format(_outward_decimal(value), ".6f")


def _format_distance(value: float | None) -> str:
    return "" if value is None else f"{value:.6f}"


def build_case(
    case_directory: Path,
    conversion_root: Path,
    output_directory: Path,
    *,
    mapping_tolerance: float,
    satisfaction_tolerance: float,
) -> dict[str, Any]:
    case_id = case_directory.name.upper()
    case_output = conversion_root / case_id
    nef_report, nef_bounds = load_exact_bounds(case_output / "nef")
    star_report, star_bounds = load_exact_bounds(case_output / "star")
    pdb_path = _single_file(case_directory, ".pdb")
    model_ids, raw_models, residue_models = load_structure(pdb_path)
    if not raw_models:
        raise ValueError(f"PDB file contains no models: {pdb_path}")
    if len(set(model_ids)) != len(model_ids):
        raise ValueError(f"PDB model identifiers are not unique in {pdb_path}.")
    if any(names != residue_models[0] for names in residue_models[1:]):
        raise ValueError(f"PDB models do not have identical residue identities: {pdb_path}")

    nef_models = _aligned_models(nef_report, raw_models, residue_models)
    star_models = _aligned_models(star_report, raw_models, residue_models)
    pairs = sorted(set(nef_bounds) | set(star_bounds))
    model_fields = [f"model_{model_id}_distance_A" for model_id in model_ids]
    if len(set(model_fields)) != len(model_fields):
        raise ValueError(f"PDB model columns are not unique in {pdb_path}.")

    rows: list[dict[str, str]] = []
    missing_distance_cells = 0
    pairs_missing_any_model = 0
    resolved_by_format = {"nef": 0, "star": 0}
    satisfied_by_format = {"nef": 0, "star": 0}
    for pair in pairs:
        nef_values = _pair_distances(pair, nef_models) if pair in nef_bounds else None
        star_values = _pair_distances(pair, star_models) if pair in star_bounds else None
        values = _merge_distances(pair, nef_values, star_values, mapping_tolerance)
        missing = sum(value is None for value in values)
        missing_distance_cells += missing
        pairs_missing_any_model += missing > 0
        row = {
            "heavy_atom_pair": _pair_text(pair),
            "nef_max_distance_A": _format_bound(nef_bounds.get(pair)),
            "star_max_distance_A": _format_bound(star_bounds.get(pair)),
        }
        row.update(
            {
                field: _format_distance(value)
                for field, value in zip(model_fields, values)
            }
        )
        rows.append(row)
        for format_name, bound in (
            ("nef", nef_bounds.get(pair)),
            ("star", star_bounds.get(pair)),
        ):
            if bound is None:
                continue
            resolved_values = [value for value in values if value is not None]
            resolved_by_format[format_name] += len(resolved_values)
            satisfied_by_format[format_name] += sum(
                value <= bound + satisfaction_tolerance for value in resolved_values
            )

    output_directory.mkdir(parents=True, exist_ok=True)
    csv_path = output_directory / f"{case_id}.csv"
    fieldnames = [
        "heavy_atom_pair",
        "nef_max_distance_A",
        "star_max_distance_A",
        *model_fields,
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    common = set(nef_bounds) & set(star_bounds)
    different_bounds = sum(
        abs(nef_bounds[pair] - star_bounds[pair]) > 1e-6 for pair in common
    )
    csv_digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    return {
        "pdb_id": case_id,
        "pdb_file": pdb_path.as_posix(),
        "csv_file": csv_path.as_posix(),
        "csv_sha256": csv_digest,
        "model_count": len(model_ids),
        "model_ids": model_ids,
        "row_count": len(rows),
        "common_pairs": len(common),
        "nef_only_pairs": len(set(nef_bounds) - set(star_bounds)),
        "star_only_pairs": len(set(star_bounds) - set(nef_bounds)),
        "common_pairs_with_different_bounds": different_bounds,
        "distance_cells": len(rows) * len(model_ids),
        "resolved_distance_cells": len(rows) * len(model_ids) - missing_distance_cells,
        "missing_distance_cells": missing_distance_cells,
        "pairs_missing_any_model": pairs_missing_any_model,
        "nef_resolved_pair_models": resolved_by_format["nef"],
        "nef_satisfied_pair_models": satisfied_by_format["nef"],
        "nef_satisfaction_fraction": (
            satisfied_by_format["nef"] / resolved_by_format["nef"]
            if resolved_by_format["nef"]
            else None
        ),
        "star_resolved_pair_models": resolved_by_format["star"],
        "star_satisfied_pair_models": satisfied_by_format["star"],
        "star_satisfaction_fraction": (
            satisfied_by_format["star"] / resolved_by_format["star"]
            if resolved_by_format["star"]
            else None
        ),
    }


def _percent(value: float | None) -> str:
    return "n/a" if value is None else f"{100.0 * value:.2f}%"


def _write_readme(
    path: Path,
    cases: list[dict[str, Any]],
    *,
    mapping_tolerance: float,
    satisfaction_tolerance: float,
) -> None:
    lines = [
        "# Exact-contact distance check",
        "",
        "Each `<PDB-ID>.csv` contains the union of the topology-validated exact",
        "heavy-atom pairs emitted from that entry's NEF and NMR-STAR inputs.",
        "A blank format-bound cell means that format did not emit the pair. A blank",
        "model-distance cell means that at least one endpoint is not observed in that",
        "deposited PDB model. Ambiguous OR/union groups are intentionally excluded:",
        "flattening their alternatives into independent restraints would change their",
        "scientific meaning.",
        "",
        "PDB author numbering is aligned independently through each format's generated",
        "sequence map to Boltz one-based residue indices. For common pairs, generation",
        f"fails if NEF- and STAR-driven aligned distances differ by more than {mapping_tolerance:g} A",
        "or if only one format resolves the coordinate. Executable YAML bounds are",
        "cross-checked against conversion provenance and retain conservative six-decimal",
        "outward rounding. Geometric PDB distances are reported to six decimals.",
        "",
        "The satisfaction columns below are descriptive comparisons using",
        f"`distance <= max_distance + {satisfaction_tolerance:g} A`; they are not model-prediction",
        "accuracy estimates.",
        "",
        "| PDB | models | rows | common | NEF only | STAR only | different bounds | missing distances | NEF satisfied | STAR satisfied |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for case in cases:
        lines.append(
            f"| {case['pdb_id']} | {case['model_count']} | {case['row_count']} | "
            f"{case['common_pairs']} | {case['nef_only_pairs']} | "
            f"{case['star_only_pairs']} | "
            f"{case['common_pairs_with_different_bounds']} | "
            f"{case['missing_distance_cells']} | "
            f"{_percent(case['nef_satisfaction_fraction'])} | "
            f"{_percent(case['star_satisfaction_fraction'])} |"
        )
    lines.extend(
        [
            "",
            "Reproduce from the repository root after running the corpus benchmark:",
            "",
            "```bash",
            "python validation/distance_check.py benchmark/input \\",
            "  --conversion-output benchmark/output \\",
            "  --output-directory benchmark/distance_check",
            "```",
            "",
            "`distance_check_summary.json` records model IDs, row counts, missing-coordinate",
            "counts, descriptive satisfaction totals, and a SHA-256 digest for every CSV.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def run(args: argparse.Namespace) -> dict[str, Any]:
    if args.expected_case_count < 1:
        raise ValueError("--expected-case-count must be positive.")
    for name in ("mapping_tolerance", "satisfaction_tolerance"):
        value = float(getattr(args, name))
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"--{name.replace('_', '-')} must be finite and non-negative.")
    case_directories = sorted(
        path for path in args.input_directory.iterdir() if path.is_dir()
    )
    if len(case_directories) != args.expected_case_count:
        raise ValueError(
            f"Expected {args.expected_case_count} benchmark entries in "
            f"{args.input_directory}, found {len(case_directories)}."
        )
    args.output_directory.mkdir(parents=True, exist_ok=True)
    cases = [
        build_case(
            case_directory,
            args.conversion_output,
            args.output_directory,
            mapping_tolerance=args.mapping_tolerance,
            satisfaction_tolerance=args.satisfaction_tolerance,
        )
        for case_directory in case_directories
    ]
    expected_csv_names = {f"{case['pdb_id']}.csv" for case in cases}
    for stale_path in args.output_directory.glob("*.csv"):
        if stale_path.name not in expected_csv_names:
            stale_path.unlink()

    summary = {
        "schema_version": 1,
        "method": {
            "scope": "exact emitted heavy-atom contacts only",
            "ambiguous_union_groups_included": False,
            "bound_serialization": "conservative outward rounding to six decimals",
            "distance_serialization": "Euclidean PDB-model distance to six decimals",
            "mapping_tolerance_angstrom": args.mapping_tolerance,
            "satisfaction_tolerance_angstrom": args.satisfaction_tolerance,
        },
        "case_count": len(cases),
        "cases": cases,
        "totals": {
            "rows": sum(case["row_count"] for case in cases),
            "distance_cells": sum(case["distance_cells"] for case in cases),
            "resolved_distance_cells": sum(
                case["resolved_distance_cells"] for case in cases
            ),
            "missing_distance_cells": sum(
                case["missing_distance_cells"] for case in cases
            ),
        },
    }
    summary_path = args.output_directory / "distance_check_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    _write_readme(
        args.output_directory / "README.md",
        cases,
        mapping_tolerance=args.mapping_tolerance,
        satisfaction_tolerance=args.satisfaction_tolerance,
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    summary = run(args)
    print(
        f"Wrote {summary['case_count']} CSV files with "
        f"{summary['totals']['rows']} exact heavy-atom pairs to "
        f"{args.output_directory}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
