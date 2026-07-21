#!/usr/bin/env python3
"""Write per-entry token-contact bounds beside PDB-model residue distances."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
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
except ModuleNotFoundError:  # Direct execution from the validation directory.
    from compare_ensemble import (  # type: ignore[no-redef]
        align_structure_to_sequence_map,
        distance,
        load_structure,
    )


TokenKey = tuple[str, int]
TokenPair = tuple[TokenKey, TokenKey]
AtomKey = tuple[str, int, str]
Point = tuple[float, float, float]
ResiduePoints = dict[TokenKey, list[Point]]


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
        default=Path("benchmark/distance_check_token"),
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


def _token_key(payload: dict[str, Any]) -> TokenKey:
    return str(payload["chain"]), int(payload["residue_index"])


def _pair_key(token1: dict[str, Any], token2: dict[str, Any]) -> TokenPair:
    pair = tuple(sorted((_token_key(token1), _token_key(token2))))
    return pair  # type: ignore[return-value]


def _pair_text(pair: TokenPair) -> str:
    return "--".join(":".join(map(str, token)) for token in pair)


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}.")
    return payload


def _validate_bound(pair: TokenPair, bound: float, path: Path) -> None:
    if not math.isfinite(bound) or not 4.0 <= bound <= 20.0:
        raise ValueError(
            f"Invalid token-contact bound for {_pair_text(pair)} in {path}: {bound}."
        )


def _report_bounds(report: dict[str, Any], path: Path) -> dict[TokenPair, float]:
    if "token_constraints" not in report:
        raise ValueError(f"Missing token_constraints in {path}.")
    constraints = report["token_constraints"]
    if not isinstance(constraints, list):
        raise ValueError(f"Expected token_constraints to be a list in {path}.")

    bounds: dict[TokenPair, float] = {}
    for constraint in constraints:
        token1 = _token_key(constraint["token1"])
        token2 = _token_key(constraint["token2"])
        pair = _pair_key(constraint["token1"], constraint["token2"])
        if token1 == token2:
            raise ValueError(f"Self token pair {_pair_text(pair)} in {path}.")
        if (token1, token2) != pair:
            raise ValueError(f"Non-canonical token pair {_pair_text(pair)} in {path}.")
        bound = float(constraint["max_distance"])
        _validate_bound(pair, bound, path)
        if pair in bounds:
            raise ValueError(f"Duplicate token pair {_pair_text(pair)} in {path}.")
        bounds[pair] = bound
    return bounds


def _yaml_token(value: Any, path: Path) -> TokenKey:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"Invalid token endpoint in {path}: {value!r}.")
    return str(value[0]), int(value[1])


def _yaml_bounds(path: Path) -> dict[TokenPair, float]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or set(payload) != {"constraints"}:
        raise ValueError(f"Unexpected token-contact YAML root in {path}.")
    constraints = payload["constraints"]
    if not isinstance(constraints, list):
        raise ValueError(f"Expected a constraints list in {path}.")

    bounds: dict[TokenPair, float] = {}
    for wrapper in constraints:
        if not isinstance(wrapper, dict) or set(wrapper) != {"contact"}:
            raise ValueError(f"Unexpected token-contact record in {path}.")
        contact = wrapper["contact"]
        if not isinstance(contact, dict) or set(contact) != {
            "token1",
            "token2",
            "max_distance",
            "force",
        }:
            raise ValueError(f"Unexpected token-contact schema in {path}.")
        if contact["force"] is not False:
            raise ValueError(f"Token contact is forced in {path}.")
        token1 = _yaml_token(contact["token1"], path)
        token2 = _yaml_token(contact["token2"], path)
        pair = tuple(sorted((token1, token2)))
        if token1 == token2:
            raise ValueError(f"Self token pair {_pair_text(pair)} in {path}.")
        if (token1, token2) != pair:
            raise ValueError(f"Non-canonical token pair {_pair_text(pair)} in {path}.")
        bound = float(contact["max_distance"])
        _validate_bound(pair, bound, path)
        if pair in bounds:
            raise ValueError(f"Duplicate token pair {_pair_text(pair)} in {path}.")
        bounds[pair] = bound
    return bounds


def load_token_bounds(
    format_directory: Path,
) -> tuple[dict[str, Any], dict[TokenPair, float]]:
    report_path = format_directory / "conversion_report.json"
    yaml_path = format_directory / "token_constraints.yaml"
    report = _load_json(report_path)
    report_bounds = _report_bounds(report, report_path)
    yaml_bounds = _yaml_bounds(yaml_path)
    if set(report_bounds) != set(yaml_bounds):
        missing_yaml = sorted(set(report_bounds) - set(yaml_bounds))
        extra_yaml = sorted(set(yaml_bounds) - set(report_bounds))
        raise ValueError(
            f"Token YAML/report pair mismatch in {format_directory}: "
            f"missing_yaml={len(missing_yaml)}, extra_yaml={len(extra_yaml)}."
        )
    for pair, raw_bound in report_bounds.items():
        serialized = float(_outward_decimal(raw_bound))
        if yaml_bounds[pair] != serialized:
            raise ValueError(
                f"Token YAML/report bound mismatch for {_pair_text(pair)} in "
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


def _residue_points(models: list[dict[AtomKey, Point]]) -> list[ResiduePoints]:
    indexed_models: list[ResiduePoints] = []
    for atoms in models:
        indexed: defaultdict[TokenKey, list[Point]] = defaultdict(list)
        for (chain, residue_index, _atom_name), point in atoms.items():
            indexed[(chain, residue_index)].append(point)
        indexed_models.append(dict(indexed))
    return indexed_models


def _pair_distances(
    pair: TokenPair,
    models: list[ResiduePoints],
) -> list[float | None]:
    token1, token2 = pair
    values: list[float | None] = []
    for residues in models:
        points1 = residues.get(token1)
        points2 = residues.get(token2)
        if not points1 or not points2:
            values.append(None)
            continue
        values.append(
            min(
                distance(point1, point2)
                for point1 in points1
                for point2 in points2
            )
        )
    return values


def _merge_distances(
    pair: TokenPair,
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
    nef_report, nef_bounds = load_token_bounds(case_output / "nef")
    star_report, star_bounds = load_token_bounds(case_output / "star")
    pdb_path = _single_file(case_directory, ".pdb")
    model_ids, raw_models, residue_models = load_structure(
        pdb_path, heavy_atoms_only=True
    )
    if not raw_models:
        raise ValueError(f"PDB file contains no models: {pdb_path}")
    if len(set(model_ids)) != len(model_ids):
        raise ValueError(f"PDB model identifiers are not unique in {pdb_path}.")
    if any(names != residue_models[0] for names in residue_models[1:]):
        raise ValueError(
            f"PDB models do not have identical residue identities: {pdb_path}"
        )

    nef_models = _residue_points(
        _aligned_models(nef_report, raw_models, residue_models)
    )
    star_models = _residue_points(
        _aligned_models(star_report, raw_models, residue_models)
    )
    pairs = sorted(set(nef_bounds) | set(star_bounds))
    model_fields = [
        f"model_{model_id}_min_heavy_atom_distance_A" for model_id in model_ids
    ]
    if len(set(model_fields)) != len(model_fields):
        raise ValueError(f"PDB model columns are not unique in {pdb_path}.")

    rows: list[dict[str, str]] = []
    missing_distance_cells = 0
    pairs_missing_any_model = 0
    resolved_by_format = {"nef": 0, "star": 0}
    satisfied_by_format = {"nef": 0, "star": 0}
    for pair in pairs:
        nef_values = _pair_distances(pair, nef_models) if pair in nef_bounds else None
        star_values = (
            _pair_distances(pair, star_models) if pair in star_bounds else None
        )
        values = _merge_distances(pair, nef_values, star_values, mapping_tolerance)
        missing = sum(value is None for value in values)
        missing_distance_cells += missing
        pairs_missing_any_model += missing > 0
        row = {
            "residue_pair": _pair_text(pair),
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
        "residue_pair",
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
        "# Token-contact distance check",
        "",
        "Each `<PDB-ID>.csv` contains the union of the ordinary token contacts",
        "emitted from that entry's NEF and NMR-STAR conversions. NEF and STAR",
        "bounds remain separate; a blank bound means that format did not emit the",
        "residue pair. Contacts may originate from exact atom restraints or from",
        "safely collapsed single-token-pair union restraints.",
        "",
        "For every deposited PDB model, the reported geometric value is the minimum",
        "Euclidean distance between any non-hydrogen atom in the first residue and",
        "any non-hydrogen atom in the second residue. A blank model cell means that",
        "one or both residues have no observed heavy atom in that model.",
        "",
        "PDB author numbering is aligned independently through each format's generated",
        "sequence map to Boltz one-based residue indices. For common pairs, generation",
        "fails if NEF- and STAR-driven distances differ by more than "
        f"{mapping_tolerance:g} A",
        "or if only one format resolves the coordinate. Token YAML is cross-checked",
        "against conversion-report provenance, including canonical unique pairs,",
        "4-20 A bounds, conservative six-decimal outward rounding, and `force: false`.",
        "Geometric PDB distances are reported to six decimals.",
        "",
        "The satisfaction columns below are descriptive comparisons using",
        "`minimum heavy-atom distance <= max_distance + "
        f"{satisfaction_tolerance:g} A`;",
        "they are not model-prediction accuracy estimates.",
        "",
        "| PDB | models | rows | common | NEF only | STAR only | different bounds "
        "| missing distances | NEF satisfied | STAR satisfied |",
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
            "mkdir -p benchmark/distance_check_token",
            "docker run --rm --network none --read-only \\",
            "  -v \"$PWD:/work:ro\" \\",
            "  -v \"$PWD/benchmark/distance_check_token:/output\" \\",
            "  --entrypoint python -w /work nmr2boltz:latest \\",
            "  validation/distance_check_token.py benchmark/input \\",
            "  --conversion-output benchmark/output --output-directory /output",
            "```",
            "",
            "`distance_check_token_summary.json` records model IDs, row counts,",
            "missing-coordinate counts, descriptive satisfaction totals, and a",
            "SHA-256 digest for every CSV.",
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
            raise ValueError(
                f"--{name.replace('_', '-')} must be finite and non-negative."
            )
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
            "scope": "all emitted ordinary token contacts",
            "source_kinds": ["exact", "collapsed_union"],
            "bound_serialization": "conservative outward rounding to six decimals",
            "distance_definition": (
                "minimum Euclidean distance between non-hydrogen atoms in each "
                "residue, per PDB model"
            ),
            "distance_serialization": "PDB-model distance to six decimals",
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
    summary_path = args.output_directory / "distance_check_token_summary.json"
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
        f"{summary['totals']['rows']} token residue pairs to "
        f"{args.output_directory}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
