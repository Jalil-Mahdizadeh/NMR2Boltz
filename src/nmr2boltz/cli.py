from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .benchmark import BenchmarkManifestError, run_benchmark
from .output import write_outputs
from .project import ProjectionSettings, project_document
from .star import StarDataError, parse_star_document
from .target import TargetValidationError, validate_report_against_target
from .topology import TopologyLibrary, TopologyResolutionError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nmr2boltz",
        description=(
            "Parse NEF or NMR-STAR distance restraints and conservatively project proton "
            "upper bounds onto directly bonded heavy atoms for Boltz-2 atom_contact constraints."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)
    convert = subparsers.add_parser("convert", help="Convert one NEF/NMR-STAR file.")
    convert.add_argument("input", type=Path, help="Input NEF or NMR-STAR file.")
    convert.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="Output directory (default: <input_stem>_nmr2boltz).",
    )
    convert.add_argument(
        "--format",
        choices=["auto", "nef", "nmr-star"],
        default="auto",
        help="Input format; auto-detected from loop categories by default.",
    )
    convert.add_argument(
        "--averaging",
        choices=["sum-r6", "mean-r6", "hard-or"],
        default="sum-r6",
        help=(
            "Atom-set interpretation. sum-r6 is conservative for unnormalized NOE r^-6 sums; "
            "mean-r6 and hard-or use factor 1."
        ),
    )
    convert.add_argument(
        "--projection-margin",
        type=float,
        default=0.0,
        metavar="ANGSTROM",
        help="Additional non-negative slack added after the triangle-inequality projection.",
    )
    convert.add_argument(
        "--missing-upper-policy",
        choices=["reject", "upper-linear", "target-plus-uncertainty", "target"],
        default="reject",
        help="How to handle rows without an explicit upper bound. Default is reject.",
    )
    convert.add_argument(
        "--pseudoatom-policy",
        choices=["reject", "atomset"],
        default="reject",
        help=(
            "Geometric Q/M pseudoatoms are rejected by default. atomset is an explicit approximation "
            "and must be reviewed by an NMR expert."
        ),
    )
    convert.add_argument(
        "--residue-map",
        type=Path,
        help=(
            "CSV/TSV/JSON mapping source_chain, source_sequence_code, source_residue_name to "
            "boltz_chain and boltz_residue_index."
        ),
    )
    convert.add_argument(
        "--target-yaml",
        type=Path,
        help=(
            "Validate every mapped chain, residue index, residue identity, exact endpoint, and "
            "union-alternative endpoint against the Boltz input YAML before writing constraints."
        ),
    )
    convert.add_argument(
        "--chain-map",
        action="append",
        default=[],
        metavar="SOURCE=BOLTZ",
        help="Map a source chain code to a Boltz chain ID; repeat as needed.",
    )
    convert.add_argument(
        "--allow-inferred-sequence-map",
        action="store_true",
        help=(
            "Infer residue order from restraint identifiers only when no sequence loop exists. "
            "This is off by default because a wrong index can silently target the wrong Boltz residue."
        ),
    )
    convert.add_argument(
        "--ccd",
        action="append",
        default=[],
        type=Path,
        metavar="PATH",
        help=(
            "External wwPDB Chemical Component Dictionary mmCIF file or directory for modified residues; "
            "repeat as needed."
        ),
    )
    convert.add_argument(
        "--bond-length-config",
        type=Path,
        help=(
            "Optional JSON overrides for conservative X-H bond-length upper envelopes. "
            "See docs/SCIENTIFIC_METHOD.md for the schema."
        ),
    )
    convert.add_argument(
        "--origin",
        action="append",
        default=[],
        help=(
            "Include only a restraint origin/type (case-insensitive), for example NOE. Repeat for multiple. "
            "By default all general distance restraints are processed."
        ),
    )
    convert.add_argument(
        "--boltz-min-distance",
        type=float,
        default=2.0,
        metavar="ANGSTROM",
        help=(
            "Minimum accepted by BoltzUI atom_contact; lower exact and union bounds are "
            "raised conservatively (default: 2.0)."
        ),
    )
    convert.add_argument(
        "--boltz-max-distance",
        type=float,
        default=20.0,
        metavar="ANGSTROM",
        help=(
            "Maximum accepted by BoltzUI atom_contact; larger exact bounds and complete "
            "union groups containing a larger alternative are rejected (default: 20.0)."
        ),
    )
    convert.add_argument(
        "--min-sequence-separation",
        type=int,
        default=0,
        metavar="N",
        help="Explicitly filter same-chain projected contacts with residue separation < N.",
    )
    convert.add_argument(
        "--exclude-intraresidue",
        action="store_true",
        help=(
            "Do not emit projected pairs within the same mapped residue. "
            "All-local and mixed local/nonlocal OR groups are excluded in full."
        ),
    )
    convert.add_argument(
        "--exclude-intrachain",
        action="store_true",
        help=(
            "Emit only contacts between different mapped Boltz polymer chain IDs. "
            "Mixed intrachain/inter-chain OR groups are excluded in full."
        ),
    )
    convert.add_argument(
        "--hypotheses",
        type=int,
        default=0,
        metavar="N",
        help=(
            "Generate up to N deterministic assignment-hypothesis YAML files by choosing one alternative "
            "from each unresolved OR group. These are hypotheses, not conservative conversions."
        ),
    )
    convert.add_argument("--seed", type=int, default=0, help="Random seed for hypothesis selection.")
    convert.add_argument(
        "--strict",
        action="store_true",
        help="Exit with status 3 if any ambiguity or rejection remains after writing the audit outputs.",
    )
    benchmark = subparsers.add_parser(
        "benchmark",
        help="Run a versioned, checksum-aware conversion benchmark manifest.",
    )
    benchmark.add_argument("manifest", type=Path, help="Benchmark manifest YAML (schema version 1).")
    benchmark.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="Output directory (default: <manifest_stem>_results).",
    )
    benchmark.add_argument(
        "--case",
        action="append",
        default=[],
        metavar="ID",
        help="Run only one named case; repeat to select multiple cases.",
    )
    return parser


def parse_chain_map(values: list[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Invalid --chain-map {value!r}; expected SOURCE=BOLTZ.")
        source, target = value.split("=", 1)
        source = source.strip()
        target = target.strip()
        if not source or not target:
            raise ValueError(f"Invalid --chain-map {value!r}; both sides must be non-empty.")
        result[source] = target
    return result


def command_convert(args: argparse.Namespace) -> int:
    if not args.input.is_file():
        raise ValueError(f"Input file does not exist: {args.input}")
    if args.hypotheses < 0:
        raise ValueError("--hypotheses must be non-negative.")
    chain_map = parse_chain_map(args.chain_map)
    origins = {value.strip().lower() for value in args.origin if value.strip()} or None
    parsed = parse_star_document(
        args.input,
        format_hint=args.format,
        missing_upper_policy=args.missing_upper_policy,
        residue_map_path=args.residue_map,
        chain_map=chain_map,
        allow_inferred_sequence_map=args.allow_inferred_sequence_map,
        origins=origins,
    )
    topology = TopologyLibrary(args.ccd, bond_length_config=args.bond_length_config)
    settings = ProjectionSettings(
        averaging_policy=args.averaging,
        projection_margin=args.projection_margin,
        pseudoatom_policy=args.pseudoatom_policy,
        boltz_min_distance=args.boltz_min_distance,
        boltz_max_distance=args.boltz_max_distance,
        min_sequence_separation=args.min_sequence_separation,
        include_intraresidue=not args.exclude_intraresidue,
        include_intrachain=not args.exclude_intrachain,
    )
    parser_settings = {
        "format_hint": args.format,
        "missing_upper_policy": args.missing_upper_policy,
        "origin_filter": sorted(origins) if origins else "all",
        "residue_map": str(args.residue_map) if args.residue_map else None,
        "chain_map": chain_map,
        "allow_inferred_sequence_map": args.allow_inferred_sequence_map,
        "external_ccd_paths": [str(path) for path in args.ccd],
        "bond_length_config": str(args.bond_length_config) if args.bond_length_config else None,
    }
    report = project_document(
        parsed,
        input_file=str(args.input),
        topology_library=topology,
        settings=settings,
        parser_settings=parser_settings,
    )
    if args.target_yaml:
        target_validation = validate_report_against_target(report, args.target_yaml)
        target_validation.require_valid()
        report.target_validation = target_validation.to_dict()
    output_dir = args.output_dir or args.input.with_name(f"{args.input.stem}_nmr2boltz")
    written = write_outputs(
        report,
        output_dir,
        hypothesis_count=args.hypotheses,
        random_seed=args.seed,
    )
    print(f"Detected format: {report.detected_format}")
    print(f"Restraint groups read: {report.statistics['restraint_groups_read']}")
    print(f"Exact atom constraints emitted: {len(report.emitted_constraints)}")
    print(f"Atom-contact union groups emitted: {len(report.ambiguous_groups)}")
    print(f"Token contacts emitted: {len(report.token_constraints)}")
    if args.exclude_intrachain:
        print(
            "Intrachain restraint groups excluded: "
            f"{report.statistics['intrachain_groups_filtered']}"
        )
    print(f"Rejection records: {len(report.rejections)}")
    if report.target_validation:
        print(
            "Boltz target validation: PASS "
            f"({report.target_validation['checked_sequence_records']} sequence records; "
            f"{report.target_validation['warning_count']} warning(s))"
        )
    print(f"Output directory: {output_dir.resolve()}")
    print(f"Files written: {len(written)}")
    if args.strict and (report.ambiguous_groups or report.rejections):
        return 3
    return 0


def command_benchmark(args: argparse.Namespace) -> int:
    output_dir = args.output_dir or args.manifest.with_name(f"{args.manifest.stem}_results")
    run = run_benchmark(
        args.manifest,
        output_dir,
        selected_cases=set(args.case) or None,
    )
    for case in run.cases:
        details = case.error or "; ".join(case.mismatches)
        suffix = f" - {details}" if details else ""
        print(f"{case.case_id}: {case.status.upper()}{suffix}")
    print(f"Benchmark: {run.passed} passed, {run.failed} failed")
    print(f"Summary: {(Path(run.output_directory) / 'benchmark_summary.json').resolve()}")
    return 0 if run.failed == 0 else 4


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "convert":
            return command_convert(args)
        if args.command == "benchmark":
            return command_benchmark(args)
        parser.error(f"Unknown command: {args.command}")
    except (
        BenchmarkManifestError,
        ValueError,
        StarDataError,
        TargetValidationError,
        TopologyResolutionError,
    ) as exc:
        print(f"nmr2boltz: error: {exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
