from __future__ import annotations

import hashlib
import json
import re
import shutil
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .output import write_outputs
from .project import ProjectionSettings, project_document
from .star import parse_star_document
from .target import TargetValidationResult, validate_report_against_target
from .topology import TopologyLibrary


class BenchmarkManifestError(ValueError):
    """Raised when a benchmark manifest does not satisfy the versioned schema."""


_CASE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_OPTION_KEYS = {
    "format",
    "missing_upper_policy",
    "residue_map",
    "chain_map",
    "allow_inferred_sequence_map",
    "origins",
    "ccd",
    "bond_length_config",
    "averaging",
    "projection_margin",
    "pseudoatom_policy",
    "boltz_min_distance",
    "boltz_max_distance",
    "min_sequence_separation",
    "include_intraresidue",
}
_METRIC_KEYS = {
    "restraint_groups_read",
    "source_alternatives_read",
    "safe_groups_before_pair_deduplication",
    "emitted_constraints",
    "ambiguous_groups",
    "rejection_records",
    "sequence_records",
    "target_validation_errors",
    "target_validation_warnings",
}


@dataclass(frozen=True)
class FileSpec:
    path: Path
    sha256: str | None = None


@dataclass
class BenchmarkCaseResult:
    case_id: str
    status: str
    input_file: str | None = None
    target_file: str | None = None
    output_directory: str | None = None
    metrics: dict[str, int] = field(default_factory=dict)
    expected: dict[str, int] = field(default_factory=dict)
    mismatches: list[str] = field(default_factory=list)
    target_validation: dict[str, Any] | None = None
    error: str | None = None


@dataclass
class BenchmarkRunResult:
    manifest: str
    schema_version: int
    output_directory: str
    cases: list[BenchmarkCaseResult]

    @property
    def passed(self) -> int:
        return sum(case.status == "pass" for case in self.cases)

    @property
    def failed(self) -> int:
        return len(self.cases) - self.passed

    def to_dict(self) -> dict[str, Any]:
        return {
            "manifest": self.manifest,
            "schema_version": self.schema_version,
            "output_directory": self.output_directory,
            "case_count": len(self.cases),
            "passed": self.passed,
            "failed": self.failed,
            "status": "pass" if self.failed == 0 else "fail",
            "cases": [asdict(case) for case in self.cases],
        }


def _mapping(value: Any, description: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise BenchmarkManifestError(f"{description} must be a mapping.")
    return value


def _resolve_path(base: Path, value: Any, description: str) -> Path:
    text = str(value).strip() if value is not None else ""
    if not text:
        raise BenchmarkManifestError(f"{description} path must be non-empty.")
    path = Path(text)
    return path if path.is_absolute() else base / path


def _file_spec(base: Path, value: Any, description: str) -> FileSpec:
    payload = _mapping(value, description)
    unknown = set(payload) - {"path", "sha256"}
    if unknown:
        raise BenchmarkManifestError(
            f"{description} has unsupported field(s): {', '.join(sorted(unknown))}."
        )
    path = _resolve_path(base, payload.get("path"), description)
    digest = payload.get("sha256")
    if digest is not None:
        digest = str(digest).strip().lower()
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise BenchmarkManifestError(f"{description} sha256 must contain 64 hexadecimal digits.")
    return FileSpec(path=path, sha256=digest)


def _verify_file(spec: FileSpec, description: str) -> None:
    if not spec.path.is_file():
        raise BenchmarkManifestError(f"{description} does not exist: {spec.path}")
    if spec.sha256 is None:
        return
    digest = hashlib.sha256()
    with spec.path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    actual = digest.hexdigest()
    if actual != spec.sha256:
        raise BenchmarkManifestError(
            f"{description} checksum mismatch for {spec.path}: expected {spec.sha256}, found {actual}."
        )


def _optional_path(base: Path, options: dict[str, Any], key: str) -> Path | None:
    value = options.get(key)
    return _resolve_path(base, value, f"Option {key}") if value is not None else None


def _path_list(base: Path, options: dict[str, Any], key: str) -> list[Path]:
    value = options.get(key, [])
    if not isinstance(value, list):
        raise BenchmarkManifestError(f"Option {key} must be a list.")
    return [_resolve_path(base, item, f"Option {key}") for item in value]


def _string_list(options: dict[str, Any], key: str) -> list[str]:
    value = options.get(key, [])
    if not isinstance(value, list):
        raise BenchmarkManifestError(f"Option {key} must be a list.")
    result = [str(item).strip() for item in value]
    if any(not item for item in result):
        raise BenchmarkManifestError(f"Option {key} cannot contain empty values.")
    return result


def _bool_option(options: dict[str, Any], key: str, default: bool) -> bool:
    value = options.get(key, default)
    if not isinstance(value, bool):
        raise BenchmarkManifestError(f"Option {key} must be true or false.")
    return value


def _case_metrics(report: Any, target: TargetValidationResult | None) -> dict[str, int]:
    statistics = report.statistics
    return {
        "restraint_groups_read": int(statistics["restraint_groups_read"]),
        "source_alternatives_read": int(statistics["source_alternatives_read"]),
        "safe_groups_before_pair_deduplication": int(
            statistics["safe_groups_before_pair_deduplication"]
        ),
        "emitted_constraints": len(report.emitted_constraints),
        "ambiguous_groups": len(report.ambiguous_groups),
        "rejection_records": len(report.rejections),
        "sequence_records": len(report.sequence_map),
        "target_validation_errors": len(target.errors) if target else 0,
        "target_validation_warnings": len(target.warnings) if target else 0,
    }


def _expected_metrics(value: Any, case_id: str) -> dict[str, int]:
    payload = _mapping(value, f"Benchmark case {case_id} expected")
    unknown = set(payload) - _METRIC_KEYS
    if unknown:
        raise BenchmarkManifestError(
            f"Benchmark case {case_id} has unsupported expected metric(s): "
            f"{', '.join(sorted(unknown))}."
        )
    result: dict[str, int] = {}
    for key, raw in payload.items():
        if isinstance(raw, bool):
            raise BenchmarkManifestError(f"Expected metric {key} in case {case_id} must be an integer.")
        try:
            parsed = int(raw)
        except (TypeError, ValueError) as exc:
            raise BenchmarkManifestError(
                f"Expected metric {key} in case {case_id} must be an integer."
            ) from exc
        if parsed < 0 or parsed != raw:
            raise BenchmarkManifestError(
                f"Expected metric {key} in case {case_id} must be a non-negative integer."
            )
        result[key] = parsed
    return result


def _run_case(
    case: dict[str, Any],
    *,
    base: Path,
    output_root: Path,
) -> BenchmarkCaseResult:
    case_id = str(case.get("id", "")).strip()
    result = BenchmarkCaseResult(case_id=case_id, status="fail")
    try:
        unknown = set(case) - {"id", "restraints", "target", "options", "expected"}
        if unknown:
            raise BenchmarkManifestError(
                f"Benchmark case {case_id!r} has unsupported field(s): {', '.join(sorted(unknown))}."
            )
        if not _CASE_ID.fullmatch(case_id):
            raise BenchmarkManifestError(
                f"Invalid benchmark case ID {case_id!r}; use letters, digits, dot, underscore, or hyphen."
            )
        case_output = output_root / case_id
        result.output_directory = str(case_output)
        if case_output.exists():
            shutil.rmtree(case_output)
        restraints = _file_spec(base, case.get("restraints"), f"Case {case_id} restraints")
        target = (
            _file_spec(base, case["target"], f"Case {case_id} target")
            if case.get("target") is not None
            else None
        )
        _verify_file(restraints, f"Case {case_id} restraints")
        if target:
            _verify_file(target, f"Case {case_id} target")
        result.input_file = str(restraints.path)
        result.target_file = str(target.path) if target else None

        options = _mapping(case.get("options", {}), f"Benchmark case {case_id} options")
        option_unknown = set(options) - _OPTION_KEYS
        if option_unknown:
            raise BenchmarkManifestError(
                f"Benchmark case {case_id} has unsupported option(s): "
                f"{', '.join(sorted(option_unknown))}."
            )
        chain_map_raw = options.get("chain_map", {})
        chain_map_payload = _mapping(chain_map_raw, f"Benchmark case {case_id} chain_map")
        chain_map = {str(key): str(value) for key, value in chain_map_payload.items()}
        origins = {item.lower() for item in _string_list(options, "origins")} or None
        residue_map = _optional_path(base, options, "residue_map")
        ccd_paths = _path_list(base, options, "ccd")
        bond_length_config = _optional_path(base, options, "bond_length_config")

        parsed = parse_star_document(
            restraints.path,
            format_hint=str(options.get("format", "auto")),
            missing_upper_policy=str(options.get("missing_upper_policy", "reject")),
            residue_map_path=residue_map,
            chain_map=chain_map,
            allow_inferred_sequence_map=_bool_option(
                options, "allow_inferred_sequence_map", False
            ),
            origins=origins,
        )
        topology = TopologyLibrary(ccd_paths, bond_length_config=bond_length_config)
        settings = ProjectionSettings(
            averaging_policy=str(options.get("averaging", "sum-r6")),
            projection_margin=float(options.get("projection_margin", 0.0)),
            pseudoatom_policy=str(options.get("pseudoatom_policy", "reject")),
            boltz_min_distance=float(options.get("boltz_min_distance", 2.0)),
            boltz_max_distance=float(options.get("boltz_max_distance", 20.0)),
            min_sequence_separation=int(options.get("min_sequence_separation", 0)),
            include_intraresidue=_bool_option(options, "include_intraresidue", True),
        )
        parser_settings = {
            "benchmark_case": case_id,
            "format_hint": options.get("format", "auto"),
            "missing_upper_policy": options.get("missing_upper_policy", "reject"),
            "origin_filter": sorted(origins) if origins else "all",
            "residue_map": str(residue_map) if residue_map else None,
            "chain_map": chain_map,
            "allow_inferred_sequence_map": _bool_option(
                options, "allow_inferred_sequence_map", False
            ),
            "external_ccd_paths": [str(path) for path in ccd_paths],
            "bond_length_config": str(bond_length_config) if bond_length_config else None,
        }
        report = project_document(
            parsed,
            input_file=str(restraints.path),
            topology_library=topology,
            settings=settings,
            parser_settings=parser_settings,
        )
        target_result = validate_report_against_target(report, target.path) if target else None
        if target_result:
            report.target_validation = target_result.to_dict()

        expected = _expected_metrics(case.get("expected", {}), case_id)
        metrics = _case_metrics(report, target_result)
        mismatches = [
            f"{key}: expected {expected_value}, observed {metrics[key]}"
            for key, expected_value in expected.items()
            if metrics[key] != expected_value
        ]
        if target_result and target_result.errors:
            mismatches.append(
                f"target validation produced {len(target_result.errors)} error(s)"
            )

        result.metrics = metrics
        result.expected = expected
        result.mismatches = mismatches
        result.target_validation = target_result.to_dict() if target_result else None
        if not mismatches:
            write_outputs(report, case_output)
            result.status = "pass"
    except Exception as exc:  # The suite must record a failed case and continue with the corpus.
        result.error = f"{type(exc).__name__}: {exc}"
    return result


def run_benchmark(
    manifest_path: str | Path,
    output_directory: str | Path,
    *,
    selected_cases: set[str] | None = None,
) -> BenchmarkRunResult:
    manifest = Path(manifest_path)
    try:
        payload = yaml.safe_load(manifest.read_text(encoding="utf-8"))
    except OSError as exc:
        raise BenchmarkManifestError(f"Unable to read benchmark manifest {manifest}: {exc}") from exc
    except yaml.YAMLError as exc:
        raise BenchmarkManifestError(f"Unable to parse benchmark manifest {manifest}: {exc}") from exc
    root = _mapping(payload, "Benchmark manifest")
    unknown = set(root) - {"schema_version", "cases"}
    if unknown:
        raise BenchmarkManifestError(
            f"Benchmark manifest has unsupported field(s): {', '.join(sorted(unknown))}."
        )
    if root.get("schema_version") != 1:
        raise BenchmarkManifestError("Benchmark manifest schema_version must be 1.")
    cases = root.get("cases")
    if not isinstance(cases, list) or not cases:
        raise BenchmarkManifestError("Benchmark manifest requires a non-empty cases list.")
    if any(not isinstance(case, dict) for case in cases):
        raise BenchmarkManifestError("Every benchmark case must be a mapping.")

    case_ids = [str(case.get("id", "")).strip() for case in cases]
    duplicates = sorted({case_id for case_id in case_ids if case_ids.count(case_id) > 1})
    if duplicates:
        raise BenchmarkManifestError(
            f"Benchmark case IDs must be unique; duplicates: {', '.join(duplicates)}."
        )
    if selected_cases:
        unknown_selection = selected_cases - set(case_ids)
        if unknown_selection:
            raise BenchmarkManifestError(
                f"Selected benchmark case(s) not found: {', '.join(sorted(unknown_selection))}."
            )
        cases = [case for case in cases if str(case.get("id", "")).strip() in selected_cases]

    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    results = [_run_case(case, base=manifest.parent, output_root=output) for case in cases]
    run = BenchmarkRunResult(
        manifest=str(manifest),
        schema_version=1,
        output_directory=str(output),
        cases=results,
    )
    (output / "benchmark_summary.json").write_text(
        json.dumps(run.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return run
