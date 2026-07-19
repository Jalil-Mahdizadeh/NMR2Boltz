import hashlib
import json
from pathlib import Path

import yaml

from nmr2boltz.benchmark import run_benchmark


FIXTURES = Path(__file__).parent / "fixtures"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _manifest(tmp_path: Path, checksum: str) -> Path:
    target = tmp_path / "target.yaml"
    target.write_text(
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        "      sequence: VAYLG\n"
        "      msa: empty\n",
        encoding="utf-8",
    )
    restraints = FIXTURES / "example.nef"
    manifest = tmp_path / "benchmark.yaml"
    manifest.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "cases": [
                    {
                        "id": "fixture-nef",
                        "restraints": {"path": str(restraints), "sha256": checksum},
                        "target": {"path": "target.yaml", "sha256": _sha256(target)},
                        "options": {"averaging": "sum-r6"},
                        "expected": {
                            "sequence_records": 5,
                            "target_validation_errors": 0,
                        },
                    }
                ],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return manifest


def test_benchmark_manifest_runs_conversion_and_writes_auditable_summary(tmp_path):
    restraints = FIXTURES / "example.nef"
    manifest = _manifest(tmp_path, _sha256(restraints))
    output = tmp_path / "results"

    run = run_benchmark(manifest, output)

    assert run.passed == 1
    assert run.failed == 0
    assert (output / "fixture-nef" / "conversion_report.json").is_file()
    assert (output / "fixture-nef" / "atom_constraints_exact.yaml").is_file()
    assert (output / "fixture-nef" / "atom_constraints_union.yaml").is_file()
    assert not (output / "fixture-nef" / "boltz_constraints.yaml").exists()
    assert not (output / "fixture-nef" / "proposed_atom_contact_unions.yaml").exists()
    summary = json.loads((output / "benchmark_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "pass"
    assert summary["cases"][0]["metrics"]["emitted_constraints"] > 0


def test_benchmark_checksum_mismatch_is_recorded_as_case_failure(tmp_path):
    manifest = _manifest(tmp_path, "0" * 64)
    stale = tmp_path / "results" / "fixture-nef" / "atom_constraints_exact.yaml"
    stale.parent.mkdir(parents=True)
    stale.write_text("stale executable output\n", encoding="utf-8")

    run = run_benchmark(manifest, tmp_path / "results")

    assert run.failed == 1
    assert "checksum mismatch" in (run.cases[0].error or "")
    assert not stale.exists()


def test_benchmark_manifest_supports_exclude_intrachain(tmp_path):
    restraints = FIXTURES / "multichain.nef"
    manifest = tmp_path / "benchmark.yaml"
    manifest.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "cases": [
                    {
                        "id": "interchain-only",
                        "restraints": {
                            "path": str(restraints),
                            "sha256": _sha256(restraints),
                        },
                        "options": {"exclude_intrachain": True},
                        "expected": {
                            "emitted_constraints": 2,
                            "ambiguous_groups": 1,
                            "intrachain_groups_filtered": 3,
                            "mixed_chain_scope_groups_filtered": 1,
                        },
                    }
                ],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    output = tmp_path / "results"

    run = run_benchmark(manifest, output)

    assert run.passed == 1
    conversion = json.loads(
        (output / "interchain-only" / "conversion_report.json").read_text()
    )
    assert conversion["settings"]["include_intrachain"] is False
    assert conversion["settings"]["exclude_intrachain"] is True
    for constraint in conversion["emitted_constraints"]:
        assert constraint["atom1"]["chain"] != constraint["atom2"]["chain"]
    for group in conversion["ambiguous_groups"]:
        for alternative in group["alternatives"]:
            assert alternative["atom1"]["chain"] != alternative["atom2"]["chain"]
