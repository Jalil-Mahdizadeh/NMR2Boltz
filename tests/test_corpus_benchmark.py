from pathlib import Path
from types import SimpleNamespace

from nmr2boltz.model import BoltzAtom, EmittedConstraint
from validation.benchmark_corpus import (
    _baseline_payload,
    _evaluate_gate,
    _portable_path,
    compare_reports,
)
from validation.discrepancy_audit import audit_reports, audit_summary


def _constraint(atom1: BoltzAtom, atom2: BoltzAtom, distance: float) -> EmittedConstraint:
    return EmittedConstraint(
        atom1=atom1,
        atom2=atom2,
        max_distance=distance,
        source_groups=["g1"],
        raw_projected_distance=distance,
    )


def test_compare_reports_detects_pair_and_bound_differences():
    a = BoltzAtom("A", 1, "CA")
    b = BoltzAtom("A", 2, "CB")
    c = BoltzAtom("A", 3, "CG")
    nef = SimpleNamespace(emitted_constraints=[_constraint(a, b, 5.0)])
    star = SimpleNamespace(
        emitted_constraints=[_constraint(b, a, 5.5), _constraint(a, c, 7.0)]
    )

    result = compare_reports(nef, star)

    assert result["common_atom_pairs"] == 1
    assert result["star_only_atom_pairs"] == 1
    assert result["common_pairs_with_different_bounds"] == 1
    assert result["maximum_common_bound_delta_angstrom"] == 0.5
    assert result["exact_pair_and_bound_parity"] is False


def _endpoint(expression: str, hint: str | None = None) -> dict:
    return {
        "chain_code": "A",
        "sequence_code": "1",
        "residue_name": "ALA",
        "atom_expression": expression,
        "canonical_chain_code": "1" if hint else None,
        "canonical_sequence_code": "1" if hint else None,
        "canonical_residue_name": "ALA" if hint else None,
        "canonical_atom_hint": hint,
    }


def _alternative(row: str, atom1: str, atom2: str, hint2: str | None = None) -> dict:
    return {
        "endpoint1": _endpoint(atom1, atom1 if hint2 else None),
        "endpoint2": _endpoint(atom2, hint2),
        "upper_bound": 4.0,
        "row_ids": [row],
    }


def _report(constraint: dict | None, group: dict, rejections: list[dict] | None = None) -> dict:
    return {
        "sequence_map": [
            {
                "source_chain": "A",
                "source_sequence_code": "1",
                "residue_name": "ALA",
                "boltz_chain": "A",
                "boltz_residue_index": 1,
                "aliases": [["1", "1", "ALA"], ["A", "1", "ALA"]],
            }
        ],
        "emitted_constraints": [] if constraint is None else [constraint],
        "source_restraint_groups": [group],
        "rejections": rejections or [],
    }


def test_discrepancy_audit_preserves_atomset_vs_explicit_or_evidence():
    atom1 = {"chain": "A", "residue_index": 1, "atom_name": "CA"}
    atom2 = {"chain": "A", "residue_index": 1, "atom_name": "CB"}
    nef_group = {
        "list_name": "nef_distance_restraint_list_test",
        "restraint_id": "1",
        "alternatives": [_alternative("n1", "HA", "HB%")],
    }
    star_group = {
        "list_name": "test",
        "restraint_id": "1",
        "alternatives": [
            _alternative(str(index), "HA", "HB", f"HB{index}")
            for index in (1, 2, 3)
        ],
    }
    nef_constraint = {
        "atom1": atom1,
        "atom2": atom2,
        "max_distance": 7.04,
        "source_groups": ["nef_distance_restraint_list_test:1"],
        "provenance": [
            {
                "explicit_pair_count": 3,
                "averaging_policy": "sum-r6",
                "projected_upper_bound": 7.04,
            }
        ],
    }
    star_constraint = {
        "atom1": atom1,
        "atom2": atom2,
        "max_distance": 6.24,
        "source_groups": ["test:1"],
        "provenance": [
            {
                "explicit_pair_count": 1,
                "averaging_policy": "sum-r6",
                "projected_upper_bound": 6.24,
            }
        ],
    }

    rows = audit_reports(
        "TEST",
        _report(nef_constraint, nef_group),
        _report(star_constraint, star_group),
    )

    assert len(rows) == 1
    assert rows[0]["classification"] == "expected_format_difference"
    assert rows[0]["nef_pair_counts_N"] == [3]
    assert rows[0]["star_pair_counts_N"] == [1]
    assert "HB1" in " ".join(rows[0]["star_canonical_expansions"])
    assert "physical_atom_set" in rows[0]["nef_pseudoatom_handling"]


def test_fail_closed_gate_rejects_missing_coordinate_resolution(tmp_path):
    run = {
        "cases": [],
        "aggregate": {
            "formats": {
                "nef": {
                    "implication_failures": 0,
                    "constraints_missing_in_any_model": 1,
                },
                "star": {
                    "implication_failures": 0,
                    "constraints_missing_in_any_model": 0,
                },
            }
        },
        "discrepancy_audit": audit_summary([]),
    }
    baseline = tmp_path / "reviewed.json"
    baseline.write_text(__import__("json").dumps(_baseline_payload(run)), encoding="utf-8")

    gate = _evaluate_gate(run, baseline)

    assert gate["status"] == "fail"
    assert {failure["reason"] for failure in gate["failures"]} == {
        "missing_coordinate_resolution"
    }


def test_benchmark_artifact_paths_are_portable():
    assert _portable_path(Path("benchmark") / "input") == "benchmark/input"
