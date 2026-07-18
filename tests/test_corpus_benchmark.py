import copy
import json
from pathlib import Path
from types import SimpleNamespace

from nmr2boltz.model import BoltzAtom, EmittedConstraint
from validation.benchmark_corpus import (
    _baseline_payload,
    _evaluate_gate,
    _missing_coordinate_review,
    _portable_path,
    compare_reports,
)
from validation.discrepancy_audit import (
    _rejected_geometric_pseudoatom,
    _stereo_assignment_vs_physical_set,
    _strictly_equivalent,
    _verified_canonical_alias,
    _wildcard_set_vs_explicit_or,
    audit_reports,
    audit_summary,
)


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


def _semantic(
    kinds: set[str],
    alternatives: list[set[str]],
    *,
    raw: str,
    rejections: set[str] | None = None,
) -> dict:
    pairs = frozenset(pair for alternative in alternatives for pair in alternative)
    return {
        "alternatives": tuple({"physical_pairs": tuple(sorted(item))} for item in alternatives),
        "raw_signature": (raw,),
        "physical_alternative_signature": tuple(
            tuple(sorted(item)) for item in alternatives
        ),
        "physical_pairs": pairs,
        "kinds": frozenset(kinds),
        "rejection_reasons": frozenset(rejections or set()),
        "unverified_heavy_atoms": frozenset(),
    }


def test_wildcard_vs_explicit_or_predicate_is_structural_and_fail_closed():
    wildcard = _semantic(
        {"physical_atom_set", "explicit_atom_or_alias"},
        [{"p1", "p2", "p3"}],
        raw="HB%",
    )
    members = _semantic(
        {"explicit_atom_or_alias"}, [{"p1"}, {"p2"}], raw="HB1|HB2"
    )
    disjoint = _semantic({"explicit_atom_or_alias"}, [{"p4"}], raw="HG1")
    rejected_members = _semantic(
        {"explicit_atom_or_alias"},
        [{"p1"}, {"p2"}],
        raw="HB1|HB2",
        rejections={"missing_upper_bound"},
    )

    assert _wildcard_set_vs_explicit_or(wildcard, members)
    assert not _wildcard_set_vs_explicit_or(wildcard, disjoint)
    assert not _wildcard_set_vs_explicit_or(wildcard, rejected_members)


def test_stereo_assignment_vs_physical_set_predicate_requires_compatible_atoms():
    stereo = _semantic(
        {"stereospecific_assignment_alternative", "explicit_atom_or_alias"},
        [{"p1", "p2"}],
        raw="HBx",
    )
    physical = _semantic({"physical_atom_set"}, [{"p2"}], raw="HB2")
    disjoint = _semantic({"physical_atom_set"}, [{"p3"}], raw="HG%")
    same_parent_branch = _semantic(
        {"stereospecific_assignment_alternative", "physical_atom_set"},
        [{"p1", "p2"}],
        raw="HDx%--HDy%",
        rejections={"partially_unprojectable_or_group", "same_heavy_parent_atom"},
    )
    arbitrary_rejection = _semantic(
        {"stereospecific_assignment_alternative"},
        [{"p1", "p2"}],
        raw="HBx",
        rejections={"missing_upper_bound"},
    )

    assert _stereo_assignment_vs_physical_set(stereo, physical)
    assert _stereo_assignment_vs_physical_set(same_parent_branch, physical)
    assert not _stereo_assignment_vs_physical_set(stereo, disjoint)
    assert not _stereo_assignment_vs_physical_set(arbitrary_rejection, physical)


def test_geometric_pseudoatom_predicate_requires_exact_rejection_policy():
    pseudoatom = _semantic(
        {"rejected_geometric_pseudoatom"},
        [],
        raw="QB",
        rejections={"unresolved_atom_topology"},
    )
    physical = _semantic({"physical_atom_set"}, [{"p1", "p2"}], raw="HB%")
    arbitrary_rejection = _semantic(
        {"rejected_geometric_pseudoatom"},
        [],
        raw="QB",
        rejections={"missing_upper_bound"},
    )

    assert _rejected_geometric_pseudoatom(pseudoatom, physical)
    assert not _rejected_geometric_pseudoatom(arbitrary_rejection, physical)


def test_verified_alias_predicate_requires_identical_physical_or_structure():
    hn = _semantic({"explicit_atom_or_alias"}, [{"amide-H"}], raw="HN")
    h = _semantic({"explicit_atom_or_alias"}, [{"amide-H"}], raw="H")
    different = _semantic({"explicit_atom_or_alias"}, [{"HA"}], raw="HA")

    assert _verified_canonical_alias(hn, h)
    assert not _verified_canonical_alias(hn, different)


def test_semantic_equivalence_accepts_topology_proven_canonical_reconstruction():
    nef_set = _semantic(
        {"explicit_atom_or_alias", "physical_atom_set"},
        [{"HA--HB2", "HA--HB3"}],
        raw="HA--HB%",
    )
    star_set = _semantic(
        {"explicit_atom_or_alias", "physical_atom_set"},
        [{"HA--HB2", "HA--HB3"}],
        raw="HA--HB[HB2,HB3]",
    )

    assert _strictly_equivalent(nef_set, star_set)


def test_semantic_equivalence_rejects_superficially_matching_pair_sets():
    physical_set = _semantic(
        {"explicit_atom_or_alias", "physical_atom_set"},
        [{"p1", "p2"}],
        raw="HB%",
    )
    independent_or = _semantic(
        {"explicit_atom_or_alias", "physical_atom_set"},
        [{"p1"}, {"p2"}],
        raw="HB1|HB2",
    )
    wrong_semantics = _semantic(
        {"stereospecific_assignment_alternative"},
        [{"p1", "p2"}],
        raw="HBx",
    )
    rejected = _semantic(
        {"explicit_atom_or_alias", "physical_atom_set"},
        [{"p1", "p2"}],
        raw="HB",
        rejections={"missing_upper_bound"},
    )

    assert not _strictly_equivalent(physical_set, independent_or)
    assert not _strictly_equivalent(physical_set, wrong_semantics)
    assert not _strictly_equivalent(physical_set, rejected)


def test_equivalent_source_evidence_with_different_bound_is_parser_bug():
    atom1 = {"chain": "A", "residue_index": 1, "atom_name": "CA"}
    atom2 = {"chain": "A", "residue_index": 1, "atom_name": "CB"}
    group = {
        "list_name": "test",
        "restraint_id": "1",
        "alternatives": [_alternative("1", "HA", "HB1")],
    }
    first = {
        "atom1": atom1,
        "atom2": atom2,
        "max_distance": 5.0,
        "source_groups": ["test:1"],
        "provenance": [],
    }
    second = {**first, "max_distance": 5.5}

    rows = audit_reports("TEST", _report(first, group), _report(second, group))

    assert rows[0]["classification"] == "parser_projection_bug"


def test_arbitrary_rejection_is_not_an_expected_format_difference():
    atom1 = {"chain": "A", "residue_index": 1, "atom_name": "CA"}
    atom2 = {"chain": "A", "residue_index": 1, "atom_name": "CB"}
    group = {
        "list_name": "test",
        "restraint_id": "1",
        "alternatives": [_alternative("1", "HA", "HB1")],
    }
    constraint = {
        "atom1": atom1,
        "atom2": atom2,
        "max_distance": 5.0,
        "source_groups": ["test:1"],
        "provenance": [],
    }
    rejected = [{"group_id": "test:1", "reason": "missing_upper_bound"}]

    rows = audit_reports(
        "TEST",
        _report(constraint, group, rejected),
        _report({**constraint, "max_distance": 5.5}, group),
    )

    assert rows[0]["classification"] == "unresolved"


def _gate_run(missing_contacts: list[dict]) -> dict:
    def conversion(contacts: list[dict]) -> dict:
        return {
            "restraint_groups": 1,
            "source_alternatives": 1,
            "emitted_constraints": 1,
            "ambiguous_groups": 0,
            "rejection_reasons": {},
            "atom_topology_validation": {
                "checked_constraints": 1,
                "violation_count": 0,
            },
            "coordinates": {
                "heavy_atom_constraints": {
                    "resolved_model_cases": 0,
                    "satisfied_model_cases": 0,
                    "violated_model_cases": 0,
                    "missing_model_cases": len(contacts),
                    "constraints_missing_in_any_model": len(contacts),
                    "missing_coordinate_contacts": contacts,
                },
                "projection_implication_audit": {
                    "antecedent_cases": 0,
                    "failure_count": 0,
                },
            },
        }

    return {
        "cases": [
            {
                "case_id": "TEST",
                "status": "pass",
                "format_parity": {},
                "formats": {
                    "nef": conversion(missing_contacts),
                    "star": conversion([]),
                },
            }
        ],
        "aggregate": {
            "formats": {
                "nef": {
                    "implication_failures": 0,
                },
                "star": {
                    "implication_failures": 0,
                },
            }
        },
        "discrepancy_audit": audit_summary([]),
    }


def test_fail_closed_gate_accepts_only_exact_reviewed_coordinate_set(tmp_path):
    contact = {
        "constraint_id": "C00001",
        "atom1": "A:1:CA",
        "atom2": "A:2:CB",
        "upper_bound_angstrom": 6.0,
        "models_resolved": 0,
        "source_groups": "test:1",
    }
    run = _gate_run([contact])
    baseline = tmp_path / "reviewed.json"
    baseline.write_text(json.dumps(_baseline_payload(run)), encoding="utf-8")

    gate = _evaluate_gate(run, baseline)
    changed_bound = copy.deepcopy(run)
    changed_bound["cases"][0]["formats"]["nef"]["coordinates"][
        "heavy_atom_constraints"
    ]["missing_coordinate_contacts"][0]["upper_bound_angstrom"] = 6.1
    removed = _gate_run([])
    added = _gate_run(
        [
            contact,
            {
                **contact,
                "constraint_id": "C00002",
                "atom2": "A:3:CG",
            },
        ]
    )

    assert gate["status"] == "pass"
    assert _missing_coordinate_review(run)["contact_count"] == 1
    for changed_run in (changed_bound, removed, added):
        changed_gate = _evaluate_gate(changed_run, baseline)
        assert changed_gate["status"] == "fail"
        assert "changed_reviewed_missing_coordinate_set" in {
            failure["reason"] for failure in changed_gate["failures"]
        }


def test_gate_remains_fail_closed_for_every_scientific_failure(tmp_path):
    run = _gate_run([])
    baseline = tmp_path / "reviewed.json"
    baseline.write_text(json.dumps(_baseline_payload(run)), encoding="utf-8")

    implication = copy.deepcopy(run)
    implication["aggregate"]["formats"]["nef"]["implication_failures"] = 1

    unresolved = copy.deepcopy(run)
    unresolved["discrepancy_audit"] = audit_summary(
        [
            {
                "audit_id": "TEST|nef_only|A:1:CA--A:2:CB",
                "case_id": "TEST",
                "discrepancy_type": "nef_only",
                "classification": "unresolved",
            }
        ]
    )

    parser_bug = copy.deepcopy(run)
    parser_bug["discrepancy_audit"] = audit_summary(
        [
            {
                "audit_id": "TEST|different_bound|A:1:CA--A:2:CB",
                "case_id": "TEST",
                "discrepancy_type": "different_bound",
                "classification": "parser_projection_bug",
            }
        ]
    )

    changed_audit = copy.deepcopy(run)
    changed_audit["discrepancy_audit"]["digest_sha256"] = "changed"

    changed_metric = copy.deepcopy(run)
    changed_metric["cases"][0]["formats"]["nef"]["emitted_constraints"] = 2

    topology_violation = copy.deepcopy(run)
    topology_violation["cases"][0]["formats"]["nef"]["atom_topology_validation"][
        "violation_count"
    ] = 1

    expected = (
        (implication, "projection_implication_failure"),
        (unresolved, "unresolved_format_discrepancy"),
        (parser_bug, "parser_projection_bug"),
        (changed_audit, "unreviewed_discrepancy_change"),
        (changed_metric, "unreviewed_metric_change"),
        (topology_violation, "emitted_atom_topology_violation"),
    )
    for changed_run, reason in expected:
        gate = _evaluate_gate(changed_run, baseline)
        assert gate["status"] == "fail"
        assert reason in {failure["reason"] for failure in gate["failures"]}


def test_benchmark_artifact_paths_are_portable():
    assert _portable_path(Path("benchmark") / "input") == "benchmark/input"
