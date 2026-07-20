from __future__ import annotations

from itertools import product
from pathlib import Path

import pytest
import yaml

from nmr2boltz.model import Endpoint, RawAlternative, RestraintGroup, SequenceRecord
from nmr2boltz.output import write_outputs
from nmr2boltz.project import ProjectionSettings, project_document
from nmr2boltz.star import (
    ParsedStarDocument,
    SequenceResolver,
    normalize_nmrstar_canonical_expansions,
    parse_star_document,
)
from nmr2boltz.topology import TopologyLibrary


CORPUS = Path(__file__).parents[1] / "benchmark" / "input"
UPPER = 4.32


def _endpoint(sequence: str, component: str, author_atom: str, canonical: str) -> Endpoint:
    return Endpoint(
        chain_code="A",
        sequence_code=sequence,
        residue_name=component,
        atom_expression=author_atom,
        canonical_chain_code="A",
        canonical_sequence_code=sequence,
        canonical_residue_name=component,
        canonical_atom_hint=canonical,
    )


def _star_group(
    pairs: list[tuple[str, str, str, str]],
    *,
    component1: str,
    component2: str,
    restraint_id: str = "1",
) -> RestraintGroup:
    alternatives = []
    for index, (author1, canonical1, author2, canonical2) in enumerate(pairs, start=1):
        alternatives.append(
            RawAlternative(
                source_format="nmr-star",
                list_name="synthetic",
                restraint_id=restraint_id,
                endpoint1=_endpoint("1", component1, author1, canonical1),
                endpoint2=_endpoint("2", component2, author2, canonical2),
                upper_bound=UPPER,
                weight=1.0,
                member_id=str(index),
                member_logic_code="OR",
                row_ids=[str(index)],
            )
        )
    return RestraintGroup(
        source_format="nmr-star",
        list_name="synthetic",
        restraint_id=restraint_id,
        alternatives=alternatives,
        origin="NOE",
    )


def _nef_group(
    atom1: str,
    atom2: str,
    *,
    component1: str,
    component2: str,
) -> RestraintGroup:
    return RestraintGroup(
        source_format="nef",
        list_name="synthetic",
        restraint_id="1",
        alternatives=[
            RawAlternative(
                source_format="nef",
                list_name="synthetic",
                restraint_id="1",
                endpoint1=Endpoint("A", "1", component1, atom1),
                endpoint2=Endpoint("A", "2", component2, atom2),
                upper_bound=UPPER,
                weight=1.0,
                row_ids=["1"],
            )
        ],
        origin="NOE",
    )


def _report(group: RestraintGroup, component1: str, component2: str):
    resolver = SequenceResolver()
    resolver.add(SequenceRecord("A", "1", component1, "A", 1, "synthetic"))
    resolver.add(SequenceRecord("A", "2", component2, "A", 2, "synthetic"))
    parsed = ParsedStarDocument(
        entry=None,
        detected_format=group.source_format,
        sequence_resolver=resolver,
        restraint_groups=[group],
        embedded_topologies=[],
        warnings=[],
    )
    return project_document(
        parsed,
        input_file="synthetic",
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(averaging_policy="sum-r6"),
    )


def _cartesian_rows(
    author1: str,
    atoms1: list[str],
    author2: str,
    atoms2: list[str],
) -> list[tuple[str, str, str, str]]:
    return [
        (author1, atom1, author2, atom2)
        for atom1, atom2 in product(atoms1, atoms2)
    ]


@pytest.mark.parametrize(
    ("component1", "author1", "atoms1", "component2", "author2", "atoms2", "count"),
    [
        ("ALA", "HB", ["HB1", "HB2", "HB3"], "ALA", "HA", ["HA"], 3),
        ("ALA", "HB", ["HB1", "HB2", "HB3"], "ALA", "HB", ["HB1", "HB2", "HB3"], 9),
        ("GLY", "HA", ["HA2", "HA3"], "ALA", "HB", ["HB1", "HB2", "HB3"], 6),
    ],
)
def test_topology_verified_expansions_use_full_pair_multiplicity(
    component1, author1, atoms1, component2, author2, atoms2, count
):
    group = _star_group(
        _cartesian_rows(author1, atoms1, author2, atoms2),
        component1=component1,
        component2=component2,
    )
    report = _report(group, component1, component2)

    assert len(report.emitted_constraints) == 1
    provenance = report.emitted_constraints[0].provenance[0]
    assert provenance["explicit_pair_count"] == count
    assert provenance["averaging_factor"] == pytest.approx(count ** (1 / 6))
    assert provenance["source_rows"] == [str(index) for index in range(1, count + 1)]
    assert len(provenance["source_observation"]["canonical_expansions"]) == count


def test_incomplete_or_expansion_is_rejected_fail_closed():
    group = _star_group(
        _cartesian_rows("HB", ["HB1", "HB2"], "HA", ["HA"]),
        component1="ALA",
        component2="ALA",
    )
    report = _report(group, "ALA", "ALA")

    assert not report.emitted_constraints
    rejection = next(
        item
        for item in report.rejections
        if item.reason == "inconsistent_nmrstar_canonical_expansion"
    )
    assert "incomplete" in rejection.details
    assert rejection.row_ids == ["1", "2"]


def test_incomplete_two_sided_cartesian_expansion_is_rejected():
    rows = [
        ("HB", atom1, "HB", atom2)
        for atom1, atom2 in zip(
            ["HB1", "HB2", "HB3"], ["HB1", "HB2", "HB3"]
        )
    ]
    report = _report(
        _star_group(rows, component1="ALA", component2="ALA"), "ALA", "ALA"
    )

    assert not report.emitted_constraints
    assert any(
        item.reason == "inconsistent_nmrstar_canonical_expansion"
        and "Cartesian product" in item.details
        for item in report.rejections
    )


def test_inconsistent_expansion_bounds_are_rejected_not_split_into_n1_rows():
    group = _star_group(
        _cartesian_rows("HB", ["HB1", "HB2", "HB3"], "HA", ["HA"]),
        component1="ALA",
        component2="ALA",
    )
    group.alternatives[1].upper_bound = UPPER + 0.1
    report = _report(group, "ALA", "ALA")

    assert not report.emitted_constraints
    rejection = next(
        item
        for item in report.rejections
        if item.reason == "inconsistent_nmrstar_canonical_expansion"
    )
    assert "inconsistent bounds" in rejection.details
    assert rejection.row_ids == ["1", "2", "3"]


def test_unavailable_component_topology_rejects_expansion_in_full():
    report = _report(
        _star_group(
            _cartesian_rows("HX", ["HX1", "HX2"], "HA", ["HA"]),
            component1="ZZZ",
            component2="ALA",
        ),
        "ZZZ",
        "ALA",
    )

    assert not report.emitted_constraints
    rejection = next(
        item
        for item in report.rejections
        if item.reason == "inconsistent_nmrstar_canonical_expansion"
    )
    assert "topology-verified" in rejection.details
    assert rejection.row_ids == ["1", "2"]


def test_alias_equivalent_canonical_spellings_are_not_reconstructed_as_atom_set():
    group = _star_group(
        [
            ("HN", "HN", "HN", "HN"),
            ("HN", "H", "HN", "H"),
        ],
        component1="ALA",
        component2="ALA",
    )

    issues = normalize_nmrstar_canonical_expansions(group, TopologyLibrary())

    assert issues == []
    assert [alternative.row_ids for alternative in group.alternatives] == [["1"], ["2"]]
    assert [
        (
            alternative.endpoint1.canonical_atom_hint,
            alternative.endpoint2.canonical_atom_hint,
        )
        for alternative in group.alternatives
    ] == [("HN", "HN"), ("H", "H")]

    report = _report(group, "ALA", "ALA")
    assert len(report.emitted_constraints) == 1
    assert not report.ambiguous_groups
    assert not report.rejections
    assert report.emitted_constraints[0].provenance[0]["source_rows"] == ["1", "2"]


def test_leu_hd1_hd2_parent_ambiguity_remains_union_with_n3_per_branch():
    rows = _cartesian_rows("HD1", ["HD11", "HD12", "HD13"], "HA", ["HA"])
    rows += _cartesian_rows("HD2", ["HD21", "HD22", "HD23"], "HA", ["HA"])
    report = _report(
        _star_group(rows, component1="LEU", component2="ALA"), "LEU", "ALA"
    )

    assert not report.emitted_constraints
    assert len(report.ambiguous_groups) == 1
    alternatives = report.ambiguous_groups[0].alternatives
    assert {
        frozenset((item.atom1.atom_name, item.atom2.atom_name))
        for item in alternatives
    } == {frozenset(("CA", "CD1")), frozenset(("CA", "CD2"))}
    assert {item.explicit_pair_count for item in alternatives} == {3}


def test_equivalent_nef_and_star_exact_and_union_yaml_match(tmp_path):
    exact_nef = _report(_nef_group("HB%", "HA", component1="ALA", component2="ALA"), "ALA", "ALA")
    exact_star = _report(
        _star_group(
            _cartesian_rows("HB", ["HB1", "HB2", "HB3"], "HA", ["HA"]),
            component1="ALA",
            component2="ALA",
        ),
        "ALA",
        "ALA",
    )
    union_nef = _report(_nef_group("HDx%", "HA", component1="LEU", component2="ALA"), "LEU", "ALA")
    union_rows = _cartesian_rows("HD1", ["HD11", "HD12", "HD13"], "HA", ["HA"])
    union_rows += _cartesian_rows("HD2", ["HD21", "HD22", "HD23"], "HA", ["HA"])
    union_star = _report(
        _star_group(union_rows, component1="LEU", component2="ALA"),
        "LEU",
        "ALA",
    )

    reports = {
        "exact-nef": exact_nef,
        "exact-star": exact_star,
        "union-nef": union_nef,
        "union-star": union_star,
    }
    for name, report in reports.items():
        write_outputs(report, tmp_path / name)

    for filename in ("atom_constraints_exact.yaml", "atom_constraints_union.yaml"):
        assert yaml.safe_load((tmp_path / "exact-nef" / filename).read_text()) == yaml.safe_load(
            (tmp_path / "exact-star" / filename).read_text()
        )
        assert yaml.safe_load((tmp_path / "union-nef" / filename).read_text()) == yaml.safe_load(
            (tmp_path / "union-star" / filename).read_text()
        )


def test_43jx_restraint_468_star_expansion_matches_nef_bound():
    nef_path = CORPUS / "43JX" / "43jx_nmr-data.nef"
    star_path = CORPUS / "43JX" / "43jx_nmr-data.str"
    reports = {}
    groups = {}
    for label, path in (("nef", nef_path), ("star", star_path)):
        parsed = parse_star_document(path)
        report = project_document(
            parsed,
            input_file=str(path),
            topology_library=TopologyLibrary(),
            settings=ProjectionSettings(),
        )
        reports[label] = report
        groups[label] = next(
            group
            for group in parsed.restraint_groups
            if group.restraint_id == "468" and "CYANA_distance_restraints_2" in group.list_name
        )

    projected = {}
    for label, report in reports.items():
        group_id = groups[label].group_id
        projected[label] = next(
            provenance["projected_upper_bound"]
            for constraint in report.emitted_constraints
            for provenance in constraint.provenance
            if provenance["group_id"] == group_id
        )
    assert projected["star"] == pytest.approx(projected["nef"])
    assert projected["star"] == pytest.approx(7.428047646360333)
    star_alternative = groups["star"].alternatives[0]
    assert star_alternative.row_ids == ["796", "797", "798"]
    assert star_alternative.endpoint2.canonical_atom_set == ["HB1", "HB2", "HB3"]
