from __future__ import annotations

import copy
import re
from pathlib import Path

import yaml
import pytest

from nmr2boltz.model import BoltzAtom
from nmr2boltz.output import OutputBoundValidationError, write_outputs
from nmr2boltz.project import ProjectionSettings, project_document
from nmr2boltz.star import parse_star_document
from nmr2boltz.topology import AtomTopologyValidationError, TopologyLibrary


FIXTURES = Path(__file__).parent / "fixtures"
FORBIDDEN_METADATA = {
    "group_id",
    "restraint_id",
    "source_rows",
    "x_nmr2boltz_schema",
    "compatible_with_current_boltzui_atom_contact_parser",
    "note",
    "label",
    "labels",
}


def _report():
    path = FIXTURES / "example.nef"
    return project_document(
        parse_star_document(path),
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )


def _all_keys(value):
    if isinstance(value, dict):
        for key, child in value.items():
            yield key
            yield from _all_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from _all_keys(child)


def test_exact_and_union_files_have_disjoint_minimal_schemas(tmp_path):
    report = _report()
    write_outputs(report, tmp_path)

    exact = yaml.safe_load((tmp_path / "atom_constraints_exact.yaml").read_text())
    union = yaml.safe_load((tmp_path / "atom_constraints_union.yaml").read_text())

    assert set(exact) == {"constraints"}
    assert exact["constraints"]
    assert all(set(item) == {"atom_contact"} for item in exact["constraints"])
    assert all(
        set(item["atom_contact"]) == {"atom1", "atom2", "max_distance", "force"}
        for item in exact["constraints"]
    )

    assert set(union) == {"constraints"}
    assert union["constraints"]
    assert all(set(item) == {"atom_contact_union"} for item in union["constraints"])
    for item in union["constraints"]:
        contact_union = item["atom_contact_union"]
        assert set(contact_union) == {"alternatives", "force"}
        assert len(contact_union["alternatives"]) >= 2
        assert all(
            set(alternative) == {"atom1", "atom2", "max_distance"}
            for alternative in contact_union["alternatives"]
        )
    assert FORBIDDEN_METADATA.isdisjoint(_all_keys(union))


def test_empty_constraint_files_are_explicit_empty_lists(tmp_path):
    report = _report()
    report.emitted_constraints = []
    report.ambiguous_groups = []

    write_outputs(report, tmp_path)

    assert (tmp_path / "atom_constraints_exact.yaml").read_text() == "constraints: []\n"
    assert (tmp_path / "atom_constraints_union.yaml").read_text() == "constraints: []\n"


def test_union_preserves_per_alternative_bounds_with_outward_six_decimal_format(tmp_path):
    report = _report()
    group = report.ambiguous_groups[0]
    assert len(group.alternatives) >= 2
    group.alternatives[0].max_distance = 6.4312691
    group.alternatives[1].max_distance = 7.0000001

    write_outputs(report, tmp_path)

    text = (tmp_path / "atom_constraints_union.yaml").read_text()
    payload = yaml.safe_load(text)
    bounds = sorted(
        alternative["max_distance"]
        for wrapped in payload["constraints"]
        for alternative in wrapped["atom_contact_union"]["alternatives"]
    )
    assert 6.431270 in bounds
    assert 7.000001 in bounds
    serialized = re.findall(r"max_distance: ([0-9]+\.[0-9]+)$", text, re.MULTILINE)
    assert serialized
    assert all(re.fullmatch(r"[0-9]+\.[0-9]{6}", value) for value in serialized)


def test_projected_union_bounds_are_raised_to_configured_minimum_before_serialization(
    tmp_path,
):
    path = FIXTURES / "example.nef"
    report = project_document(
        parse_star_document(path),
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(boltz_min_distance=10.0),
    )

    adjusted = [
        alternative
        for group in report.ambiguous_groups
        for alternative in group.alternatives
        if alternative.raw_projected_distance is not None
    ]
    assert adjusted
    assert all(alternative.max_distance == 10.0 for alternative in adjusted)
    assert all(alternative.raw_projected_distance < 10.0 for alternative in adjusted)

    write_outputs(report, tmp_path)
    payload = yaml.safe_load((tmp_path / "atom_constraints_union.yaml").read_text())
    assert all(
        alternative["max_distance"] >= 10.0
        for wrapped in payload["constraints"]
        for alternative in wrapped["atom_contact_union"]["alternatives"]
    )


def test_project_document_quarantines_complete_union_above_configured_maximum():
    path = FIXTURES / "example.nef"
    report = project_document(
        parse_star_document(path),
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(boltz_max_distance=7.5),
    )

    assert not report.ambiguous_groups
    rejection = next(
        item
        for item in report.rejections
        if item.reason == "ambiguous_group_bound_exceeds_boltz_maximum"
    )
    assert rejection.group_id == "distance_restraints:2"
    assert rejection.provenance["complete_union_group_quarantined"] is True
    assert report.statistics["union_groups_quarantined_above_boltz_maximum"] == 1


@pytest.mark.parametrize("bound", [1.999999, 20.000001, float("nan")])
def test_final_output_bound_invariant_rejects_corrupt_union_report(
    tmp_path, bound
):
    report = _report()
    report.ambiguous_groups[0].alternatives[0].max_distance = bound
    destination = tmp_path / "invalid-union-bound"

    with pytest.raises(OutputBoundValidationError, match="bounds must be finite"):
        write_outputs(report, destination)

    assert not destination.exists()


def test_exact_bounds_use_outward_six_decimal_format(tmp_path):
    report = _report()
    report.emitted_constraints[0].max_distance = 7.2
    write_outputs(report, tmp_path)

    text = (tmp_path / "atom_constraints_exact.yaml").read_text()
    assert "max_distance: 7.200000" in text
    assert yaml.safe_load(text)["constraints"][0]["atom_contact"]["max_distance"] == 7.2


def test_union_order_is_deterministic_and_old_outputs_are_removed(tmp_path):
    report = _report()
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    (first / "boltz_constraints.yaml").write_text("stale\n")
    (first / "proposed_atom_contact_unions.yaml").write_text("stale\n")
    write_outputs(report, first)

    reordered = copy.deepcopy(report)
    reordered.ambiguous_groups.reverse()
    for group in reordered.ambiguous_groups:
        group.alternatives.reverse()
    write_outputs(reordered, second)

    assert (first / "atom_constraints_union.yaml").read_bytes() == (
        second / "atom_constraints_union.yaml"
    ).read_bytes()
    assert not (first / "boltz_constraints.yaml").exists()
    assert not (first / "proposed_atom_contact_unions.yaml").exists()


def test_invalid_union_atom_fails_before_any_output_is_written(tmp_path):
    report = _report()
    report.ambiguous_groups[0].alternatives[0].atom1 = BoltzAtom("A", 1, "ZN")
    destination = tmp_path / "invalid-union"

    with pytest.raises(AtomTopologyValidationError, match="topology validation failed"):
        write_outputs(report, destination)

    assert not destination.exists()
