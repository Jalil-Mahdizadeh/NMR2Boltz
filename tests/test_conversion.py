from pathlib import Path

import pytest

from nmr2boltz.output import write_outputs
from nmr2boltz.project import ProjectionSettings, project_document
from nmr2boltz.star import parse_star_document
from nmr2boltz.topology import TopologyLibrary

FIXTURES = Path(__file__).parent / "fixtures"


def convert(name: str, **settings):
    path = FIXTURES / name
    parsed = parse_star_document(path, missing_upper_policy=settings.pop("missing_upper_policy", "reject"))
    return project_document(
        parsed,
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(**settings),
    )


def test_nef_projection_and_ambiguity():
    report = convert("example.nef")
    assert report.detected_format == "nef"
    assert [record.source_sequence_code for record in report.sequence_map] == ["10", "10A", "12", "13", "14"]
    emitted = {(c.atom1.atom_name, c.atom2.atom_name): c for c in report.emitted_constraints}
    assert ("CG1", "CB") in emitted
    expected = (9 ** (1 / 6)) * 4.0 + 1.12 + 1.12
    assert emitted[("CG1", "CB")].max_distance == pytest.approx(expected)
    assert any(group.restraint_id == "2" and len(group.alternatives) == 2 for group in report.ambiguous_groups)
    reasons = {item.reason for item in report.rejections}
    assert "missing_upper_bound" in reasons
    assert "same_heavy_parent_atom" in reasons
    assert "projected_bound_exceeds_boltz_maximum" in reasons


def test_missing_upper_can_be_derived_only_by_explicit_policy():
    path = FIXTURES / "example.nef"
    parsed = parse_star_document(path, missing_upper_policy="target-plus-uncertainty")
    report = project_document(
        parsed,
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(averaging_policy="hard-or"),
    )
    derived = next(
        provenance
        for constraint in report.emitted_constraints
        for provenance in constraint.provenance
        if provenance["group_id"].endswith(":4")
    )
    assert derived["source_observation"]["bound_source"] == "target_plus_uncertainty"
    assert derived["source_observation"]["target_value"] == 3.5
    assert derived["source_observation"]["target_uncertainty"] == 0.5
    assert not any(item.group_id.endswith(":4") and item.reason == "missing_upper_bound" for item in report.rejections)


def test_nmrstar_author_wildcard_rows_are_deduplicated():
    path = FIXTURES / "example.str"
    parsed = parse_star_document(path)
    assert parsed.detected_format == "nmr-star"
    assert parsed.sequence_resolver.records[0].aliases == [
        ("1", "1", "VAL"),
        ("1", "10", "VAL"),
        ("A", "1", "VAL"),
        ("A", "10", "VAL"),
    ]
    group1 = next(group for group in parsed.restraint_groups if group.restraint_id == "1")
    assert len(group1.alternatives) == 1
    assert group1.alternatives[0].endpoint1.atom_expression == "HG1%"
    assert set(group1.alternatives[0].row_ids) == {"1", "2", "3"}
    report = project_document(
        parsed,
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )
    assert len(report.emitted_constraints) == 1
    assert len(report.ambiguous_groups) == 1


def test_nmrstar_legacy_xplor_wildcard_rows_are_deduplicated(tmp_path):
    source = (FIXTURES / "example.str").read_text(encoding="utf-8")
    path = tmp_path / "legacy-wildcard.str"
    path.write_text(source.replace("HG1%", "HG1#").replace("HB%", "HB#"), encoding="utf-8")

    parsed = parse_star_document(path)
    group1 = next(group for group in parsed.restraint_groups if group.restraint_id == "1")

    assert len(group1.alternatives) == 1
    assert group1.alternatives[0].endpoint1.atom_expression == "HG1#"
    assert set(group1.alternatives[0].row_ids) == {"1", "2", "3"}


def test_embedded_custom_component_topology():
    path = FIXTURES / "custom_component.nef"
    parsed = parse_star_document(path)
    assert any(topology.comp_id == "LIG" for topology in parsed.embedded_topologies)
    report = project_document(
        parsed,
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )
    assert len(report.emitted_constraints) == 1
    atoms = {report.emitted_constraints[0].atom1.atom_name, report.emitted_constraints[0].atom2.atom_name}
    assert atoms == {"C1", "CA"}


def test_output_bundle(tmp_path):
    report = convert("example.nef")
    paths = write_outputs(report, tmp_path, hypothesis_count=3, random_seed=7)
    assert (tmp_path / "boltz_constraints.yaml").is_file()
    assert (tmp_path / "conversion_report.json").is_file()
    assert (tmp_path / "proposed_atom_contact_unions.yaml").is_file()
    assert (tmp_path / "sequences.fasta").read_text(encoding="utf-8") == ">A\nVAYLG\n"
    assert (tmp_path / "hypotheses" / "manifest.json").is_file()
    assert len(paths) >= 10


def test_sequence_only_nef_produces_empty_conversion(tmp_path):
    source = (FIXTURES / "example.nef").read_text(encoding="utf-8")
    path = tmp_path / "sequence-only.nef"
    path.write_text(source.split("save_distance_restraints", 1)[0], encoding="utf-8")

    parsed = parse_star_document(path)
    report = project_document(
        parsed,
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )

    assert report.detected_format == "nef"
    assert report.emitted_constraints == []
    assert len(report.sequence_map) == 5
    assert any("empty distance-constraint conversion" in warning for warning in report.warnings)


def test_residue_name_mismatch_is_rejected_before_topology(tmp_path):
    source = (FIXTURES / "example.nef").read_text(encoding="utf-8")
    path = tmp_path / "mismatched.nef"
    path.write_text(
        source.replace("A 10 VAL HG1% A 10A", "A 10 ALA HG1% A 10A", 1),
        encoding="utf-8",
    )

    report = project_document(
        parse_star_document(path),
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )

    mismatch = next(item for item in report.rejections if item.reason == "sequence_residue_mismatch")
    assert "restraint endpoint A:10:ALA:HG1%" in mismatch.details
    assert "A:1 (VAL)" in mismatch.details
