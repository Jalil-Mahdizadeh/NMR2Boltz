import copy
import json
import math
from pathlib import Path

import pytest
import yaml

import nmr2boltz.output as output_module
from nmr2boltz.model import Endpoint, SequenceRecord
from nmr2boltz.output import _rounded, write_outputs
from nmr2boltz.project import ProjectionSettings, project_document
from nmr2boltz.star import (
    SequenceResolver,
    StarDataError,
    extract_restraint_groups,
    parse_star_document,
)
from nmr2boltz.topology import TopologyLibrary, TopologyResolutionError


FIXTURES = Path(__file__).parent / "fixtures"


class _StarLoop:
    category = "_Gen_dist_constraint"
    tags = [
        "_Gen_dist_constraint.ID",
        "_Gen_dist_constraint.Member_ID",
        "_Gen_dist_constraint.Member_logic_code",
        "_Gen_dist_constraint.Auth_asym_ID_1",
        "_Gen_dist_constraint.Auth_seq_ID_1",
        "_Gen_dist_constraint.Auth_comp_ID_1",
        "_Gen_dist_constraint.Auth_atom_ID_1",
        "_Gen_dist_constraint.Auth_asym_ID_2",
        "_Gen_dist_constraint.Auth_seq_ID_2",
        "_Gen_dist_constraint.Auth_comp_ID_2",
        "_Gen_dist_constraint.Auth_atom_ID_2",
        "_Gen_dist_constraint.Distance_upper_bound_val",
    ]

    def __init__(self, logic_code):
        self.data = [["1", "1", logic_code, "A", "1", "ALA", "HA", "A", "2", "GLY", "HA2", "4.0"]]


class _Saveframe(list):
    name = "constraints"
    tags = [["_Gen_dist_constraint_list.Constraint_type", "NOE"]]


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("projection_margin", math.nan),
        ("projection_margin", math.inf),
        ("boltz_min_distance", math.nan),
        ("boltz_max_distance", math.inf),
    ],
)
def test_projection_settings_reject_non_finite_values(field, value):
    settings = ProjectionSettings(**{field: value})
    with pytest.raises(ValueError, match="finite"):
        settings.validate()


@pytest.mark.parametrize(("logic_code", "complex_logic"), [("OR", False), ("AND", True), ("XOR", True)])
def test_nmrstar_member_logic_is_only_flattened_for_explicit_or(logic_code, complex_logic):
    warnings = []
    groups = extract_restraint_groups(
        [_Saveframe([_StarLoop(logic_code)])],
        "nmr-star",
        missing_upper_policy="reject",
        origins=None,
        warnings=warnings,
    )
    assert groups[0].complex_logic is complex_logic
    if complex_logic:
        assert "only explicit OR is flattened" in groups[0].warnings[0]


def test_invalid_missing_upper_policy_is_rejected():
    with pytest.raises(StarDataError, match="Unsupported missing upper-bound policy"):
        parse_star_document(FIXTURES / "example.nef", missing_upper_policy="invented")


def test_negative_target_uncertainty_cannot_derive_upper_bound(tmp_path):
    source = (FIXTURES / "example.nef").read_text(encoding="utf-8")
    path = tmp_path / "negative-uncertainty.nef"
    path.write_text(source.replace("3.5 0.5 . . . .", "3.5 -0.5 . . . ."), encoding="utf-8")
    parsed = parse_star_document(path, missing_upper_policy="target-plus-uncertainty")
    alternative = next(group for group in parsed.restraint_groups if group.restraint_id == "4").alternatives[0]
    assert alternative.upper_bound is None
    assert alternative.bound_source == "missing"
    assert "must be non-negative" in alternative.warnings[0]


def test_yaml_rounding_is_outward_for_upper_bounds():
    value = 3.0000001
    rounded = _rounded(value)
    assert rounded == 3.000001
    assert rounded >= value


def test_partial_or_group_is_not_emitted_when_one_alternative_fails():
    parsed = parse_star_document(FIXTURES / "example.nef")
    group = next(item for item in parsed.restraint_groups if item.restraint_id == "1")
    invalid = copy.deepcopy(group.alternatives[0])
    invalid.endpoint1.atom_expression = "HZ9"
    invalid.row_ids = ["synthetic-invalid-alternative"]
    group.alternatives.append(invalid)

    report = project_document(
        parsed,
        input_file=str(FIXTURES / "example.nef"),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )

    assert not any(group.group_id in constraint.source_groups for constraint in report.emitted_constraints)
    assert any(
        item.group_id == group.group_id and item.reason == "partially_unprojectable_or_group"
        for item in report.rejections
    )


def test_atom_set_with_same_parent_branch_is_not_partially_emitted():
    parsed = parse_star_document(FIXTURES / "example.nef")
    group = next(item for item in parsed.restraint_groups if item.restraint_id == "1")
    alternative = group.alternatives[0]
    alternative.endpoint1.atom_expression = "HG%"
    alternative.endpoint2 = copy.deepcopy(alternative.endpoint1)
    alternative.endpoint2.atom_expression = "HG1%"
    alternative.row_ids = ["synthetic-mixed-parent-set"]
    group.alternatives = [alternative]
    parsed.restraint_groups = [group]

    report = project_document(
        parsed,
        input_file=str(FIXTURES / "example.nef"),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )

    assert not report.emitted_constraints
    assert any(
        item.reason == "atom_set_contains_same_heavy_parent_pair"
        for item in report.rejections
    )


def test_residue_map_replaces_records_and_requires_positive_indices(tmp_path):
    valid_map = tmp_path / "valid.tsv"
    valid_map.write_text(
        "source_chain\tsource_sequence_code\tsource_residue_name\tboltz_chain\tboltz_residue_index\n"
        "A\t10\tVAL\tZ\t101\n"
        "A\t10A\tALA\tZ\t102\n"
        "A\t12\tTYR\tZ\t103\n"
        "A\t13\tLEU\tZ\t104\n"
        "A\t14\tGLY\tZ\t105\n",
        encoding="utf-8",
    )
    parsed = parse_star_document(FIXTURES / "example.nef", residue_map_path=valid_map)
    assert len(parsed.sequence_resolver.records) == 5
    assert {record.boltz_chain for record in parsed.sequence_resolver.records} == {"Z"}

    invalid_map = tmp_path / "invalid.tsv"
    invalid_map.write_text(
        "source_chain\tsource_sequence_code\tsource_residue_name\tboltz_chain\tboltz_residue_index\n"
        "A\t10\tVAL\tZ\t0\n",
        encoding="utf-8",
    )
    with pytest.raises(StarDataError, match="positive"):
        parse_star_document(FIXTURES / "example.nef", residue_map_path=invalid_map)


def test_missing_chain_conflicting_author_and_canonical_identifiers_are_rejected():
    resolver = SequenceResolver()
    resolver.add(SequenceRecord("A", "1", "ALA", "A", 1, "synthetic"))
    resolver.add(SequenceRecord("B", "2", "ALA", "B", 1, "synthetic"))
    endpoint = Endpoint(
        chain_code=None,
        sequence_code="1",
        residue_name="ALA",
        atom_expression="HA",
        canonical_chain_code=None,
        canonical_sequence_code="2",
        canonical_residue_name="ALA",
        canonical_atom_hint="HA",
    )

    with pytest.raises(
        StarDataError,
        match=r"Conflicting missing-chain identifiers.*'1'->A:1, '2'->B:2",
    ):
        resolver.resolve(endpoint)


def test_missing_chain_equivalent_identifiers_resolve_once():
    resolver = SequenceResolver()
    record = SequenceRecord("A", "1", "ALA", "A", 1, "synthetic")
    resolver.add(record)
    endpoint = Endpoint(
        None, "1", "ALA", "HA", None, "1", "ALA", "HA"
    )

    resolved, warnings = resolver.resolve(endpoint)

    assert resolved is record
    assert warnings == [
        "missing/unresolved chain for ?:1:ALA:HA inferred as A"
    ]


def test_non_finite_bond_length_override_is_rejected(tmp_path):
    config = tmp_path / "bond-lengths.json"
    config.write_text(json.dumps({"element_upper": {"C": math.nan}}), encoding="utf-8")
    with pytest.raises(TopologyResolutionError, match="finite"):
        TopologyLibrary(bond_length_config=config)


def test_reused_output_directory_removes_stale_hypotheses(tmp_path):
    parsed = parse_star_document(FIXTURES / "example.nef")
    report = project_document(
        parsed,
        input_file=str(FIXTURES / "example.nef"),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )
    write_outputs(report, tmp_path, hypothesis_count=3, random_seed=7)
    write_outputs(report, tmp_path, hypothesis_count=1, random_seed=7)
    hypothesis_files = sorted((tmp_path / "hypotheses").glob("hypothesis_*.yaml"))
    assert [path.name for path in hypothesis_files] == ["hypothesis_0001.yaml"]
    assert yaml.safe_load(hypothesis_files[0].read_text(encoding="utf-8"))["constraints"]


def test_output_bundle_failure_leaves_prior_directory_unchanged(
    tmp_path, monkeypatch
):
    parsed = parse_star_document(FIXTURES / "example.nef")
    report = project_document(
        parsed,
        input_file=str(FIXTURES / "example.nef"),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )
    destination = tmp_path / "bundle"
    destination.mkdir()
    old_exact = destination / "atom_constraints_exact.yaml"
    old_exact.write_text("old exact bundle\n", encoding="utf-8")
    marker = destination / "keep.me"
    marker.write_text("prior output\n", encoding="utf-8")

    def fail_tsv(*_args, **_kwargs):
        raise OSError("simulated filesystem failure")

    monkeypatch.setattr(output_module, "_write_tsv", fail_tsv)

    with pytest.raises(OSError, match="simulated filesystem failure"):
        write_outputs(report, destination)

    assert old_exact.read_text(encoding="utf-8") == "old exact bundle\n"
    assert marker.read_text(encoding="utf-8") == "prior output\n"
    assert not (destination / "conversion_report.json").exists()
    assert not list(tmp_path.glob(".bundle.staging-*"))
    assert not list(tmp_path.glob(".bundle.backup-*"))


def test_output_bundle_commit_failure_restores_prior_directory(
    tmp_path, monkeypatch
):
    parsed = parse_star_document(FIXTURES / "example.nef")
    report = project_document(
        parsed,
        input_file=str(FIXTURES / "example.nef"),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )
    destination = tmp_path / "bundle"
    destination.mkdir()
    marker = destination / "prior.txt"
    marker.write_text("prior output\n", encoding="utf-8")
    original_replace = Path.replace

    def fail_staging_commit(path, target):
        if path.name.startswith(".bundle.staging-") and Path(target) == destination:
            raise OSError("simulated commit failure")
        return original_replace(path, target)

    monkeypatch.setattr(Path, "replace", fail_staging_commit)

    with pytest.raises(
        output_module.OutputBundleError,
        match="Unable to commit staged output bundle",
    ):
        write_outputs(report, destination)

    assert marker.read_text(encoding="utf-8") == "prior output\n"
    assert not (destination / "conversion_report.json").exists()
    assert not list(tmp_path.glob(".bundle.staging-*"))
    assert not list(tmp_path.glob(".bundle.backup-*"))
