from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from nmr2boltz.cli import main
from nmr2boltz.model import BoltzAtom
from nmr2boltz.output import write_outputs
from nmr2boltz.project import ProjectionSettings, project_document
from nmr2boltz.star import parse_star_document
from nmr2boltz.topology import TopologyLibrary


FIXTURE = Path(__file__).parent / "fixtures" / "intraresidue_filter.nef"


def _convert(*, include_intraresidue: bool):
    return project_document(
        parse_star_document(FIXTURE),
        input_file=str(FIXTURE),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(
            include_intraresidue=include_intraresidue
        ),
    )


def _filtered_rejection(report, group_id: str):
    return next(
        rejection
        for rejection in report.rejections
        if rejection.group_id == group_id
        and rejection.reason == "intraresidue_filtered"
    )


def test_exclude_intraresidue_preserves_exact_contact_behavior_and_provenance():
    report = _convert(include_intraresidue=False)

    rejection = _filtered_rejection(report, "restraints:1")
    assert rejection.details == (
        "The projected restraint is intraresidue and was removed by policy."
    )
    assert rejection.row_ids == ["1"]
    assert rejection.provenance["filter"] == "exclude_intraresidue"
    assert rejection.provenance["residue_scope"] == "intraresidue"
    assert len(rejection.provenance["projected_alternatives"]) == 1
    assert not any(
        constraint.atom1.chain == constraint.atom2.chain
        and constraint.atom1.residue_index
        == constraint.atom2.residue_index
        for constraint in report.emitted_constraints
    )


def test_exclude_intraresidue_rejects_complete_all_intraresidue_union():
    report = _convert(include_intraresidue=False)

    rejection = _filtered_rejection(report, "restraints:2")
    assert rejection.row_ids == ["2", "3"]
    assert rejection.provenance["residue_scope"] == "intraresidue"
    assert {
        alternative["scope"]
        for alternative in rejection.provenance["projected_alternatives"]
    } == {"intraresidue"}
    assert len(rejection.provenance["projected_alternatives"]) == 2
    assert "every alternative is intraresidue" in rejection.details
    assert not any(
        group.group_id == "restraints:2"
        for group in report.ambiguous_groups
    )


def test_exclude_intraresidue_rejects_complete_mixed_scope_union():
    report = _convert(include_intraresidue=False)

    rejection = _filtered_rejection(report, "restraints:3")
    assert rejection.row_ids == ["4", "5"]
    assert (
        rejection.provenance["residue_scope"]
        == "mixed_intraresidue_interresidue"
    )
    assert {
        alternative["scope"]
        for alternative in rejection.provenance["projected_alternatives"]
    } == {"intraresidue", "interresidue"}
    assert "narrow and strengthen" in rejection.details
    assert report.statistics["intraresidue_groups_filtered"] == 3
    assert (
        report.statistics[
            "projected_alternatives_removed_by_intraresidue_filter"
        ]
        == 5
    )
    assert (
        report.statistics[
            "mixed_intraresidue_interresidue_groups_filtered"
        ]
        == 1
    )


def test_exclude_intraresidue_leaves_interresidue_union_unaffected():
    report = _convert(include_intraresidue=False)

    assert [group.group_id for group in report.ambiguous_groups] == [
        "restraints:4"
    ]
    alternatives = report.ambiguous_groups[0].alternatives
    assert len(alternatives) == 2
    assert all(
        alternative.atom1.chain != alternative.atom2.chain
        or alternative.atom1.residue_index
        != alternative.atom2.residue_index
        for alternative in alternatives
    )
    assert [constraint.source_groups for constraint in report.emitted_constraints] == [
        ["restraints:5"]
    ]


def test_cli_flag_filters_exact_and_union_outputs(tmp_path):
    output = tmp_path / "output"
    status = main(
        [
            "convert",
            str(FIXTURE),
            "--exclude-intraresidue",
            "--output-dir",
            str(output),
        ]
    )

    assert status == 0
    report = json.loads((output / "conversion_report.json").read_text())
    assert report["settings"]["include_intraresidue"] is False
    assert report["statistics"]["intraresidue_groups_filtered"] == 3
    exact = yaml.safe_load(
        (output / "atom_constraints_exact.yaml").read_text()
    )
    union = yaml.safe_load(
        (output / "atom_constraints_union.yaml").read_text()
    )
    assert len(exact["constraints"]) == 1
    assert len(union["constraints"]) == 1


@pytest.mark.parametrize("output_kind", ["exact", "union"])
def test_writer_fails_closed_if_intraresidue_excluded_report_is_mutated(
    tmp_path, output_kind
):
    report = _convert(include_intraresidue=False)
    if output_kind == "exact":
        item = report.emitted_constraints[0]
    else:
        item = report.ambiguous_groups[0].alternatives[0]
    item.atom2 = BoltzAtom(
        item.atom1.chain,
        item.atom1.residue_index,
        "CA",
    )
    destination = tmp_path / "invalid"

    with pytest.raises(
        ValueError,
        match="Intraresidue-excluded output validation failed",
    ):
        write_outputs(report, destination)

    assert not destination.exists()
