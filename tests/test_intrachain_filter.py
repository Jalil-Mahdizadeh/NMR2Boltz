from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from nmr2boltz.cli import build_parser, main
from nmr2boltz.model import BoltzAtom
from nmr2boltz.output import write_outputs
from nmr2boltz.project import ProjectionSettings, project_document
from nmr2boltz.star import parse_star_document
from nmr2boltz.topology import TopologyLibrary


FIXTURES = Path(__file__).parent / "fixtures"


def _convert(name: str, *, include_intrachain: bool):
    path = FIXTURES / name
    return project_document(
        parse_star_document(path),
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(include_intrachain=include_intrachain),
    )


def _assert_only_interchain(report) -> None:
    assert report.emitted_constraints
    assert all(
        item.atom1.chain != item.atom2.chain
        for item in report.emitted_constraints
    )
    assert all(
        alternative.atom1.chain != alternative.atom2.chain
        for group in report.ambiguous_groups
        for alternative in group.alternatives
    )


@pytest.mark.parametrize("name", ["multichain.nef", "multichain.str"])
def test_exclude_intrachain_filters_exact_and_union_outputs_for_both_formats(name):
    report = _convert(name, include_intrachain=False)

    _assert_only_interchain(report)
    assert len(report.emitted_constraints) == 2
    assert len(report.ambiguous_groups) == 1
    assert report.settings["include_intrachain"] is False
    assert report.settings["exclude_intrachain"] is True
    assert report.statistics["intrachain_groups_filtered"] == 3
    assert (
        report.statistics["projected_alternatives_removed_by_intrachain_filter"]
        == 4
    )
    assert report.statistics["mixed_chain_scope_groups_filtered"] == 1

    filtered = [
        rejection
        for rejection in report.rejections
        if rejection.reason == "intrachain_filtered"
    ]
    assert len(filtered) == 3
    mixed = next(
        rejection
        for rejection in filtered
        if rejection.provenance["chain_scope"] == "mixed_intrachain_interchain"
    )
    assert mixed.row_ids
    assert {item["scope"] for item in mixed.provenance["projected_alternatives"]} == {
        "intrachain",
        "interchain",
    }
    assert "narrow and strengthen" in mixed.details

    nucleotide = next(
        constraint
        for constraint in report.emitted_constraints
        if {constraint.atom1.chain, constraint.atom2.chain} == {"D", "R"}
    )
    assert {nucleotide.atom1.atom_name, nucleotide.atom2.atom_name} == {"C8"}


@pytest.mark.parametrize("name", ["multichain.nef", "multichain.str"])
def test_default_keeps_intrachain_constraints(name):
    report = _convert(name, include_intrachain=True)

    assert any(
        item.atom1.chain == item.atom2.chain for item in report.emitted_constraints
    )
    assert any(
        alternative.atom1.chain == alternative.atom2.chain
        for group in report.ambiguous_groups
        for alternative in group.alternatives
    )
    assert report.statistics["intrachain_groups_filtered"] == 0


def test_cli_flag_writes_only_interchain_exact_and_union_yaml(tmp_path):
    output = tmp_path / "output"
    status = main(
        [
            "convert",
            str(FIXTURES / "multichain.nef"),
            "--exclude-intrachain",
            "--output-dir",
            str(output),
        ]
    )

    assert status == 0
    report = json.loads((output / "conversion_report.json").read_text())
    assert report["settings"]["include_intrachain"] is False
    assert report["settings"]["exclude_intrachain"] is True
    assert report["statistics"]["intrachain_groups_filtered"] == 3

    exact = yaml.safe_load((output / "atom_constraints_exact.yaml").read_text())
    union = yaml.safe_load((output / "atom_constraints_union.yaml").read_text())
    exact_pairs = [
        item["atom_contact"]["atom1"][0:1] + item["atom_contact"]["atom2"][0:1]
        for item in exact["constraints"]
    ]
    union_pairs = [
        alternative["atom1"][0:1] + alternative["atom2"][0:1]
        for item in union["constraints"]
        for alternative in item["atom_contact_union"]["alternatives"]
    ]
    assert exact_pairs and union_pairs
    assert all(chain1 != chain2 for chain1, chain2 in exact_pairs + union_pairs)


def test_cli_help_exposes_exact_flag_name():
    options = build_parser().parse_args(
        ["convert", str(FIXTURES / "multichain.nef"), "--exclude-intrachain"]
    )
    assert options.exclude_intrachain is True


def test_all_intrachain_input_writes_empty_constraint_files_and_strict_audit(tmp_path):
    output = tmp_path / "output"
    status = main(
        [
            "convert",
            str(FIXTURES / "example.nef"),
            "--exclude-intrachain",
            "--strict",
            "--output-dir",
            str(output),
        ]
    )

    assert status == 3
    for filename in ("atom_constraints_exact.yaml", "atom_constraints_union.yaml"):
        assert yaml.safe_load((output / filename).read_text()) == {"constraints": []}
    rejections = (output / "rejections.tsv").read_text()
    assert "intrachain_filtered" in rejections


def test_filter_uses_mapped_boltz_chain_not_source_chain(tmp_path):
    mapping = tmp_path / "residue-map.tsv"
    mapping.write_text(
        "source_chain\tsource_sequence_code\tsource_residue_name\t"
        "boltz_chain\tboltz_residue_index\n"
        "A\t1\tALA\tX\t1\n"
        "A\t2\tLEU\tX\t2\n"
        "B\t1\tGLY\tX\t3\n"
        "R\t1\tA\tR\t1\n"
        "D\t1\tDA\tD\t1\n",
        encoding="utf-8",
    )
    path = FIXTURES / "multichain.nef"
    report = project_document(
        parse_star_document(path, residue_map_path=mapping),
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(include_intrachain=False),
    )

    _assert_only_interchain(report)
    assert len(report.emitted_constraints) == 1
    assert not report.ambiguous_groups
    assert {report.emitted_constraints[0].atom1.chain,
            report.emitted_constraints[0].atom2.chain} == {"D", "R"}
    assert report.statistics["intrachain_groups_filtered"] == 5


@pytest.mark.parametrize("output_kind", ["exact", "union"])
def test_writer_fails_closed_if_interchain_report_is_mutated(
    tmp_path, output_kind
):
    report = _convert("multichain.nef", include_intrachain=False)
    if output_kind == "exact":
        item = report.emitted_constraints[0]
    else:
        item = report.ambiguous_groups[0].alternatives[0]
    item.atom2 = BoltzAtom(
        item.atom1.chain,
        item.atom2.residue_index,
        item.atom2.atom_name,
    )
    destination = tmp_path / "invalid"

    with pytest.raises(ValueError, match="Inter-chain-only output validation failed"):
        write_outputs(report, destination)

    assert not destination.exists()
