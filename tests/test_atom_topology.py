from dataclasses import asdict
from pathlib import Path

import pytest
import yaml

from nmr2boltz.model import BoltzAtom, Endpoint, RawAlternative, RestraintGroup, SequenceRecord
from nmr2boltz.output import write_outputs
from nmr2boltz.project import ProjectionSettings, project_document
from nmr2boltz.star import ParsedStarDocument, SequenceResolver, parse_star_document
from nmr2boltz.topology import (
    AtomTopologyValidationError,
    ComponentTopology,
    TopologyLibrary,
    TopologyResolutionError,
    emitted_atom_topology_violations,
)


def _document(
    residues: list[tuple[str, str, str]],
    contacts: list[tuple[str, str, str, str, str, str, str, str]],
) -> ParsedStarDocument:
    resolver = SequenceResolver()
    for chain, sequence_code, residue_name in residues:
        resolver.add(
            SequenceRecord(
                source_chain=chain,
                source_sequence_code=sequence_code,
                residue_name=residue_name,
                boltz_chain=chain,
                boltz_residue_index=int(sequence_code),
                source="synthetic-sequence",
            )
        )
    groups = []
    for index, contact in enumerate(contacts, start=1):
        chain1, seq1, residue1, atom1, chain2, seq2, residue2, atom2 = contact
        groups.append(
            RestraintGroup(
                source_format="nef",
                list_name="topology_regression",
                restraint_id=str(index),
                alternatives=[
                    RawAlternative(
                        source_format="nef",
                        list_name="topology_regression",
                        restraint_id=str(index),
                        endpoint1=Endpoint(chain1, seq1, residue1, atom1),
                        endpoint2=Endpoint(chain2, seq2, residue2, atom2),
                        upper_bound=4.0,
                        lower_bound=1.8,
                        row_ids=[f"row-{index}"],
                    )
                ],
            )
        )
    return ParsedStarDocument(
        entry=None,
        detected_format="nef",
        sequence_resolver=resolver,
        restraint_groups=groups,
        embedded_topologies=[],
        warnings=[],
    )


def _project(parsed: ParsedStarDocument, library: TopologyLibrary | None = None):
    return project_document(
        parsed,
        input_file="synthetic.nef",
        topology_library=library or TopologyLibrary(),
        settings=ProjectionSettings(),
    )


def test_glutamine_zinc_atom_is_quarantined_with_complete_provenance():
    parsed = _document(
        [("A", "1", "ALA"), ("B", "48", "GLN")],
        [("A", "1", "ALA", "CA", "B", "48", "GLN", "ZN")],
    )

    report = _project(parsed)

    assert not report.emitted_constraints
    rejection = next(
        item for item in report.rejections
        if item.reason == "atom_not_present_in_mapped_residue"
    )
    assert rejection.group_id == "topology_regression:1"
    assert rejection.row_ids == ["row-1"]
    assert rejection.provenance["invalid_endpoints"] == [
        {
            "chain": "B",
            "residue_number": 48,
            "residue_name": "GLN",
            "atom_name": "ZN",
            "topology_source": "builtin-standard",
            "reason": "atom_absent_from_component_topology",
        }
    ]
    assert rejection.provenance["original_bounds"]["upper_bound"] == 4.0
    assert rejection.provenance["original_bounds"]["lower_bound"] == 1.8


@pytest.mark.parametrize(
    ("residue", "atom"),
    [("LEU", "N3"), ("LEU", "N4"), ("LEU", "O2"),
     ("TRP", "N1"), ("TRP", "N2"), ("TRP", "O6")],
)
def test_protein_residue_rejects_projected_nucleotide_atoms(residue, atom):
    parsed = _document(
        [("A", "14", residue), ("B", "1", "ALA")],
        [("A", "14", residue, atom, "B", "1", "ALA", "CA")],
    )

    report = _project(parsed)

    assert not report.emitted_constraints
    assert any(
        item.reason == "atom_not_present_in_mapped_residue"
        and item.provenance["invalid_endpoints"][0]["atom_name"] == atom
        for item in report.rejections
    )


def test_standard_amino_acid_and_nucleotide_heavy_atoms_still_emit():
    parsed = _document(
        [("P", "1", "GLN"), ("P", "2", "TRP"), ("R", "1", "A"), ("R", "2", "G")],
        [
            ("P", "1", "GLN", "OE1", "P", "2", "TRP", "NE1"),
            ("R", "1", "A", "N6", "R", "2", "G", "O6"),
        ],
    )

    report = _project(parsed)

    assert len(report.emitted_constraints) == 2
    assert emitted_atom_topology_violations(report) == []


@pytest.mark.parametrize(("comp_id", "atom_name"), [("MSE", "SE"), ("ZN", "ZN"), ("LIG", "C1")])
def test_available_modified_ligand_and_ion_topology_emits(comp_id, atom_name):
    library = TopologyLibrary()
    if comp_id not in {"MSE"}:
        topology = ComponentTopology(comp_id=comp_id, source="declared-test-ccd")
        topology.add_atom(atom_name, "ZN" if comp_id == "ZN" else "C")
        library.register(topology)
    parsed = _document(
        [("L", "1", comp_id), ("A", "1", "ALA")],
        [("L", "1", comp_id, atom_name, "A", "1", "ALA", "CA")],
    )

    report = _project(parsed, library)

    assert len(report.emitted_constraints) == 1
    assert emitted_atom_topology_violations(report) == []


def test_unknown_component_topology_fails_closed():
    parsed = _document(
        [("X", "1", "NOCCD"), ("A", "1", "ALA")],
        [("X", "1", "NOCCD", "C1", "A", "1", "ALA", "CA")],
    )

    report = _project(parsed)

    assert not report.emitted_constraints
    rejection = next(
        item for item in report.rejections
        if item.reason == "atom_not_present_in_mapped_residue"
    )
    assert rejection.provenance["invalid_endpoints"][0]["reason"] == "component_topology_unavailable"


def test_distinct_same_component_source_residues_cannot_share_boltz_position():
    resolver = SequenceResolver()
    for source_sequence, boltz_index, residue in (
        ("1", 1, "ALA"),
        ("2", 1, "ALA"),
        ("3", 2, "GLY"),
    ):
        resolver.add(
            SequenceRecord(
                source_chain="A",
                source_sequence_code=source_sequence,
                residue_name=residue,
                boltz_chain="A",
                boltz_residue_index=boltz_index,
                source="synthetic-collision",
            )
        )
    parsed = ParsedStarDocument(
        entry=None,
        detected_format="nef",
        sequence_resolver=resolver,
        restraint_groups=[],
        embedded_topologies=[],
        warnings=[],
    )

    with pytest.raises(
        AtomTopologyValidationError,
        match=r"Distinct source residues A:1 \(ALA\) and A:2 \(ALA\).*A:1",
    ):
        _project(parsed)


def test_hydrogen_with_multiple_heavy_parents_fails_closed_deterministically():
    def topology(bonds):
        component = ComponentTopology("SYN", source="synthetic-ccd")
        for atom, element in (("H1", "H"), ("C1", "C"), ("N1", "N")):
            component.add_atom(atom, element)
        for atom1, atom2 in bonds:
            component.add_bond(atom1, atom2)
        return component

    forward = topology([("H1", "C1"), ("H1", "N1")])
    reverse = topology([("H1", "N1"), ("H1", "C1")])

    assert forward.hydrogen_parent_conflicts["H1"] == ("C1", "N1")
    assert reverse.hydrogen_parent_conflicts["H1"] == ("C1", "N1")
    for component in (forward, reverse):
        assert "H1" not in component.hydrogen_parent
        with pytest.raises(
            TopologyResolutionError,
            match=r"H1 has multiple heavy-atom parents.*C1, N1",
        ):
            component.atom_choice("H1", "test")


def test_canonical_component_alias_topology_is_valid_at_mapped_position():
    resolver = SequenceResolver()
    resolver.add(
        SequenceRecord(
            source_chain="A",
            source_sequence_code="1",
            residue_name="HSD",
            boltz_chain="A",
            boltz_residue_index=1,
            source="synthetic-nmrstar",
            aliases=[("A", "1", "HIS")],
        )
    )
    resolver.add(
        SequenceRecord("A", "2", "ALA", "A", 2, "synthetic-nmrstar")
    )
    group = RestraintGroup(
        source_format="nmr-star",
        list_name="canonical-component",
        restraint_id="1",
        alternatives=[
            RawAlternative(
                source_format="nmr-star",
                list_name="canonical-component",
                restraint_id="1",
                endpoint1=Endpoint(
                    "A", "1", "HSD", "HD2", "A", "1", "HIS", "HD2"
                ),
                endpoint2=Endpoint(
                    "A", "2", "ALA", "HA", "A", "2", "ALA", "HA"
                ),
                upper_bound=4.0,
                row_ids=["1"],
            )
        ],
    )
    parsed = ParsedStarDocument(
        entry=None,
        detected_format="nmr-star",
        sequence_resolver=resolver,
        restraint_groups=[group],
        embedded_topologies=[],
        warnings=[],
    )

    report = _project(parsed)

    assert len(report.emitted_constraints) == 1
    assert not any(
        rejection.reason == "atom_not_present_in_mapped_residue"
        for rejection in report.rejections
    )
    assert emitted_atom_topology_violations(report) == []


def test_external_ccd_without_readable_atom_table_is_not_silently_ignored(
    tmp_path,
):
    ccd = tmp_path / "bad.cif"
    ccd.write_text("data_BAD\n_chem_comp.id BAD\n", encoding="utf-8")
    library = TopologyLibrary(external_ccd_paths=[ccd])

    with pytest.raises(
        TopologyResolutionError,
        match=r"BAD.*no readable _chem_comp_atom",
    ):
        library.get("BAD")


def test_external_ccd_two_column_bond_table_still_loads(tmp_path):
    ccd = tmp_path / "lig.cif"
    ccd.write_text(
        "data_LIG\n"
        "_chem_comp.id LIG\n"
        "loop_\n"
        "_chem_comp_atom.atom_id\n"
        "_chem_comp_atom.type_symbol\n"
        "H1 H\n"
        "C1 C\n"
        "loop_\n"
        "_chem_comp_bond.atom_id_1\n"
        "_chem_comp_bond.atom_id_2\n"
        "H1 C1\n",
        encoding="utf-8",
    )
    topology = TopologyLibrary(external_ccd_paths=[ccd]).get("LIG")

    assert topology is not None
    assert topology.hydrogen_parent["H1"] == "C1"


def test_final_output_validator_fails_before_writing_invalid_yaml(tmp_path):
    parsed = _document(
        [("A", "1", "ALA"), ("A", "2", "GLY")],
        [("A", "1", "ALA", "CB", "A", "2", "GLY", "CA")],
    )
    report = _project(parsed)
    report.emitted_constraints[0].atom2 = BoltzAtom("A", 2, "ZN")

    with pytest.raises(AtomTopologyValidationError, match="topology validation failed"):
        write_outputs(report, tmp_path / "should-not-exist")

    assert not (tmp_path / "should-not-exist").exists()


def test_rejection_provenance_and_yaml_atom_membership_are_deterministic(tmp_path):
    bad = _document(
        [("A", "1", "ALA"), ("B", "48", "GLN")],
        [("A", "1", "ALA", "CA", "B", "48", "GLN", "ZN")],
    )
    first = _project(bad)
    second = _project(bad)
    assert [asdict(item) for item in first.rejections] == [
        asdict(item) for item in second.rejections
    ]

    good = _document(
        [("A", "1", "ALA"), ("A", "2", "GLY")],
        [("A", "1", "ALA", "CB", "A", "2", "GLY", "CA")],
    )
    report = _project(good)
    write_outputs(report, tmp_path / "valid")
    payload = yaml.safe_load(
        (tmp_path / "valid" / "atom_constraints_exact.yaml").read_text()
    )
    assert payload["constraints"][0]["atom_contact"]["atom1"] == ["A", 1, "CB"]
    assert emitted_atom_topology_violations(report) == []


@pytest.mark.parametrize(
    ("case_id", "suffix", "expected_rejections"),
    [("8R1X", "nef", 4), ("8R1X", "str", 4), ("9CCH", "nef", 8), ("9CCH", "str", 8)],
)
def test_deposited_corpus_atom_defects_are_quarantined(case_id, suffix, expected_rejections):
    corpus = Path(__file__).parents[1] / "benchmark" / "input" / case_id
    path = next(corpus.glob(f"*.{suffix}"))

    report = _project(parse_star_document(path))
    quarantined = [
        item for item in report.rejections
        if item.reason == "atom_not_present_in_mapped_residue"
    ]

    assert len(quarantined) == expected_rejections
    assert emitted_atom_topology_violations(report) == []
    for constraint in report.emitted_constraints:
        for atom in (constraint.atom1, constraint.atom2):
            assert not (
                case_id == "9CCH"
                and atom.chain == "B"
                and atom.residue_index == 48
                and atom.atom_name.startswith("ZN")
            )
            assert not (
                case_id == "8R1X"
                and (
                    (atom.chain, atom.residue_index, atom.atom_name)
                    in {
                        ("A", 14, "N3"), ("A", 14, "N4"), ("A", 14, "O2"),
                        ("A", 15, "N1"), ("A", 15, "N2"), ("A", 15, "O6"),
                    }
                )
            )
