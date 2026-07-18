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
