import gzip
from pathlib import Path

from nmr2boltz.model import BoltzAtom, ProjectedAlternative
from nmr2boltz.project import (
    ProjectionSettings,
    _merge_independent_constraints,
    _merge_or_alternatives,
    project_document,
)
from nmr2boltz.star import parse_star_document
from nmr2boltz.topology import TopologyLibrary

FIXTURES = Path(__file__).parent / "fixtures"


def projected(group: str, bound: float) -> ProjectedAlternative:
    return ProjectedAlternative(
        atom1=BoltzAtom("A", 1, "CA"),
        atom2=BoltzAtom("A", 8, "CB"),
        max_distance=bound,
        source_upper_bound=bound,
        averaging_policy="hard-or",
        averaging_factor=1.0,
        explicit_pair_count=1,
        bond_offset=0.0,
        group_id=group,
    )


def test_duplicate_pair_inside_or_uses_larger_bound():
    merged = _merge_or_alternatives([projected("g", 5.0), projected("g", 7.0)])
    assert len(merged) == 1
    assert merged[0].max_distance == 7.0


def test_duplicate_pair_across_independent_groups_uses_smaller_bound():
    emitted, rejected = _merge_independent_constraints(
        [projected("g1", 5.0), projected("g2", 7.0)], ProjectionSettings()
    )
    assert not rejected
    assert len(emitted) == 1
    assert emitted[0].max_distance == 5.0
    assert emitted[0].source_groups == ["g1", "g2"]


def test_boltz_minimum_weakens_and_maximum_is_not_clipped():
    emitted, rejected = _merge_independent_constraints(
        [projected("low", 1.5), projected("high-other-pair", 25.0)], ProjectionSettings()
    )
    # Both helper inputs target the same pair, so the independent conjunction is 1.5 A.
    assert len(emitted) == 1
    assert emitted[0].max_distance == 2.0
    assert not rejected

    over, over_rejected = _merge_independent_constraints([projected("over", 25.0)], ProjectionSettings())
    assert not over
    assert over_rejected[0].reason == "projected_bound_exceeds_boltz_maximum"
    assert "not clipped" in over_rejected[0].details


def test_compressed_nef_input(tmp_path):
    compressed = tmp_path / "example.nef.gz"
    with (FIXTURES / "example.nef").open("rb") as source, gzip.open(compressed, "wb") as target:
        target.write(source.read())
    parsed = parse_star_document(compressed)
    assert parsed.detected_format == "nef"
    assert len(parsed.restraint_groups) == 7


def test_user_residue_map_overrides_boltz_chain_and_index(tmp_path):
    mapping = tmp_path / "map.tsv"
    mapping.write_text(
        "source_chain\tsource_sequence_code\tsource_residue_name\tboltz_chain\tboltz_residue_index\n"
        "A\t10\tVAL\tZ\t101\n"
        "A\t10A\tALA\tZ\t102\n"
        "A\t12\tTYR\tZ\t103\n"
        "A\t13\tLEU\tZ\t104\n"
        "A\t14\tGLY\tZ\t105\n",
        encoding="utf-8",
    )
    parsed = parse_star_document(FIXTURES / "example.nef", residue_map_path=mapping)
    report = project_document(
        parsed,
        input_file=str(FIXTURES / "example.nef"),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )
    assert report.emitted_constraints
    assert {item.atom1.chain for item in report.emitted_constraints} == {"Z"}
    assert min(
        atom.residue_index
        for item in report.emitted_constraints
        for atom in (item.atom1, item.atom2)
    ) >= 101


def test_pseudoatom_atomset_mode_is_explicitly_labeled():
    choices = TopologyLibrary().resolve_expression("ALA", "QB", pseudoatom_policy="atomset")
    assert {atom.parent_atom for atom in choices[0].atoms} == {"CB"}
    assert "pseudoatom" in " ".join(choices[0].warnings)
