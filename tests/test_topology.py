from nmr2boltz.topology import TopologyLibrary, TopologyResolutionError


def test_methyl_wildcard_collapses_to_one_parent():
    library = TopologyLibrary()
    choices = library.resolve_expression("VAL", "HG1%")
    assert len(choices) == 1
    assert {atom.atom_name for atom in choices[0].atoms} == {"HG11", "HG12", "HG13"}
    assert {atom.parent_atom for atom in choices[0].atoms} == {"CG1"}


def test_aromatic_wildcard_has_two_parents():
    library = TopologyLibrary()
    choices = library.resolve_expression("TYR", "HD%")
    assert len(choices) == 1
    assert {atom.parent_atom for atom in choices[0].atoms} == {"CD1", "CD2"}


def test_xy_branch_stays_assignment_ambiguous():
    library = TopologyLibrary()
    choices = library.resolve_expression("VAL", "HGx%")
    assert len(choices) == 2
    assert [{atom.parent_atom for atom in choice.atoms} for choice in choices] == [{"CG1"}, {"CG2"}]


def test_nucleotide_xy_does_not_consume_prime_notation():
    library = TopologyLibrary()

    choices = library.resolve_expression("DC", "H4y")

    assert {atom.atom_name for choice in choices for atom in choice.atoms} == {
        "H41",
        "H42",
    }
    assert {atom.parent_atom for choice in choices for atom in choice.atoms} == {"N4"}


def test_pseudoatom_rejected_by_default():
    library = TopologyLibrary()
    try:
        library.resolve_expression("ALA", "QB")
    except TopologyResolutionError as exc:
        assert "pseudoatom" in str(exc)
    else:
        raise AssertionError("QB should not be silently treated as an atom set")
