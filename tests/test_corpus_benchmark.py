from types import SimpleNamespace

from nmr2boltz.model import BoltzAtom, EmittedConstraint
from validation.benchmark_corpus import compare_reports


def _constraint(atom1: BoltzAtom, atom2: BoltzAtom, distance: float) -> EmittedConstraint:
    return EmittedConstraint(
        atom1=atom1,
        atom2=atom2,
        max_distance=distance,
        source_groups=["g1"],
        raw_projected_distance=distance,
    )


def test_compare_reports_detects_pair_and_bound_differences():
    a = BoltzAtom("A", 1, "CA")
    b = BoltzAtom("A", 2, "CB")
    c = BoltzAtom("A", 3, "CG")
    nef = SimpleNamespace(emitted_constraints=[_constraint(a, b, 5.0)])
    star = SimpleNamespace(
        emitted_constraints=[_constraint(b, a, 5.5), _constraint(a, c, 7.0)]
    )

    result = compare_reports(nef, star)

    assert result["common_atom_pairs"] == 1
    assert result["star_only_atom_pairs"] == 1
    assert result["common_pairs_with_different_bounds"] == 1
    assert result["maximum_common_bound_delta_angstrom"] == 0.5
    assert result["exact_pair_and_bound_parity"] is False
