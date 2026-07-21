from __future__ import annotations

import copy
import json
import re
from pathlib import Path

import pytest
import yaml

from nmr2boltz.model import (
    AmbiguousGroup,
    BoltzAtom,
    BoltzToken,
    ConversionReport,
    EmittedConstraint,
    ProjectedAlternative,
    SequenceRecord,
)
from nmr2boltz.output import (
    _constraint_payload,
    _dump_yaml,
    _token_constraint_payload,
    _union_payload,
    write_outputs,
)
from nmr2boltz.project import ProjectionSettings, project_document
from nmr2boltz.star import parse_star_document
from nmr2boltz.token import (
    TokenProjectionValidationError,
    project_token_constraints,
)
from nmr2boltz.topology import TopologyLibrary


FIXTURES = Path(__file__).parent / "fixtures"


def _sequence(*tokens: tuple[str, int]) -> list[SequenceRecord]:
    return [
        SequenceRecord(
            source_chain=chain,
            source_sequence_code=str(index),
            residue_name="ALA",
            boltz_chain=chain,
            boltz_residue_index=index,
            source="test",
        )
        for chain, index in tokens
    ]


def _exact(
    group: str,
    bound: float,
    atom1: BoltzAtom | None = None,
    atom2: BoltzAtom | None = None,
) -> EmittedConstraint:
    atom1 = atom1 or BoltzAtom("A", 1, "CA")
    atom2 = atom2 or BoltzAtom("B", 2, "CB")
    return EmittedConstraint(
        atom1=atom1,
        atom2=atom2,
        max_distance=bound,
        source_groups=[group],
        raw_projected_distance=bound,
        provenance=[{"group_id": group}],
    )


def _alternative(
    group: str,
    bound: float,
    atom1: BoltzAtom,
    atom2: BoltzAtom,
    row: str,
) -> ProjectedAlternative:
    pair = tuple(sorted((atom1, atom2)))
    return ProjectedAlternative(
        atom1=pair[0],
        atom2=pair[1],
        max_distance=bound,
        source_upper_bound=bound,
        averaging_policy="hard-or",
        averaging_factor=1.0,
        explicit_pair_count=1,
        bond_offset=0.0,
        group_id=group,
        source_rows=[row],
    )


def _union(group: str, alternatives: list[ProjectedAlternative]) -> AmbiguousGroup:
    list_name, restraint_id = group.split(":", 1)
    return AmbiguousGroup(
        group_id=group,
        restraint_id=restraint_id,
        list_name=list_name,
        alternatives=alternatives,
        reason="test OR",
        source_format="nef",
    )


def _report(
    *,
    exact: list[EmittedConstraint] | None = None,
    unions: list[AmbiguousGroup] | None = None,
    tokens: tuple[tuple[str, int], ...] = (("A", 1), ("A", 2), ("B", 2), ("C", 3)),
) -> ConversionReport:
    return ConversionReport(
        input_file="test.nef",
        detected_format="nef",
        settings={},
        sequence_map=_sequence(*tokens),
        emitted_constraints=exact or [],
        ambiguous_groups=unions or [],
        rejections=[],
        warnings=[],
        statistics={},
    )


def _project(report: ConversionReport):
    result = project_token_constraints(report)
    report.token_constraints = result.constraints
    report.token_projection_omissions = result.omissions
    report.token_projection_statistics = result.statistics
    return result


def test_exact_atom_pair_projects_to_canonical_unordered_token_pair():
    result = _project(
        _report(
            exact=[
                _exact(
                    "g1",
                    6.0,
                    BoltzAtom("B", 2, "CB"),
                    BoltzAtom("A", 1, "CA"),
                )
            ]
        )
    )

    assert result.constraints[0].pair_key == (
        BoltzToken("A", 1),
        BoltzToken("B", 2),
    )
    assert result.constraints[0].source_kinds == ["exact"]


def test_same_token_exact_candidate_is_omitted_without_general_rejection():
    report = _report(
        exact=[
            _exact(
                "g1",
                6.0,
                BoltzAtom("A", 1, "CA"),
                BoltzAtom("A", 1, "CB"),
            )
        ]
    )
    result = _project(report)

    assert not result.constraints
    assert result.omissions[0].reason == "same_token"
    assert report.rejections == []


def test_independent_candidates_for_same_token_pair_use_minimum_bound():
    report = _report(
        exact=[
            _exact("g7", 7.0),
            _exact(
                "g5",
                5.0,
                BoltzAtom("A", 1, "CB"),
                BoltzAtom("B", 2, "CA"),
            ),
        ]
    )
    result = _project(report)

    assert len(result.constraints) == 1
    assert result.constraints[0].max_distance == 5.0
    assert result.constraints[0].source_groups == ["g5", "g7"]
    assert "duplicate_token_pair_merged" in result.constraints[0].adjustments


def test_same_token_pair_union_collapses_with_maximum_alternative_bound():
    group = "list:1"
    result = _project(
        _report(
            unions=[
                _union(
                    group,
                    [
                        _alternative(
                            group,
                            5.0,
                            BoltzAtom("A", 1, "CA"),
                            BoltzAtom("B", 2, "CA"),
                            "r1",
                        ),
                        _alternative(
                            group,
                            7.0,
                            BoltzAtom("A", 1, "CB"),
                            BoltzAtom("B", 2, "CB"),
                            "r2",
                        ),
                    ],
                )
            ]
        )
    )

    assert result.constraints[0].max_distance == 7.0
    assert result.constraints[0].source_kinds == ["collapsed_union"]
    assert result.statistics["collapsed_union_candidates"] == 1


def test_union_spanning_multiple_token_pairs_is_omitted_intact():
    group = "list:1"
    result = _project(
        _report(
            unions=[
                _union(
                    group,
                    [
                        _alternative(
                            group,
                            5.0,
                            BoltzAtom("A", 1, "CA"),
                            BoltzAtom("B", 2, "CA"),
                            "r1",
                        ),
                        _alternative(
                            group,
                            6.0,
                            BoltzAtom("A", 1, "CB"),
                            BoltzAtom("C", 3, "CA"),
                            "r2",
                        ),
                    ],
                )
            ]
        )
    )

    assert not result.constraints
    assert [item.reason for item in result.omissions] == [
        "union_spans_multiple_token_pairs"
    ]


def test_mixed_self_and_nonself_union_is_omitted_intact():
    group = "list:1"
    result = _project(
        _report(
            unions=[
                _union(
                    group,
                    [
                        _alternative(
                            group,
                            5.0,
                            BoltzAtom("A", 1, "CA"),
                            BoltzAtom("A", 1, "CB"),
                            "r1",
                        ),
                        _alternative(
                            group,
                            6.0,
                            BoltzAtom("A", 1, "CB"),
                            BoltzAtom("B", 2, "CA"),
                            "r2",
                        ),
                    ],
                )
            ]
        )
    )

    assert not result.constraints
    assert result.omissions[0].reason == "union_contains_self_token_alternative"


def test_exact_and_collapsed_union_candidates_merge_with_minimum_bound():
    group = "list:1"
    result = _project(
        _report(
            exact=[_exact("exact:1", 8.0)],
            unions=[
                _union(
                    group,
                    [
                        _alternative(
                            group,
                            5.0,
                            BoltzAtom("A", 1, "CB"),
                            BoltzAtom("B", 2, "CA"),
                            "r1",
                        ),
                        _alternative(
                            group,
                            6.0,
                            BoltzAtom("A", 1, "CA"),
                            BoltzAtom("B", 2, "CB"),
                            "r2",
                        ),
                    ],
                )
            ],
        )
    )

    assert result.constraints[0].max_distance == 6.0
    assert result.constraints[0].source_kinds == ["collapsed_union", "exact"]


def test_subminimum_bound_is_raised_and_audited():
    constraint = _exact("g1", 3.25)
    constraint.max_distance = 10.0
    result = _project(_report(exact=[constraint]))

    assert result.constraints[0].max_distance == 4.0
    assert result.constraints[0].raw_candidate_bounds == [3.25]
    assert result.constraints[0].contributions[0].adjustments == [
        "raised_to_token_minimum"
    ]
    assert result.statistics["subminimum_adjustments"] == 1


def test_above_maximum_bound_is_omitted_and_never_clipped():
    result = _project(_report(exact=[_exact("g1", 20.000001)]))

    assert not result.constraints
    assert result.omissions[0].reason == "bound_exceeds_token_maximum"
    assert result.omissions[0].raw_bound == 20.000001


def test_collapsible_union_above_maximum_is_omitted_as_one_semantic_unit():
    group = "list:1"
    result = _project(
        _report(
            unions=[
                _union(
                    group,
                    [
                        _alternative(
                            group,
                            5.0,
                            BoltzAtom("A", 1, "CA"),
                            BoltzAtom("B", 2, "CA"),
                            "r1",
                        ),
                        _alternative(
                            group,
                            20.000001,
                            BoltzAtom("A", 1, "CB"),
                            BoltzAtom("B", 2, "CB"),
                            "r2",
                        ),
                    ],
                )
            ]
        )
    )

    assert not result.constraints
    assert len(result.omissions) == 1
    assert result.omissions[0].source_groups == [group]
    assert result.omissions[0].reason == "bound_exceeds_token_maximum"


def test_projection_is_side_effect_free_for_exact_and_union_atom_outputs():
    group = "list:1"
    report = _report(
        exact=[_exact("exact:1", 6.0)],
        unions=[
            _union(
                group,
                [
                    _alternative(
                        group,
                        5.0,
                        BoltzAtom("A", 1, "CA"),
                        BoltzAtom("B", 2, "CA"),
                        "r1",
                    ),
                    _alternative(
                        group,
                        7.0,
                        BoltzAtom("A", 1, "CB"),
                        BoltzAtom("C", 3, "CA"),
                        "r2",
                    ),
                ],
            )
        ],
    )
    exact_before = copy.deepcopy(_constraint_payload(report.emitted_constraints))
    union_before = copy.deepcopy(_union_payload(report.ambiguous_groups[0]))

    project_token_constraints(report)

    assert _constraint_payload(report.emitted_constraints) == exact_before
    assert _union_payload(report.ambiguous_groups[0]) == union_before


def test_token_yaml_is_deterministic_outward_rounded_and_force_false():
    first_report = _report(
        exact=[
            _exact(
                "g2",
                7.0,
                BoltzAtom("A", 2, "CA"),
                BoltzAtom("C", 3, "CB"),
            ),
            _exact("g1", 6.7200001),
        ]
    )
    second_report = copy.deepcopy(first_report)
    second_report.emitted_constraints.reverse()
    first = _project(first_report)
    second = _project(second_report)

    first_text = _dump_yaml(_token_constraint_payload(first.constraints))
    second_text = _dump_yaml(_token_constraint_payload(second.constraints))
    payload = yaml.safe_load(first_text)

    assert first_text == second_text
    assert "max_distance: 6.720001" in first_text
    assert re.findall(r"max_distance: (\d+\.\d{6})", first_text)
    assert all(item["contact"]["force"] is False for item in payload["constraints"])
    assert payload["constraints"][0]["contact"]["token1"] == ["A", 1]


def test_token_yaml_empty_output_is_explicit():
    assert _dump_yaml(_token_constraint_payload([])) == "constraints: []\n"


@pytest.fixture(scope="module")
def topology_library() -> TopologyLibrary:
    return TopologyLibrary()


@pytest.mark.parametrize("fixture", ["example.nef", "example.str"])
def test_nef_and_nmrstar_end_to_end_write_token_audits(
    tmp_path, fixture, topology_library
):
    path = FIXTURES / fixture
    report = project_document(
        parse_star_document(path),
        input_file=str(path),
        topology_library=topology_library,
        settings=ProjectionSettings(),
    )

    write_outputs(report, tmp_path)

    assert (tmp_path / "atom_constraints_exact.yaml").is_file()
    assert (tmp_path / "atom_constraints_union.yaml").is_file()
    assert (tmp_path / "token_constraints.yaml").is_file()
    assert (tmp_path / "token_constraints.tsv").is_file()
    tsv_lines = (tmp_path / "token_constraints.tsv").read_text().splitlines()
    assert "source_kind" in tsv_lines[0].split("\t")
    assert any(kind in "\n".join(tsv_lines[1:]) for kind in ("exact", "collapsed_union"))
    summary = (tmp_path / "summary.txt").read_text()
    assert "Token candidates:" in summary
    assert "Collapsed-union token candidates:" in summary
    audit = json.loads((tmp_path / "conversion_report.json").read_text())
    assert set(audit) >= {
        "token_constraints",
        "token_projection_omissions",
        "token_projection_statistics",
    }
    assert audit["token_projection_statistics"]["token_candidates"] >= len(
        audit["token_constraints"]
    )


def test_invalid_token_invariant_fails_before_output_is_committed(
    tmp_path, topology_library
):
    path = FIXTURES / "example.nef"
    report = project_document(
        parse_star_document(path),
        input_file=str(path),
        topology_library=topology_library,
        settings=ProjectionSettings(),
    )
    report.token_constraints[0].token1 = BoltzToken("Z", 999)
    destination = tmp_path / "invalid-token"

    with pytest.raises(TokenProjectionValidationError, match="before output commit"):
        write_outputs(report, destination)

    assert not destination.exists()
