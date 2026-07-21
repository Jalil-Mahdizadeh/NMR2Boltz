from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

from .model import (
    AmbiguousGroup,
    BoltzAtom,
    BoltzToken,
    ConversionReport,
    EmittedConstraint,
    TokenConstraint,
    TokenContribution,
    TokenProjectionOmission,
)


TOKEN_MIN_DISTANCE = 4.0
TOKEN_MAX_DISTANCE = 20.0


class TokenProjectionValidationError(ValueError):
    """Raised when a token projection is unsafe to serialize."""


@dataclass
class TokenProjectionResult:
    constraints: list[TokenConstraint]
    omissions: list[TokenProjectionOmission]
    statistics: dict[str, int | float]


@dataclass
class _TokenCandidate:
    token1: BoltzToken
    token2: BoltzToken
    contribution: TokenContribution

    @property
    def pair_key(self) -> tuple[BoltzToken, BoltzToken]:
        return (self.token1, self.token2)


def _token(atom: BoltzAtom) -> BoltzToken:
    return BoltzToken(atom.chain, atom.residue_index)


def _canonical_token_pair(
    atom1: BoltzAtom, atom2: BoltzAtom
) -> tuple[BoltzToken, BoltzToken]:
    return tuple(sorted((_token(atom1), _token(atom2))))  # type: ignore[return-value]


def _atom_dict(atom: BoltzAtom) -> dict[str, Any]:
    return {
        "chain": atom.chain,
        "residue_index": atom.residue_index,
        "atom_name": atom.atom_name,
    }


def _projected_raw_bound(alternative: Any) -> float:
    if alternative.raw_projected_distance is not None:
        return alternative.raw_projected_distance
    return alternative.max_distance


def _candidate_adjustment(raw_bound: float) -> tuple[float, list[str]]:
    if raw_bound < TOKEN_MIN_DISTANCE:
        return TOKEN_MIN_DISTANCE, ["raised_to_token_minimum"]
    return raw_bound, []


def _bound_omission(
    *,
    source_kind: str,
    source_groups: list[str],
    raw_bound: float,
    provenance: dict[str, Any],
) -> TokenProjectionOmission | None:
    if not math.isfinite(raw_bound):
        return TokenProjectionOmission(
            source_kind=source_kind,
            source_groups=source_groups,
            reason="non_finite_token_bound",
            details="The projected token bound is not finite, so the complete semantic unit was omitted.",
            raw_bound=raw_bound,
            provenance=provenance,
        )
    if raw_bound > TOKEN_MAX_DISTANCE:
        return TokenProjectionOmission(
            source_kind=source_kind,
            source_groups=source_groups,
            reason="bound_exceeds_token_maximum",
            details=(
                f"The projected token bound {raw_bound:.6g} A exceeds the native "
                f"token-contact maximum {TOKEN_MAX_DISTANCE:.6g} A. It was omitted "
                "rather than clipped downward."
            ),
            raw_bound=raw_bound,
            provenance=provenance,
        )
    return None


def _exact_candidate(
    constraint: EmittedConstraint,
) -> tuple[_TokenCandidate | None, TokenProjectionOmission | None]:
    token1, token2 = _canonical_token_pair(constraint.atom1, constraint.atom2)
    raw_bound = constraint.raw_projected_distance
    provenance = {
        "atom1": _atom_dict(constraint.pair_key[0]),
        "atom2": _atom_dict(constraint.pair_key[1]),
        "atom_constraint_bound": constraint.max_distance,
        "raw_projected_distance": constraint.raw_projected_distance,
        "atom_constraint_adjustment": constraint.boltz_adjustment,
        "atom_constraint_provenance": constraint.provenance,
        "token_raw_bound": raw_bound,
    }
    source_groups = sorted(set(constraint.source_groups))
    if token1 == token2:
        return None, TokenProjectionOmission(
            source_kind="exact",
            source_groups=source_groups,
            reason="same_token",
            details="Both heavy-atom endpoints project to the same polymer token.",
            raw_bound=raw_bound,
            provenance=provenance,
        )
    omission = _bound_omission(
        source_kind="exact",
        source_groups=source_groups,
        raw_bound=raw_bound,
        provenance=provenance,
    )
    if omission is not None:
        return None, omission
    candidate_bound, adjustments = _candidate_adjustment(raw_bound)
    return (
        _TokenCandidate(
            token1=token1,
            token2=token2,
            contribution=TokenContribution(
                source_kind="exact",
                source_groups=source_groups,
                raw_bound=raw_bound,
                candidate_bound=candidate_bound,
                adjustments=adjustments,
                provenance=provenance,
            ),
        ),
        None,
    )


def _union_candidate(
    group: AmbiguousGroup,
) -> tuple[_TokenCandidate | None, TokenProjectionOmission | None]:
    alternatives = sorted(
        group.alternatives,
        key=lambda item: (item.pair_key, item.max_distance, tuple(item.source_rows)),
    )
    projected = [
        (_canonical_token_pair(item.atom1, item.atom2), item) for item in alternatives
    ]
    source_groups = [group.group_id]
    provenance = {
        "group_id": group.group_id,
        "restraint_id": group.restraint_id,
        "list_name": group.list_name,
        "source_format": group.source_format,
        "alternatives": [
            {
                "atom1": _atom_dict(item.pair_key[0]),
                "atom2": _atom_dict(item.pair_key[1]),
                "max_distance": item.max_distance,
                "raw_projected_distance": item.raw_projected_distance,
                "source_rows": list(item.source_rows),
                "source_endpoints": list(item.source_endpoints),
            }
            for item in alternatives
        ],
    }
    if any(pair[0] == pair[1] for pair, _ in projected):
        return None, TokenProjectionOmission(
            source_kind="collapsed_union",
            source_groups=source_groups,
            reason="union_contains_self_token_alternative",
            details=(
                "At least one OR alternative projects to a self-token contact, so the "
                "complete union was omitted."
            ),
            raw_bound=max(
                (_projected_raw_bound(item) for item in alternatives),
                default=None,
            ),
            provenance=provenance,
        )
    token_pairs = {pair for pair, _ in projected}
    if len(token_pairs) != 1:
        return None, TokenProjectionOmission(
            source_kind="collapsed_union",
            source_groups=source_groups,
            reason="union_spans_multiple_token_pairs",
            details=(
                "OR alternatives project to multiple token pairs. The complete union "
                "was omitted because ordinary contacts would turn OR into AND."
            ),
            raw_bound=max(
                (_projected_raw_bound(item) for item in alternatives),
                default=None,
            ),
            provenance=provenance,
        )
    token1, token2 = next(iter(token_pairs))
    raw_bound = max(_projected_raw_bound(item) for item in alternatives)
    provenance["token_raw_bound"] = raw_bound
    omission = _bound_omission(
        source_kind="collapsed_union",
        source_groups=source_groups,
        raw_bound=raw_bound,
        provenance=provenance,
    )
    if omission is not None:
        return None, omission
    candidate_bound, adjustments = _candidate_adjustment(raw_bound)
    return (
        _TokenCandidate(
            token1=token1,
            token2=token2,
            contribution=TokenContribution(
                source_kind="collapsed_union",
                source_groups=source_groups,
                raw_bound=raw_bound,
                candidate_bound=candidate_bound,
                adjustments=adjustments,
                provenance=provenance,
            ),
        ),
        None,
    )


def _candidate_sort_key(candidate: _TokenCandidate) -> tuple[Any, ...]:
    contribution = candidate.contribution
    return (
        candidate.pair_key,
        contribution.source_kind,
        tuple(contribution.source_groups),
        contribution.raw_bound,
    )


def project_token_constraints(report: ConversionReport) -> TokenProjectionResult:
    """Project canonical resolved atom constraints into native Boltz token contacts."""
    candidates: list[_TokenCandidate] = []
    omissions: list[TokenProjectionOmission] = []
    statistics = {
        "token_min_distance_angstrom": TOKEN_MIN_DISTANCE,
        "token_max_distance_angstrom": TOKEN_MAX_DISTANCE,
        "token_candidates": 0,
        "exact_candidates": 0,
        "collapsed_union_candidates": 0,
        "unique_token_constraints_emitted": 0,
        "same_token_omissions": 0,
        "multi_token_union_omissions": 0,
        "subminimum_adjustments": 0,
        "above_maximum_omissions": 0,
        "non_finite_bound_omissions": 0,
        "duplicate_token_pair_merges": 0,
        "token_projection_omissions": 0,
    }

    for constraint in sorted(report.emitted_constraints, key=lambda item: item.pair_key):
        candidate, omission = _exact_candidate(constraint)
        if candidate is not None:
            candidates.append(candidate)
            statistics["exact_candidates"] += 1
        if omission is not None:
            omissions.append(omission)

    for group in sorted(
        report.ambiguous_groups,
        key=lambda item: (item.group_id, item.restraint_id, item.list_name, item.source_format),
    ):
        candidate, omission = _union_candidate(group)
        if candidate is not None:
            candidates.append(candidate)
            statistics["collapsed_union_candidates"] += 1
        if omission is not None:
            omissions.append(omission)

    statistics["token_candidates"] = len(candidates)
    statistics["subminimum_adjustments"] = sum(
        "raised_to_token_minimum" in candidate.contribution.adjustments
        for candidate in candidates
    )
    statistics["same_token_omissions"] = sum(
        omission.reason in {"same_token", "union_contains_self_token_alternative"}
        for omission in omissions
    )
    statistics["multi_token_union_omissions"] = sum(
        omission.reason == "union_spans_multiple_token_pairs" for omission in omissions
    )
    statistics["above_maximum_omissions"] = sum(
        omission.reason == "bound_exceeds_token_maximum" for omission in omissions
    )
    statistics["non_finite_bound_omissions"] = sum(
        omission.reason == "non_finite_token_bound" for omission in omissions
    )
    statistics["token_projection_omissions"] = len(omissions)

    grouped: dict[tuple[BoltzToken, BoltzToken], list[_TokenCandidate]] = {}
    for candidate in sorted(candidates, key=_candidate_sort_key):
        grouped.setdefault(candidate.pair_key, []).append(candidate)

    constraints: list[TokenConstraint] = []
    for pair, items in sorted(grouped.items()):
        contributions = [item.contribution for item in items]
        adjustments = sorted(
            {
                adjustment
                for contribution in contributions
                for adjustment in contribution.adjustments
            }
        )
        if len(contributions) > 1:
            adjustments.append("duplicate_token_pair_merged")
            statistics["duplicate_token_pair_merges"] += len(contributions) - 1
        constraints.append(
            TokenConstraint(
                token1=pair[0],
                token2=pair[1],
                max_distance=min(item.candidate_bound for item in contributions),
                raw_candidate_bounds=[item.raw_bound for item in contributions],
                source_groups=sorted(
                    {
                        group
                        for contribution in contributions
                        for group in contribution.source_groups
                    }
                ),
                source_kinds=sorted({item.source_kind for item in contributions}),
                adjustments=adjustments,
                contributions=contributions,
            )
        )
    statistics["unique_token_constraints_emitted"] = len(constraints)
    return TokenProjectionResult(constraints, omissions, statistics)


def require_valid_token_projection(report: ConversionReport) -> None:
    """Validate all native-token invariants before an output bundle is staged."""
    resolved_tokens = {
        BoltzToken(record.boltz_chain, record.boltz_residue_index)
        for record in report.sequence_map
    }
    violations: list[str] = []
    seen: set[tuple[BoltzToken, BoltzToken]] = set()
    prior_pair: tuple[BoltzToken, BoltzToken] | None = None
    for constraint in report.token_constraints:
        pair = constraint.pair_key
        if (constraint.token1, constraint.token2) != pair:
            violations.append(
                f"non-canonical pair {constraint.token1.display()} -- {constraint.token2.display()}"
            )
        if pair[0] == pair[1]:
            violations.append(f"self pair {pair[0].display()}")
        for token in pair:
            if token not in resolved_tokens:
                violations.append(f"unresolved token {token.display()}")
        if pair in seen:
            violations.append(f"duplicate pair {pair[0].display()} -- {pair[1].display()}")
        seen.add(pair)
        if prior_pair is not None and pair < prior_pair:
            violations.append(
                f"non-deterministic ordering at {pair[0].display()} -- {pair[1].display()}"
            )
        prior_pair = pair
        if (
            not math.isfinite(constraint.max_distance)
            or constraint.max_distance < TOKEN_MIN_DISTANCE
            or constraint.max_distance > TOKEN_MAX_DISTANCE
        ):
            violations.append(
                f"bound {constraint.max_distance!r} A for "
                f"{pair[0].display()} -- {pair[1].display()}"
            )
    if violations:
        details = "; ".join(violations[:20])
        if len(violations) > 20:
            details += f"; and {len(violations) - 20} additional violation(s)"
        raise TokenProjectionValidationError(
            "Token projection invariants failed before output commit: " + details
        )
