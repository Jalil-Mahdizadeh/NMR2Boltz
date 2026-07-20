from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import product
from typing import Any

from .model import (
    AmbiguousGroup,
    BoltzAtom,
    ConversionReport,
    EmittedConstraint,
    Endpoint,
    ProjectedAlternative,
    Rejection,
    RestraintGroup,
    SequenceRecord,
)
from .star import (
    ParsedStarDocument,
    StarDataError,
    normalize_nmrstar_canonical_expansions,
)
from .topology import (
    AtomSetChoice,
    TopologyLibrary,
    TopologyResolutionError,
    atom_topology_violations,
    component_topology_snapshot,
    mapped_residue_index,
    require_valid_emitted_atom_topology,
)


@dataclass
class ProjectionSettings:
    averaging_policy: str = "sum-r6"
    projection_margin: float = 0.0
    pseudoatom_policy: str = "reject"
    boltz_min_distance: float = 2.0
    boltz_max_distance: float = 20.0
    min_sequence_separation: int = 0
    include_intraresidue: bool = True
    include_intrachain: bool = True

    def validate(self) -> None:
        if self.averaging_policy not in {"sum-r6", "mean-r6", "hard-or"}:
            raise ValueError(f"Unsupported averaging policy: {self.averaging_policy}")
        if self.pseudoatom_policy not in {"reject", "atomset"}:
            raise ValueError(f"Unsupported pseudoatom policy: {self.pseudoatom_policy}")
        if not math.isfinite(self.projection_margin):
            raise ValueError("Projection margin must be finite.")
        if self.projection_margin < 0:
            raise ValueError("Projection margin must be non-negative.")
        if not math.isfinite(self.boltz_min_distance):
            raise ValueError("Boltz minimum distance must be finite.")
        if not math.isfinite(self.boltz_max_distance):
            raise ValueError("Boltz maximum distance must be finite.")
        if self.boltz_min_distance <= 0:
            raise ValueError("Boltz minimum distance must be positive.")
        if self.boltz_max_distance < self.boltz_min_distance:
            raise ValueError("Boltz maximum distance must be >= minimum distance.")
        if self.min_sequence_separation < 0:
            raise ValueError("Minimum sequence separation must be non-negative.")

    def as_dict(self) -> dict[str, Any]:
        return {
            "averaging_policy": self.averaging_policy,
            "projection_margin_angstrom": self.projection_margin,
            "pseudoatom_policy": self.pseudoatom_policy,
            "boltz_min_distance_angstrom": self.boltz_min_distance,
            "boltz_max_distance_angstrom": self.boltz_max_distance,
            "min_sequence_separation": self.min_sequence_separation,
            "include_intraresidue": self.include_intraresidue,
            "include_intrachain": self.include_intrachain,
            "exclude_intrachain": not self.include_intrachain,
        }


def _averaging_factor(policy: str, explicit_pair_count: int) -> float:
    if explicit_pair_count < 1:
        raise ValueError("Atom-set expansion produced zero explicit pairs.")
    if policy == "sum-r6":
        return explicit_pair_count ** (1.0 / 6.0)
    return 1.0


def _canonical_pair(atom1: BoltzAtom, atom2: BoltzAtom) -> tuple[BoltzAtom, BoltzAtom]:
    return tuple(sorted((atom1, atom2)))  # type: ignore[return-value]


def _residue_identity_mismatch(endpoint: Any, record: SequenceRecord) -> str | None:
    supplied = {
        str(name).upper()
        for name in (endpoint.residue_name, endpoint.canonical_residue_name)
        if name
    }
    mapped = {record.residue_name.upper()}
    mapped.update(str(alias[2]).upper() for alias in record.aliases if len(alias) >= 3)
    if not supplied or supplied & mapped:
        return None
    return (
        f"restraint endpoint {endpoint.display()} declares residue "
        f"{', '.join(sorted(supplied))}, but its resolved sequence record is "
        f"{record.boltz_chain}:{record.boltz_residue_index} "
        f"({', '.join(sorted(mapped))})"
    )


def project_document(
    parsed: ParsedStarDocument,
    *,
    input_file: str,
    topology_library: TopologyLibrary,
    settings: ProjectionSettings,
    parser_settings: dict[str, Any] | None = None,
) -> ConversionReport:
    settings.validate()
    for topology in parsed.embedded_topologies:
        topology_library.register(topology, replace=True)

    mapped_residues = mapped_residue_index(parsed.sequence_resolver.records)
    component_topologies = component_topology_snapshot(
        parsed.sequence_resolver.records, topology_library
    )

    safe_by_group: list[ProjectedAlternative] = []
    ambiguous: list[AmbiguousGroup] = []
    rejections: list[Rejection] = []
    warnings = list(parsed.warnings)
    intrachain_groups_filtered = 0
    projected_alternatives_removed_by_intrachain_filter = 0
    mixed_chain_scope_groups_filtered = 0
    intraresidue_groups_filtered = 0
    projected_alternatives_removed_by_intraresidue_filter = 0
    mixed_intraresidue_interresidue_groups_filtered = 0
    union_alternatives_raised_to_boltz_minimum = 0
    union_groups_quarantined_above_boltz_maximum = 0

    for group in parsed.restraint_groups:
        if group.complex_logic:
            rejections.append(
                Rejection(
                    group_id=group.group_id,
                    reason="complex_restraint_combination_logic",
                    details=(
                        "The restraint uses a non-null combination identifier. It was preserved in the "
                        "audit report but not flattened into a Boltz AND/OR expression."
                    ),
                    row_ids=_row_ids(group),
                )
            )
            continue
        normalization_issues = normalize_nmrstar_canonical_expansions(
            group, topology_library
        )
        if normalization_issues:
            rejections.extend(
                Rejection(
                    group_id=group.group_id,
                    reason="inconsistent_nmrstar_canonical_expansion",
                    details=issue.details,
                    row_ids=list(issue.row_ids),
                )
                for issue in normalization_issues
            )
            continue
        projected, group_rejections = _project_group(
            group,
            parsed=parsed,
            topology_library=topology_library,
            settings=settings,
        )
        projected, topology_rejections = _quarantine_invalid_projected_atoms(
            projected,
            mapped_residues=mapped_residues,
            component_topologies=component_topologies,
        )
        group_rejections.extend(topology_rejections)
        rejections.extend(group_rejections)
        if projected and group_rejections:
            rejections.append(
                Rejection(
                    group_id=group.group_id,
                    reason="partially_unprojectable_or_group",
                    details=(
                        "At least one source OR alternative or atom-set branch could not be projected. "
                        "The remaining alternatives were not emitted because doing so would strengthen "
                        "the source disjunction."
                    ),
                    row_ids=_row_ids(group),
                )
            )
            continue
        if not projected:
            if not group_rejections:
                rejections.append(
                    Rejection(
                        group_id=group.group_id,
                        reason="no_projectable_alternative",
                        details="No heavy-atom alternative could be constructed.",
                        row_ids=_row_ids(group),
                    )
                )
            continue

        merged = _merge_or_alternatives(projected)
        if not settings.include_intrachain:
            intrachain = [
                alternative
                for alternative in merged
                if alternative.atom1.chain == alternative.atom2.chain
            ]
            if intrachain:
                intrachain_groups_filtered += 1
                projected_alternatives_removed_by_intrachain_filter += len(merged)
                chain_scope = (
                    "intrachain"
                    if len(intrachain) == len(merged)
                    else "mixed_intrachain_interchain"
                )
                if chain_scope == "mixed_intrachain_interchain":
                    mixed_chain_scope_groups_filtered += 1
                rejections.append(
                    _intrachain_filter_rejection(
                        group,
                        merged,
                        chain_scope=chain_scope,
                    )
                )
                continue
        if not settings.include_intraresidue:
            intraresidue = [
                alternative
                for alternative in merged
                if (
                    alternative.atom1.chain == alternative.atom2.chain
                    and alternative.atom1.residue_index
                    == alternative.atom2.residue_index
                )
            ]
            if intraresidue:
                intraresidue_groups_filtered += 1
                projected_alternatives_removed_by_intraresidue_filter += len(
                    merged
                )
                residue_scope = (
                    "intraresidue"
                    if len(intraresidue) == len(merged)
                    else "mixed_intraresidue_interresidue"
                )
                if residue_scope == "mixed_intraresidue_interresidue":
                    mixed_intraresidue_interresidue_groups_filtered += 1
                rejections.append(
                    _intraresidue_filter_rejection(
                        group,
                        merged,
                        residue_scope=residue_scope,
                    )
                )
                continue
        if len(merged) == 1:
            candidate = merged[0]
            if (
                settings.min_sequence_separation > 0
                and candidate.atom1.chain == candidate.atom2.chain
                and abs(candidate.atom1.residue_index - candidate.atom2.residue_index)
                < settings.min_sequence_separation
            ):
                rejections.append(
                    Rejection(
                        group_id=group.group_id,
                        reason="sequence_separation_filtered",
                        details=(
                            f"Projected residue separation is below {settings.min_sequence_separation}."
                        ),
                        row_ids=candidate.source_rows,
                    )
                )
                continue
            safe_by_group.append(candidate)
        else:
            merged, bound_rejection, adjusted_count = _adapt_ambiguous_group_bounds(
                group,
                merged,
                settings,
            )
            union_alternatives_raised_to_boltz_minimum += adjusted_count
            if bound_rejection is not None:
                union_groups_quarantined_above_boltz_maximum += 1
                rejections.append(bound_rejection)
                continue
            ambiguous.append(
                AmbiguousGroup(
                    group_id=group.group_id,
                    restraint_id=group.restraint_id,
                    list_name=group.list_name,
                    alternatives=merged,
                    reason=(
                        "The original restraint is a disjunction whose alternatives map to multiple "
                        "distinct heavy-atom pairs. Emitting all pairs would convert OR into AND."
                    ),
                    source_format=group.source_format,
                    origin=group.origin,
                    warnings=list(group.warnings),
                )
            )

    emitted, final_rejections = _merge_independent_constraints(safe_by_group, settings)
    rejections.extend(final_rejections)

    statistics = {
        "restraint_groups_read": len(parsed.restraint_groups),
        "source_alternatives_read": sum(len(group.alternatives) for group in parsed.restraint_groups),
        "safe_groups_before_pair_deduplication": len(safe_by_group),
        "emitted_unique_heavy_atom_constraints": len(emitted),
        "ambiguous_or_groups_not_emitted": len(ambiguous),
        "rejection_records": len(rejections),
        "sequence_records": len(parsed.sequence_resolver.records),
        "embedded_component_topologies": len(parsed.embedded_topologies),
        "atom_topology_quarantines": sum(
            item.reason == "atom_not_present_in_mapped_residue" for item in rejections
        ),
        "intrachain_groups_filtered": intrachain_groups_filtered,
        "projected_alternatives_removed_by_intrachain_filter": (
            projected_alternatives_removed_by_intrachain_filter
        ),
        "mixed_chain_scope_groups_filtered": mixed_chain_scope_groups_filtered,
        "intraresidue_groups_filtered": intraresidue_groups_filtered,
        "projected_alternatives_removed_by_intraresidue_filter": (
            projected_alternatives_removed_by_intraresidue_filter
        ),
        "mixed_intraresidue_interresidue_groups_filtered": (
            mixed_intraresidue_interresidue_groups_filtered
        ),
        "union_alternatives_raised_to_boltz_minimum": (
            union_alternatives_raised_to_boltz_minimum
        ),
        "union_groups_quarantined_above_boltz_maximum": (
            union_groups_quarantined_above_boltz_maximum
        ),
        "emitted_atom_topology_violations": 0,
    }
    combined_settings = dict(parser_settings or {})
    combined_settings.update(settings.as_dict())
    report = ConversionReport(
        input_file=input_file,
        detected_format=parsed.detected_format,
        settings=combined_settings,
        sequence_map=parsed.sequence_resolver.records,
        emitted_constraints=emitted,
        ambiguous_groups=ambiguous,
        rejections=rejections,
        warnings=warnings,
        statistics=statistics,
        source_restraint_groups=parsed.restraint_groups,
        target_component_topologies=component_topologies,
    )
    require_valid_emitted_atom_topology(report)
    require_valid_interchain_output(report)
    require_valid_intraresidue_output(report)
    return report


def _intrachain_filter_rejection(
    group: RestraintGroup,
    alternatives: list[ProjectedAlternative],
    *,
    chain_scope: str,
) -> Rejection:
    projected = [
        {
            "atom1": {
                "chain": alternative.atom1.chain,
                "residue_index": alternative.atom1.residue_index,
                "atom_name": alternative.atom1.atom_name,
            },
            "atom2": {
                "chain": alternative.atom2.chain,
                "residue_index": alternative.atom2.residue_index,
                "atom_name": alternative.atom2.atom_name,
            },
            "scope": (
                "intrachain"
                if alternative.atom1.chain == alternative.atom2.chain
                else "interchain"
            ),
            "source_rows": list(alternative.source_rows),
        }
        for alternative in alternatives
    ]
    if chain_scope == "mixed_intrachain_interchain":
        details = (
            "The complete restraint group was removed by the inter-chain-only policy "
            "because its OR alternatives mix intrachain and inter-chain contacts. "
            "Keeping only the inter-chain alternatives would narrow and strengthen the "
            "source disjunction."
        )
    else:
        details = (
            "The complete projected restraint group was removed by the "
            "inter-chain-only policy because every alternative is intrachain."
        )
    return Rejection(
        group_id=group.group_id,
        reason="intrachain_filtered",
        details=details,
        row_ids=_row_ids(group),
        provenance={
            "filter": "exclude_intrachain",
            "chain_identity": "mapped_boltz_chain",
            "chain_scope": chain_scope,
            "projected_alternatives": projected,
        },
    )


def _intraresidue_filter_rejection(
    group: RestraintGroup,
    alternatives: list[ProjectedAlternative],
    *,
    residue_scope: str,
) -> Rejection:
    projected = [
        {
            "atom1": {
                "chain": alternative.atom1.chain,
                "residue_index": alternative.atom1.residue_index,
                "atom_name": alternative.atom1.atom_name,
            },
            "atom2": {
                "chain": alternative.atom2.chain,
                "residue_index": alternative.atom2.residue_index,
                "atom_name": alternative.atom2.atom_name,
            },
            "scope": (
                "intraresidue"
                if (
                    alternative.atom1.chain == alternative.atom2.chain
                    and alternative.atom1.residue_index
                    == alternative.atom2.residue_index
                )
                else "interresidue"
            ),
            "source_rows": list(alternative.source_rows),
        }
        for alternative in alternatives
    ]
    if residue_scope == "mixed_intraresidue_interresidue":
        details = (
            "The complete restraint group was removed by the "
            "intraresidue-exclusion policy because its OR alternatives mix "
            "intraresidue and inter-residue contacts. Keeping only the "
            "inter-residue alternatives would narrow and strengthen the "
            "source disjunction."
        )
    elif len(alternatives) == 1:
        # Preserve the established exact-contact audit wording.
        details = (
            "The projected restraint is intraresidue and was removed by policy."
        )
    else:
        details = (
            "The complete projected restraint group was removed by the "
            "intraresidue-exclusion policy because every alternative is "
            "intraresidue."
        )
    return Rejection(
        group_id=group.group_id,
        reason="intraresidue_filtered",
        details=details,
        row_ids=_row_ids(group),
        provenance={
            "filter": "exclude_intraresidue",
            "residue_identity": "mapped_boltz_chain_and_residue_index",
            "residue_scope": residue_scope,
            "projected_alternatives": projected,
        },
    )


def require_valid_interchain_output(report: ConversionReport) -> None:
    """Fail closed if an inter-chain-only report contains a same-chain endpoint pair."""
    if report.settings.get("include_intrachain", True):
        return
    violations: list[str] = []
    for constraint in report.emitted_constraints:
        if constraint.atom1.chain == constraint.atom2.chain:
            violations.append(
                f"exact {constraint.atom1.display()} -- {constraint.atom2.display()}"
            )
    for group in report.ambiguous_groups:
        for alternative in group.alternatives:
            if alternative.atom1.chain == alternative.atom2.chain:
                violations.append(
                    f"union {group.group_id} "
                    f"{alternative.atom1.display()} -- {alternative.atom2.display()}"
                )
    if violations:
        raise ValueError(
            "Inter-chain-only output validation failed: "
            + "; ".join(sorted(violations))
        )


def require_valid_intraresidue_output(report: ConversionReport) -> None:
    """Fail closed if an intraresidue-excluded report contains a local pair."""
    if report.settings.get("include_intraresidue", True):
        return
    violations: list[str] = []
    for constraint in report.emitted_constraints:
        if (
            constraint.atom1.chain == constraint.atom2.chain
            and constraint.atom1.residue_index
            == constraint.atom2.residue_index
        ):
            violations.append(
                f"exact {constraint.atom1.display()} -- "
                f"{constraint.atom2.display()}"
            )
    for group in report.ambiguous_groups:
        for alternative in group.alternatives:
            if (
                alternative.atom1.chain == alternative.atom2.chain
                and alternative.atom1.residue_index
                == alternative.atom2.residue_index
            ):
                violations.append(
                    f"union {group.group_id} "
                    f"{alternative.atom1.display()} -- "
                    f"{alternative.atom2.display()}"
                )
    if violations:
        raise ValueError(
            "Intraresidue-excluded output validation failed: "
            + "; ".join(sorted(violations))
        )


def _quarantine_invalid_projected_atoms(
    projected: list[ProjectedAlternative],
    *,
    mapped_residues: dict[tuple[str, int], str],
    component_topologies: dict[str, dict[str, object]],
) -> tuple[list[ProjectedAlternative], list[Rejection]]:
    """Remove complete projected contacts unless both atoms have topology evidence."""
    valid: list[ProjectedAlternative] = []
    rejections: list[Rejection] = []
    seen: set[tuple[Any, ...]] = set()
    for alternative in projected:
        violations = atom_topology_violations(
            (alternative.atom1, alternative.atom2),
            mapped_residues=mapped_residues,
            component_topologies=component_topologies,
        )
        if not violations:
            valid.append(alternative)
            continue
        invalid = [item.to_dict() for item in violations]
        rejection_key = (
            alternative.group_id,
            alternative.pair_key,
            tuple(alternative.source_rows),
            alternative.source_upper_bound,
            tuple(
                (
                    item["chain"],
                    item["residue_number"],
                    item["residue_name"],
                    item["atom_name"],
                    item["reason"],
                )
                for item in invalid
            ),
        )
        if rejection_key in seen:
            continue
        seen.add(rejection_key)
        invalid_labels = ", ".join(
            f"{item['chain']}:{item['residue_number']}:{item['residue_name'] or '?'}:"
            f"{item['atom_name']} ({item['reason']})"
            for item in invalid
        )
        original_bounds = {
            "lower_bound": alternative.source_observation.get("lower_bound"),
            "upper_bound": alternative.source_upper_bound,
            "target_value": alternative.source_observation.get("target_value"),
            "target_uncertainty": alternative.source_observation.get(
                "target_uncertainty"
            ),
            "upper_linear_limit": alternative.source_observation.get(
                "upper_linear_limit"
            ),
            "bound_source": alternative.source_observation.get("bound_source"),
        }
        rejections.append(
            Rejection(
                group_id=alternative.group_id,
                reason="atom_not_present_in_mapped_residue",
                details=(
                    "The complete projected contact was quarantined because exact component "
                    f"topology does not prove: {invalid_labels}."
                ),
                row_ids=list(alternative.source_rows),
                endpoint=" | ".join(alternative.source_endpoints) or None,
                provenance={
                    "restraint_group": alternative.group_id,
                    "source_row_ids": list(alternative.source_rows),
                    "source_endpoints": list(alternative.source_endpoints),
                    "mapped_contact": [
                        {
                            "chain": atom.chain,
                            "residue_number": atom.residue_index,
                            "residue_name": mapped_residues.get(
                                (atom.chain, atom.residue_index)
                            ),
                            "atom_name": atom.atom_name,
                        }
                        for atom in (alternative.atom1, alternative.atom2)
                    ],
                    "invalid_endpoints": invalid,
                    "original_bounds": original_bounds,
                    "projected_upper_bound": alternative.max_distance,
                },
            )
        )
    return valid, rejections


def _project_group(
    group: RestraintGroup,
    *,
    parsed: ParsedStarDocument,
    topology_library: TopologyLibrary,
    settings: ProjectionSettings,
) -> tuple[list[ProjectedAlternative], list[Rejection]]:
    projected: list[ProjectedAlternative] = []
    rejections: list[Rejection] = []
    for alternative in group.alternatives:
        if alternative.upper_bound is None:
            rejections.append(
                Rejection(
                    group_id=group.group_id,
                    reason="missing_upper_bound",
                    details=(
                        "No usable upper bound was available under the selected missing-upper policy."
                    ),
                    row_ids=alternative.row_ids,
                )
            )
            continue
        if not math.isfinite(alternative.upper_bound) or alternative.upper_bound <= 0:
            rejections.append(
                Rejection(
                    group_id=group.group_id,
                    reason="invalid_upper_bound",
                    details=f"Upper bound must be positive and finite; got {alternative.upper_bound!r}.",
                    row_ids=alternative.row_ids,
                )
            )
            continue
        if alternative.endpoint1.atom_expression is None or alternative.endpoint2.atom_expression is None:
            rejections.append(
                Rejection(
                    group_id=group.group_id,
                    reason="missing_atom_identifier",
                    details="One or both atom expressions are missing.",
                    row_ids=alternative.row_ids,
                )
            )
            continue
        try:
            residue1, map_warnings1 = parsed.sequence_resolver.resolve(alternative.endpoint1)
            residue2, map_warnings2 = parsed.sequence_resolver.resolve(alternative.endpoint2)
        except StarDataError as exc:
            rejections.append(
                Rejection(
                    group_id=group.group_id,
                    reason="unresolved_residue_mapping",
                    details=str(exc),
                    row_ids=alternative.row_ids,
                )
            )
            continue
        identity_errors = [
            message
            for message in (
                _residue_identity_mismatch(alternative.endpoint1, residue1),
                _residue_identity_mismatch(alternative.endpoint2, residue2),
            )
            if message is not None
        ]
        if identity_errors:
            rejections.append(
                Rejection(
                    group_id=group.group_id,
                    reason="sequence_residue_mismatch",
                    details="; ".join(identity_errors),
                    row_ids=alternative.row_ids,
                    endpoint=(
                        f"{alternative.endpoint1.display()} -- "
                        f"{alternative.endpoint2.display()}"
                    ),
                )
            )
            continue
        comp1 = alternative.endpoint1.canonical_residue_name or residue1.residue_name
        comp2 = alternative.endpoint2.canonical_residue_name or residue2.residue_name
        try:
            endpoint_choices1 = _resolve_endpoint_choices(
                topology_library,
                comp1,
                alternative.endpoint1,
                pseudoatom_policy=settings.pseudoatom_policy,
            )
            endpoint_choices2 = _resolve_endpoint_choices(
                topology_library,
                comp2,
                alternative.endpoint2,
                pseudoatom_policy=settings.pseudoatom_policy,
            )
        except TopologyResolutionError as exc:
            rejections.append(
                Rejection(
                    group_id=group.group_id,
                    reason="unresolved_atom_topology",
                    details=str(exc),
                    row_ids=alternative.row_ids,
                    endpoint=f"{alternative.endpoint1.display()} -- {alternative.endpoint2.display()}",
                )
            )
            continue

        source_observation = {
            "bound_source": alternative.bound_source,
            "explicit_or_derived_upper_bound": alternative.upper_bound,
            "lower_bound": alternative.lower_bound,
            "target_value": alternative.target_value,
            "target_uncertainty": alternative.target_uncertainty,
            "upper_linear_limit": alternative.upper_linear_limit,
            "weight": alternative.weight,
            "origin": alternative.origin,
            "combination_id": alternative.combination_id,
            "member_id": alternative.member_id,
            "member_logic_code": alternative.member_logic_code,
            "canonical_expansions": [
                dict(item) for item in alternative.canonical_expansions
            ],
        }

        for set1, set2 in product(endpoint_choices1, endpoint_choices2):
            explicit_pair_count = len(set1.atoms) * len(set2.atoms)
            factor = _averaging_factor(settings.averaging_policy, explicit_pair_count)
            by_parent_pair: dict[tuple[BoltzAtom, BoltzAtom], dict[str, Any]] = {}
            same_parent_pair_seen = False
            for atom1, atom2 in product(set1.atoms, set2.atoms):
                parent1 = BoltzAtom(
                    chain=residue1.boltz_chain,
                    residue_index=residue1.boltz_residue_index,
                    atom_name=atom1.parent_atom,
                )
                parent2 = BoltzAtom(
                    chain=residue2.boltz_chain,
                    residue_index=residue2.boltz_residue_index,
                    atom_name=atom2.parent_atom,
                )
                if parent1 == parent2:
                    same_parent_pair_seen = True
                    continue
                pair = _canonical_pair(parent1, parent2)
                offset = atom1.bond_length_upper + atom2.bond_length_upper
                threshold = (
                    factor * alternative.upper_bound
                    + offset
                    + settings.projection_margin
                )
                current = by_parent_pair.get(pair)
                detail = {
                    "threshold": threshold,
                    "offset": offset,
                    "atom_pair": f"{atom1.atom_name}--{atom2.atom_name}",
                }
                if current is None or threshold > current["threshold"]:
                    by_parent_pair[pair] = detail
            if same_parent_pair_seen:
                reason = (
                    "same_heavy_parent_atom"
                    if not by_parent_pair
                    else "atom_set_contains_same_heavy_parent_pair"
                )
                details = (
                    "Every expanded nucleus pair maps to the same heavy atom, which Boltz cannot "
                    "use as an atom-contact pair."
                    if not by_parent_pair
                    else (
                        "At least one explicit pair in this atom-set branch maps to the same heavy "
                        "parent. That branch provides no conservative constraint on the other parent "
                        "pairs, so none of them were emitted."
                    )
                )
                rejections.append(
                    Rejection(
                        group_id=group.group_id,
                        reason=reason,
                        details=details,
                        row_ids=alternative.row_ids,
                    )
                )
                continue
            if not by_parent_pair:
                continue
            common_warnings = list(dict.fromkeys(
                alternative.warnings
                + group.warnings
                + map_warnings1
                + map_warnings2
                + set1.warnings
                + set2.warnings
            ))
            source_endpoint = (
                f"{alternative.endpoint1.display()} -- {alternative.endpoint2.display()}"
            )
            for (parent1, parent2), detail in by_parent_pair.items():
                projected.append(
                    ProjectedAlternative(
                        atom1=parent1,
                        atom2=parent2,
                        max_distance=float(detail["threshold"]),
                        source_upper_bound=float(alternative.upper_bound),
                        averaging_policy=settings.averaging_policy,
                        averaging_factor=factor,
                        explicit_pair_count=explicit_pair_count,
                        bond_offset=float(detail["offset"]),
                        group_id=group.group_id,
                        source_observation=dict(source_observation),
                        source_rows=list(alternative.row_ids),
                        source_endpoints=[source_endpoint],
                        warnings=common_warnings,
                    )
                )
    return projected, rejections


def _resolve_endpoint_choices(
    topology_library: TopologyLibrary,
    comp_id: str,
    endpoint: Endpoint,
    *,
    pseudoatom_policy: str,
) -> list[AtomSetChoice]:
    if endpoint.canonical_atom_set:
        return [
            topology_library.resolve_canonical_atom_set(
                comp_id,
                endpoint.canonical_atom_set,
                author_expression=endpoint.atom_expression,
            )
        ]
    return topology_library.resolve_expression(
        comp_id,
        endpoint.atom_expression,
        canonical_hint=endpoint.canonical_atom_hint,
        pseudoatom_policy=pseudoatom_policy,
    )


def _merge_or_alternatives(
    alternatives: list[ProjectedAlternative],
) -> list[ProjectedAlternative]:
    """Merge duplicate heavy pairs inside one OR group.

    For the same pair, ``d <= u1 OR d <= u2`` is exactly ``d <= max(u1, u2)``.
    Using the minimum here would accidentally strengthen the original restraint.
    """
    merged: dict[tuple[BoltzAtom, BoltzAtom], ProjectedAlternative] = {}
    for alternative in alternatives:
        pair = alternative.pair_key
        if pair not in merged:
            merged[pair] = alternative
            continue
        current = merged[pair]
        if alternative.max_distance > current.max_distance:
            replacement = alternative
            replacement.source_rows = list(dict.fromkeys(current.source_rows + alternative.source_rows))
            replacement.source_endpoints = list(
                dict.fromkeys(current.source_endpoints + alternative.source_endpoints)
            )
            replacement.warnings = list(dict.fromkeys(current.warnings + alternative.warnings))
            merged[pair] = replacement
        else:
            current.source_rows = list(dict.fromkeys(current.source_rows + alternative.source_rows))
            current.source_endpoints = list(
                dict.fromkeys(current.source_endpoints + alternative.source_endpoints)
            )
            current.warnings = list(dict.fromkeys(current.warnings + alternative.warnings))
    return sorted(merged.values(), key=lambda item: item.pair_key)


def _adapt_ambiguous_group_bounds(
    group: RestraintGroup,
    alternatives: list[ProjectedAlternative],
    settings: ProjectionSettings,
) -> tuple[list[ProjectedAlternative], Rejection | None, int]:
    """Apply Boltz's executable interval without strengthening an OR group."""
    over_maximum = [
        alternative
        for alternative in alternatives
        if alternative.max_distance > settings.boltz_max_distance
    ]
    if over_maximum:
        maximum = max(alternative.max_distance for alternative in over_maximum)
        provenance = {
            "boltz_max_distance_angstrom": settings.boltz_max_distance,
            "complete_union_group_quarantined": True,
            "alternatives": [
                {
                    "atom1": {
                        "chain": alternative.pair_key[0].chain,
                        "residue_index": alternative.pair_key[0].residue_index,
                        "atom_name": alternative.pair_key[0].atom_name,
                    },
                    "atom2": {
                        "chain": alternative.pair_key[1].chain,
                        "residue_index": alternative.pair_key[1].residue_index,
                        "atom_name": alternative.pair_key[1].atom_name,
                    },
                    "projected_upper_bound_angstrom": alternative.max_distance,
                    "source_upper_bound_angstrom": alternative.source_upper_bound,
                    "source_rows": list(alternative.source_rows),
                }
                for alternative in alternatives
            ],
        }
        return (
            [],
            Rejection(
                group_id=group.group_id,
                reason="ambiguous_group_bound_exceeds_boltz_maximum",
                details=(
                    "At least one conservative union-alternative upper bound "
                    f"({maximum:.6g} A) exceeds the Boltz atom_contact maximum "
                    f"{settings.boltz_max_distance:.6g} A. The complete OR group "
                    "was quarantined; clipping or dropping only the incompatible "
                    "alternative would strengthen the source disjunction."
                ),
                row_ids=list(
                    dict.fromkeys(
                        row
                        for alternative in alternatives
                        for row in alternative.source_rows
                    )
                )
                or _row_ids(group),
                provenance=provenance,
            ),
            0,
        )

    adjusted_count = 0
    for alternative in alternatives:
        raw = alternative.max_distance
        if raw >= settings.boltz_min_distance:
            continue
        adjustment = (
            f"raised from {raw:.6g} to Boltz minimum "
            f"{settings.boltz_min_distance:.6g} A; this weakens rather than "
            "strengthens the union alternative"
        )
        alternative.raw_projected_distance = raw
        alternative.max_distance = settings.boltz_min_distance
        alternative.boltz_adjustment = adjustment
        alternative.warnings = list(
            dict.fromkeys([*alternative.warnings, adjustment])
        )
        adjusted_count += 1
    return alternatives, None, adjusted_count


def _merge_independent_constraints(
    projected: list[ProjectedAlternative],
    settings: ProjectionSettings,
) -> tuple[list[EmittedConstraint], list[Rejection]]:
    """Merge independent restraint groups that happen to target the same pair.

    Independent restraint groups are conjunctive. Thus ``d <= u1 AND d <= u2``
    becomes ``d <= min(u1, u2)``.
    """
    grouped: dict[tuple[BoltzAtom, BoltzAtom], list[ProjectedAlternative]] = {}
    for item in projected:
        grouped.setdefault(item.pair_key, []).append(item)
    emitted: list[EmittedConstraint] = []
    rejections: list[Rejection] = []
    for pair, items in sorted(grouped.items()):
        raw = min(item.max_distance for item in items)
        source_groups = list(dict.fromkeys(item.group_id for item in items))
        provenance = [
            {
                "group_id": item.group_id,
                "projected_upper_bound": item.max_distance,
                "source_upper_bound": item.source_upper_bound,
                "averaging_policy": item.averaging_policy,
                "averaging_factor": item.averaging_factor,
                "explicit_pair_count": item.explicit_pair_count,
                "bond_offset": item.bond_offset,
                "source_observation": dict(item.source_observation),
                "source_rows": item.source_rows,
                "source_endpoints": item.source_endpoints,
                "warnings": item.warnings,
            }
            for item in items
        ]
        adjustment: str | None = None
        if raw > settings.boltz_max_distance:
            rejections.append(
                Rejection(
                    group_id=",".join(source_groups),
                    reason="projected_bound_exceeds_boltz_maximum",
                    details=(
                        f"Conservative heavy-atom upper bound {raw:.6g} A exceeds Boltz atom_contact "
                        f"maximum {settings.boltz_max_distance:.6g} A. It was not clipped because "
                        "clipping would strengthen the restraint."
                    ),
                    row_ids=list(
                        dict.fromkeys(row for item in items for row in item.source_rows)
                    ),
                )
            )
            continue
        adjusted = raw
        if raw < settings.boltz_min_distance:
            adjusted = settings.boltz_min_distance
            adjustment = (
                f"raised from {raw:.6g} to Boltz minimum {settings.boltz_min_distance:.6g} A; "
                "this weakens rather than strengthens the restraint"
            )
        emitted.append(
            EmittedConstraint(
                atom1=pair[0],
                atom2=pair[1],
                max_distance=adjusted,
                source_groups=source_groups,
                raw_projected_distance=raw,
                boltz_adjustment=adjustment,
                provenance=provenance,
            )
        )
    return emitted, rejections


def _row_ids(group: RestraintGroup) -> list[str]:
    return list(
        dict.fromkeys(row for alternative in group.alternatives for row in alternative.row_ids)
    )
