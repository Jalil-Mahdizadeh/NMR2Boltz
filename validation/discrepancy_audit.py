"""Build row-level, fail-closed audits of paired NEF/NMR-STAR contacts."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from itertools import product
from pathlib import Path
from typing import Any

from nmr2boltz.topology import TopologyLibrary, TopologyResolutionError


BOUND_TOLERANCE = 1e-6
AUDIT_FIELDS = [
    "audit_id",
    "case_id",
    "discrepancy_type",
    "classification",
    "rationale_code",
    "rationale",
    "projected_heavy_pair",
    "nef_final_bound_angstrom",
    "star_final_bound_angstrom",
    "absolute_bound_delta_angstrom",
    "nef_source_restraint_ids",
    "star_source_restraint_ids",
    "nef_source_row_ids",
    "star_source_row_ids",
    "nef_atom_expressions",
    "star_atom_expressions",
    "nef_canonical_expansions",
    "star_canonical_expansions",
    "nef_physical_proton_sets",
    "star_physical_proton_sets",
    "nef_pseudoatom_handling",
    "star_pseudoatom_handling",
    "nef_pair_counts_N",
    "star_pair_counts_N",
    "nef_averaging_policies",
    "star_averaging_policies",
    "nef_source_upper_bounds_angstrom",
    "star_source_upper_bounds_angstrom",
    "nef_projected_bounds_angstrom",
    "star_projected_bounds_angstrom",
    "nef_rejection_reasons",
    "star_rejection_reasons",
]


def _report_dict(report: Any) -> dict[str, Any]:
    return report if isinstance(report, dict) else report.to_dict()


def _atom_key(atom: dict[str, Any]) -> tuple[str, int, str]:
    return str(atom["chain"]), int(atom["residue_index"]), str(atom["atom_name"])


def _constraint_key(constraint: dict[str, Any]) -> tuple[tuple[str, int, str], ...]:
    return tuple(sorted((_atom_key(constraint["atom1"]), _atom_key(constraint["atom2"]))))


def _pair_label(key: tuple[tuple[str, int, str], ...]) -> str:
    return "--".join(":".join(map(str, atom)) for atom in key)


def _normalize_group_id(group_id: str) -> str:
    list_name, separator, restraint_id = group_id.rpartition(":")
    if not separator:
        return group_id.casefold()
    prefix = "nef_distance_restraint_list_"
    if list_name.casefold().startswith(prefix):
        list_name = list_name[len(prefix) :]
    return f"{list_name.casefold()}:{restraint_id}"


def _group_id(group: dict[str, Any]) -> str:
    return f'{group["list_name"]}:{group["restraint_id"]}'


def _group_indexes(report: dict[str, Any]) -> tuple[dict[str, dict[str, Any]], dict[str, list[str]]]:
    groups = {
        _normalize_group_id(_group_id(group)): group
        for group in report["source_restraint_groups"]
    }
    rejections: dict[str, list[str]] = {}
    for rejection in report["rejections"]:
        key = _normalize_group_id(str(rejection["group_id"]))
        rejections.setdefault(key, []).append(str(rejection["reason"]))
    return groups, rejections


def _ordered_unique(values: list[Any]) -> list[Any]:
    return list(dict.fromkeys(values))


def _sequence_indexes(
    report: dict[str, Any],
) -> tuple[dict[tuple[str, str, str], tuple[str, int]], dict[tuple[str, str], set[tuple[str, int]]]]:
    exact: dict[tuple[str, str, str], tuple[str, int]] = {}
    loose: dict[tuple[str, str], set[tuple[str, int]]] = {}
    for record in report["sequence_map"]:
        destination = (str(record["boltz_chain"]), int(record["boltz_residue_index"]))
        aliases = [
            (
                str(record.get("source_chain") or ""),
                str(record.get("source_sequence_code") or ""),
                str(record.get("residue_name") or "").upper(),
            )
        ]
        aliases.extend(tuple(str(part) for part in alias) for alias in record.get("aliases", []))
        for chain, sequence, residue in aliases:
            # Match SequenceResolver: ambiguous aliases retain the first exact
            # record unless an explicit override was requested during parsing.
            exact.setdefault((chain, sequence, residue.upper()), destination)
            loose.setdefault((chain, sequence), set()).add(destination)
    return exact, loose


def _endpoint_destination(
    endpoint: dict[str, Any],
    exact: dict[tuple[str, str, str], tuple[str, int]],
    loose: dict[tuple[str, str], set[tuple[str, int]]],
) -> tuple[str, int] | None:
    candidates = [
        (
            str(endpoint.get("chain_code") or ""),
            str(endpoint.get("sequence_code") or ""),
            str(endpoint.get("residue_name") or "").upper(),
        ),
        (
            str(endpoint.get("canonical_chain_code") or ""),
            str(endpoint.get("canonical_sequence_code") or ""),
            str(endpoint.get("canonical_residue_name") or "").upper(),
        ),
    ]
    for candidate in candidates:
        if candidate in exact:
            return exact[candidate]
    for chain, sequence, _residue in candidates:
        destinations = loose.get((chain, sequence), set())
        if len(destinations) == 1:
            return next(iter(destinations))
    return None


def _physical_set(endpoint: dict[str, Any], topology: TopologyLibrary) -> tuple[str, str]:
    residue = str(
        endpoint.get("canonical_residue_name") or endpoint.get("residue_name") or ""
    ).upper()
    expression = str(endpoint.get("atom_expression") or "")
    if not expression:
        return "ERROR: missing atom expression", "unresolved"
    try:
        canonical_atom_set = endpoint.get("canonical_atom_set") or []
        if canonical_atom_set:
            choices = [
                topology.resolve_canonical_atom_set(
                    residue,
                    canonical_atom_set,
                    author_expression=expression,
                )
            ]
        else:
            choices = topology.resolve_expression(
                residue,
                expression,
                canonical_hint=endpoint.get("canonical_atom_hint"),
                pseudoatom_policy="reject",
            )
    except TopologyResolutionError as exc:
        message = str(exc)
        handling = (
            "rejected_geometric_pseudoatom"
            if "pseudoatom" in message.casefold()
            else "unresolved_atom_expression"
        )
        return f"ERROR: {message}", handling
    rendered = [
        {
            "assignment": choice.assignment_key,
            "semantics": choice.semantics,
            "atoms": [atom.atom_name for atom in choice.atoms],
            "parents": _ordered_unique([atom.parent_atom for atom in choice.atoms]),
        }
        for choice in choices
    ]
    if endpoint.get("canonical_atom_set"):
        handling = "physical_atom_set"
    elif any(symbol in expression for symbol in ("%", "*", "#")):
        handling = "physical_atom_set"
    elif expression.endswith(("x", "y")):
        handling = "stereospecific_assignment_alternative"
    else:
        handling = "explicit_atom_or_alias"
    return json.dumps(rendered, separators=(",", ":"), sort_keys=True), handling


def _physical_atoms(
    rendered: str, destination: tuple[str, int] | None
) -> tuple[tuple[str, int, str, str], ...]:
    """Return canonical proton/parent identities from a resolved endpoint."""
    if rendered.startswith("ERROR:") or destination is None:
        return ()
    atoms: set[tuple[str, int, str, str]] = set()
    for choice in json.loads(rendered):
        names = [str(atom).upper() for atom in choice.get("atoms", [])]
        parents = [str(atom).upper() for atom in choice.get("parents", [])]
        if len(parents) == 1:
            atoms.update(
                (destination[0], destination[1], name, parents[0]) for name in names
            )
        elif len(parents) == len(names):
            atoms.update(
                (destination[0], destination[1], name, parent)
                for name, parent in zip(names, parents)
            )
    return tuple(sorted(atoms))


def _semantic_kinds(rendered: str, handling: str) -> tuple[str, ...]:
    kinds = {handling}
    if not rendered.startswith("ERROR:"):
        for choice in json.loads(rendered):
            semantics = str(choice.get("semantics") or "")
            if "atom-set" in semantics:
                kinds.add("physical_atom_set")
            if "stereo" in semantics:
                kinds.add("stereospecific_assignment_alternative")
    return tuple(sorted(kinds))


def _physical_pairs(
    atoms1: tuple[tuple[str, int, str, str], ...],
    atoms2: tuple[tuple[str, int, str, str], ...],
) -> tuple[tuple[tuple[str, int, str, str], tuple[str, int, str, str]], ...]:
    return tuple(sorted({tuple(sorted(pair)) for pair in product(atoms1, atoms2)}))


def _format_evidence(
    report: dict[str, Any],
    constraint: dict[str, Any] | None,
    matched_keys: set[str],
    groups: dict[str, dict[str, Any]],
    rejections: dict[str, list[str]],
) -> dict[str, Any]:
    topology = TopologyLibrary()
    exact_sequence, loose_sequence = _sequence_indexes(report)
    selected = [groups[key] for key in sorted(matched_keys) if key in groups]
    restraint_ids: list[str] = []
    row_ids: list[str] = []
    expressions: list[str] = []
    canonical: list[str] = []
    physical: list[str] = []
    handling: list[str] = []
    source_bounds: list[float] = []
    representation: list[Any] = []
    semantic_groups: dict[str, dict[str, Any]] = {}
    destinations: dict[str, set[tuple[tuple[str, int] | None, tuple[str, int] | None]]] = {}
    for group in selected:
        group_label = _group_id(group)
        normalized_group = _normalize_group_id(group_label)
        restraint_ids.append(group_label)
        group_representation: list[Any] = []
        semantic_alternatives: list[dict[str, Any]] = []
        group_destinations: set[
            tuple[tuple[str, int] | None, tuple[str, int] | None]
        ] = set()
        for alternative in group["alternatives"]:
            endpoint1 = alternative["endpoint1"]
            endpoint2 = alternative["endpoint2"]
            destination1 = _endpoint_destination(
                endpoint1, exact_sequence, loose_sequence
            )
            destination2 = _endpoint_destination(
                endpoint2, exact_sequence, loose_sequence
            )
            expression = (
                f'{endpoint1.get("atom_expression")}--{endpoint2.get("atom_expression")}'
            )
            expansion_records = alternative.get("canonical_expansions") or [
                {
                    "row_id": ",".join(str(row) for row in alternative.get("row_ids", [])),
                    "atom1": endpoint1.get("canonical_atom_hint"),
                    "atom2": endpoint2.get("canonical_atom_hint"),
                }
            ]
            expansion_signature = tuple(
                sorted(
                    (
                        str(item.get("row_id") or ""),
                        str(item.get("atom1") or "."),
                        str(item.get("atom2") or "."),
                    )
                    for item in expansion_records
                )
            )
            set1, handling1 = _physical_set(endpoint1, topology)
            set2, handling2 = _physical_set(endpoint2, topology)
            atoms1 = _physical_atoms(set1, destination1)
            atoms2 = _physical_atoms(set2, destination2)
            kinds1 = _semantic_kinds(set1, handling1)
            kinds2 = _semantic_kinds(set2, handling2)
            rows = [str(row) for row in alternative.get("row_ids", [])]
            row_ids.extend(rows)
            expressions.append(f"{group_label}[{','.join(rows)}]:{expression}")
            canonical.extend(
                f'{group_label}[{item.get("row_id") or "."}]:'
                f'{item.get("atom1") or "."}--{item.get("atom2") or "."}'
                for item in expansion_records
            )
            physical.append(f"{group_label}[{','.join(rows)}]:{set1}--{set2}")
            handling.extend((handling1, handling2))
            if alternative.get("upper_bound") is not None:
                source_bounds.append(float(alternative["upper_bound"]))
            group_representation.append(
                (expression, expansion_signature, set1, set2)
            )
            semantic_alternatives.append(
                {
                    "raw_signature": tuple(
                        sorted(
                            (
                                (
                                    destination1,
                                    str(endpoint1.get("atom_expression") or "").upper(),
                                    str(endpoint1.get("canonical_atom_hint") or "").upper(),
                                    set1,
                                ),
                                (
                                    destination2,
                                    str(endpoint2.get("atom_expression") or "").upper(),
                                    str(endpoint2.get("canonical_atom_hint") or "").upper(),
                                    set2,
                                ),
                            ),
                            key=repr,
                        )
                    ),
                    "physical_pairs": _physical_pairs(atoms1, atoms2),
                    "kinds": tuple(sorted(set(kinds1) | set(kinds2))),
                    "unverified_heavy_atoms": tuple(
                        sorted(
                            atom
                            for rendered, atoms in ((set1, atoms1), (set2, atoms2))
                            if "explicit-heavy-unlisted" in rendered
                            for atom in atoms
                        )
                    ),
                }
            )
            group_destinations.add(
                tuple(
                    sorted(
                        (destination1, destination2),
                        key=lambda item: ("", 0) if item is None else item,
                    )
                )  # type: ignore[arg-type]
            )
        representation.append(
            (normalized_group, tuple(group_representation))
        )
        destinations[normalized_group] = group_destinations
        semantic_groups[normalized_group] = {
            "alternatives": tuple(semantic_alternatives),
            "raw_signature": tuple(
                sorted(
                    (alternative["raw_signature"] for alternative in semantic_alternatives),
                    key=repr,
                )
            ),
            "physical_alternative_signature": tuple(
                sorted(
                    (
                        alternative["physical_pairs"]
                        for alternative in semantic_alternatives
                    ),
                    key=repr,
                )
            ),
            "physical_pairs": frozenset(
                pair
                for alternative in semantic_alternatives
                for pair in alternative["physical_pairs"]
            ),
            "kinds": frozenset(
                kind
                for alternative in semantic_alternatives
                for kind in alternative["kinds"]
            ),
            "unverified_heavy_atoms": frozenset(
                atom
                for alternative in semantic_alternatives
                for atom in alternative["unverified_heavy_atoms"]
            ),
            "rejection_reasons": frozenset(rejections.get(normalized_group, [])),
        }

    provenance = constraint.get("provenance", []) if constraint else []
    pair_counts = [int(item["explicit_pair_count"]) for item in provenance]
    averaging = [str(item["averaging_policy"]) for item in provenance]
    projected_bounds = [float(item["projected_upper_bound"]) for item in provenance]
    return {
        "source_restraint_ids": restraint_ids,
        "source_row_ids": _ordered_unique(row_ids),
        "atom_expressions": expressions,
        "canonical_expansions": canonical,
        "physical_proton_sets": physical,
        "pseudoatom_handling": _ordered_unique(handling),
        "pair_counts_N": pair_counts,
        "averaging_policies": _ordered_unique(averaging),
        "source_upper_bounds": sorted(set(source_bounds)),
        "projected_bounds": projected_bounds,
        "rejection_reasons": sorted(
            set(reason for key in matched_keys for reason in rejections.get(key, []))
        ),
        "present_group_keys": sorted(key for key in matched_keys if key in groups),
        "representation": representation,
        "semantic_groups": semantic_groups,
        "destinations": destinations,
        "final_bound": float(constraint["max_distance"]) if constraint else None,
    }


def _bounds_by_group(
    matched_keys: set[str], groups: dict[str, dict[str, Any]]
) -> dict[str, tuple[float, ...]]:
    return {
        key: tuple(
            sorted(
                {
                    float(alt["upper_bound"])
                    for alt in groups[key]["alternatives"]
                    if alt.get("upper_bound") is not None
                }
            )
        )
        for key in matched_keys
        if key in groups
    }


def _wildcard_set_vs_explicit_or(
    first: dict[str, Any], second: dict[str, Any]
) -> bool:
    """Prove wildcard atom-set syntax versus explicit OR-member syntax."""
    for wildcard, explicit in ((first, second), (second, first)):
        if wildcard["rejection_reasons"] or explicit["rejection_reasons"]:
            continue
        if "physical_atom_set" not in wildcard["kinds"]:
            continue
        if wildcard["kinds"] & {
            "stereospecific_assignment_alternative",
            "rejected_geometric_pseudoatom",
            "unresolved_atom_expression",
        }:
            continue
        if explicit["kinds"] & {
            "stereospecific_assignment_alternative",
            "rejected_geometric_pseudoatom",
            "unresolved_atom_expression",
        }:
            continue
        wildcard_pairs = wildcard["physical_pairs"]
        explicit_pairs = explicit["physical_pairs"]
        atomset_alternatives = [
            set(pairs)
            for pairs in wildcard["physical_alternative_signature"]
            if len(pairs) > 1
        ]
        member_alternatives = [
            set(pairs) for pairs in explicit["physical_alternative_signature"]
        ]
        if (
            wildcard_pairs
            and explicit_pairs
            and explicit_pairs <= wildcard_pairs
            and atomset_alternatives
            and member_alternatives
            and all(member <= wildcard_pairs for member in member_alternatives)
            and (
                explicit["kinds"] == {"explicit_atom_or_alias"}
                or len(explicit["alternatives"]) > len(wildcard["alternatives"])
                or (
                    "explicit_atom_or_alias" in explicit["kinds"]
                    and explicit_pairs < wildcard_pairs
                )
            )
        ):
            return True
    return False


def _stereo_assignment_vs_physical_set(
    first: dict[str, Any], second: dict[str, Any]
) -> bool:
    """Prove x/y assignment possibilities versus a resolved physical atom set."""
    allowed_stereo_rejections = {
        frozenset(),
        frozenset(
            {"partially_unprojectable_or_group", "same_heavy_parent_atom"}
        ),
    }
    for stereo, physical in ((first, second), (second, first)):
        if (
            stereo["rejection_reasons"] not in allowed_stereo_rejections
            or physical["rejection_reasons"]
        ):
            continue
        if "stereospecific_assignment_alternative" not in stereo["kinds"]:
            continue
        if "rejected_geometric_pseudoatom" in stereo["kinds"]:
            continue
        if physical["kinds"] & {
            "stereospecific_assignment_alternative",
            "rejected_geometric_pseudoatom",
            "unresolved_atom_expression",
        }:
            continue
        stereo_pairs = stereo["physical_pairs"]
        physical_pairs = physical["physical_pairs"]
        if (
            stereo_pairs
            and physical_pairs
            and (physical_pairs <= stereo_pairs or stereo_pairs <= physical_pairs)
        ):
            return True
    return False


def _rejected_geometric_pseudoatom(
    first: dict[str, Any], second: dict[str, Any]
) -> bool:
    """Prove a conservative Q/M pseudoatom rejection against resolved evidence."""
    for pseudoatom, resolved in ((first, second), (second, first)):
        if "rejected_geometric_pseudoatom" not in pseudoatom["kinds"]:
            continue
        if "rejected_geometric_pseudoatom" in resolved["kinds"]:
            continue
        if pseudoatom["rejection_reasons"] != {"unresolved_atom_topology"}:
            continue
        if resolved["rejection_reasons"] or not resolved["physical_pairs"]:
            continue
        return True
    return False


def _verified_canonical_alias(
    first: dict[str, Any], second: dict[str, Any]
) -> bool:
    """Prove different explicit names resolve to the same atoms and OR structure."""
    explicit = {"explicit_atom_or_alias"}
    return bool(
        first["kinds"] == explicit
        and second["kinds"] == explicit
        and first["physical_pairs"]
        and first["physical_pairs"] == second["physical_pairs"]
        and first["physical_alternative_signature"]
        == second["physical_alternative_signature"]
        and first["raw_signature"] != second["raw_signature"]
        and not first["rejection_reasons"]
        and not second["rejection_reasons"]
    )


ALLOWLIST_PREDICATES = (
    ("wildcard_set_vs_explicit_or", _wildcard_set_vs_explicit_or),
    ("stereo_assignment_vs_physical_set", _stereo_assignment_vs_physical_set),
    ("rejected_geometric_pseudoatom", _rejected_geometric_pseudoatom),
    ("verified_canonical_alias", _verified_canonical_alias),
)


def _strictly_equivalent(first: dict[str, Any], second: dict[str, Any]) -> bool:
    if first["rejection_reasons"] != second["rejection_reasons"]:
        return False
    if first["raw_signature"] == second["raw_signature"]:
        return True
    # A NEF atom-set expression and a topology-verified reconstruction of the
    # corresponding NMR-STAR canonical rows are semantically identical even
    # though their raw spellings differ. Require the complete OR structure,
    # physical atoms, and semantic kinds to agree; pair-set equality alone
    # would incorrectly accept superficially similar representations.
    return bool(
        not first["rejection_reasons"]
        and first["physical_pairs"]
        and first["physical_pairs"] == second["physical_pairs"]
        and first["physical_alternative_signature"]
        == second["physical_alternative_signature"]
        and first["kinds"] == second["kinds"]
    )


def _classify(
    matched_keys: set[str],
    nef_groups: dict[str, dict[str, Any]],
    star_groups: dict[str, dict[str, Any]],
    nef: dict[str, Any],
    star: dict[str, Any],
) -> tuple[str, str, str]:
    rejection_reasons = set(nef["rejection_reasons"]) | set(star["rejection_reasons"])
    if rejection_reasons & {"sequence_residue_mismatch", "unresolved_residue_mapping"}:
        return (
            "deposition_inconsistency",
            "sequence_or_residue_mapping_conflict",
            "A corresponding deposited restraint conflicts with its deposited sequence or residue mapping.",
        )
    missing_nef = matched_keys - set(nef["present_group_keys"])
    missing_star = matched_keys - set(star["present_group_keys"])
    if missing_nef or missing_star:
        return (
            "deposition_inconsistency",
            "missing_corresponding_source_restraint",
            "A contributing source restraint ID is absent from the counterpart deposition.",
        )
    if _bounds_by_group(matched_keys, nef_groups) != _bounds_by_group(
        matched_keys, star_groups
    ):
        return (
            "deposition_inconsistency",
            "source_upper_bound_difference",
            "Corresponding restraint IDs contain different deposited upper bounds.",
        )
    if any(nef["destinations"].get(key) != star["destinations"].get(key) for key in matched_keys):
        return (
            "deposition_inconsistency",
            "source_residue_numbering_difference",
            "Corresponding restraint IDs resolve to different target residues in the two deposited sequence maps.",
        )
    unverified_heavy_differences = [
        key
        for key in matched_keys
        if nef["semantic_groups"][key]["unverified_heavy_atoms"]
        != star["semantic_groups"][key]["unverified_heavy_atoms"]
        and (
            nef["semantic_groups"][key]["unverified_heavy_atoms"]
            or star["semantic_groups"][key]["unverified_heavy_atoms"]
        )
    ]
    if unverified_heavy_differences:
        return (
            "deposition_inconsistency",
            "unverified_heavy_atom_name_difference",
            "Corresponding source rows use different unverified heavy-atom names in groups "
            + ", ".join(sorted(unverified_heavy_differences))
            + "; topology cannot prove them equivalent.",
        )
    allowed: set[str] = set()
    unrecognized: list[str] = []
    equivalent: list[str] = []
    for key in sorted(matched_keys):
        nef_semantics = nef["semantic_groups"][key]
        star_semantics = star["semantic_groups"][key]
        if _strictly_equivalent(nef_semantics, star_semantics):
            equivalent.append(key)
            continue
        matched_predicates = [
            code
            for code, predicate in ALLOWLIST_PREDICATES
            if predicate(nef_semantics, star_semantics)
        ]
        if matched_predicates:
            allowed.update(matched_predicates)
        else:
            unrecognized.append(key)
    if unrecognized:
        return (
            "unresolved",
            "unrecognized_semantic_difference",
            "Source groups "
            + ", ".join(unrecognized)
            + " differ but satisfy no tested semantic allowlist predicate.",
        )
    if allowed:
        return (
            "expected_format_difference",
            "+".join(sorted(allowed)),
            "Every non-equivalent source group is proven by an explicit semantic allowlist predicate: "
            + ", ".join(sorted(allowed))
            + ".",
        )
    if equivalent and len(equivalent) == len(matched_keys):
        return (
            "parser_projection_bug",
            "equivalent_source_evidence_changed_projection",
            "Equivalent source evidence produced a different projected pair or final bound.",
        )
    return (
        "unresolved",
        "unclassified_projection_difference",
        "The discrepancy satisfies neither a deposition rule nor a tested semantic predicate.",
    )


def audit_reports(case_id: str, nef_report: Any, star_report: Any) -> list[dict[str, Any]]:
    """Return one complete evidence row per pair or bound discrepancy."""
    nef_dict = _report_dict(nef_report)
    star_dict = _report_dict(star_report)
    nef_constraints = {
        _constraint_key(item): item for item in nef_dict["emitted_constraints"]
    }
    star_constraints = {
        _constraint_key(item): item for item in star_dict["emitted_constraints"]
    }
    nef_groups, nef_rejections = _group_indexes(nef_dict)
    star_groups, star_rejections = _group_indexes(star_dict)
    rows: list[dict[str, Any]] = []
    for pair in sorted(set(nef_constraints) | set(star_constraints)):
        nef_constraint = nef_constraints.get(pair)
        star_constraint = star_constraints.get(pair)
        if nef_constraint is None:
            discrepancy_type = "star_only"
        elif star_constraint is None:
            discrepancy_type = "nef_only"
        elif abs(
            float(nef_constraint["max_distance"])
            - float(star_constraint["max_distance"])
        ) > BOUND_TOLERANCE:
            discrepancy_type = "different_bound"
        else:
            continue
        source_group_ids = [
            str(group)
            for constraint in (nef_constraint, star_constraint)
            if constraint is not None
            for group in constraint.get("source_groups", [])
        ]
        matched_keys = {_normalize_group_id(group) for group in source_group_ids}
        nef = _format_evidence(
            nef_dict, nef_constraint, matched_keys, nef_groups, nef_rejections
        )
        star = _format_evidence(
            star_dict, star_constraint, matched_keys, star_groups, star_rejections
        )
        classification, rationale_code, rationale = _classify(
            matched_keys, nef_groups, star_groups, nef, star
        )
        pair_label = _pair_label(pair)
        audit_id = f"{case_id.upper()}|{discrepancy_type}|{pair_label}"
        nef_bound = nef["final_bound"]
        star_bound = star["final_bound"]
        row: dict[str, Any] = {
            "audit_id": audit_id,
            "case_id": case_id.upper(),
            "discrepancy_type": discrepancy_type,
            "classification": classification,
            "rationale_code": rationale_code,
            "rationale": rationale,
            "projected_heavy_pair": pair_label,
            "nef_final_bound_angstrom": nef_bound,
            "star_final_bound_angstrom": star_bound,
            "absolute_bound_delta_angstrom": (
                abs(nef_bound - star_bound)
                if nef_bound is not None and star_bound is not None
                else None
            ),
        }
        for prefix, evidence in (("nef", nef), ("star", star)):
            for source, destination in (
                ("source_restraint_ids", "source_restraint_ids"),
                ("source_row_ids", "source_row_ids"),
                ("atom_expressions", "atom_expressions"),
                ("canonical_expansions", "canonical_expansions"),
                ("physical_proton_sets", "physical_proton_sets"),
                ("pseudoatom_handling", "pseudoatom_handling"),
                ("pair_counts_N", "pair_counts_N"),
                ("averaging_policies", "averaging_policies"),
                ("source_upper_bounds", "source_upper_bounds_angstrom"),
                ("projected_bounds", "projected_bounds_angstrom"),
                ("rejection_reasons", "rejection_reasons"),
            ):
                row[f"{prefix}_{destination}"] = evidence[source]
        rows.append(row)
    return rows


def audit_digest(rows: list[dict[str, Any]]) -> str:
    reviewed = [
        {field: row.get(field) for field in AUDIT_FIELDS}
        for row in sorted(rows, key=lambda item: item["audit_id"])
    ]
    payload = json.dumps(reviewed, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def audit_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "row_count": len(rows),
        "digest_sha256": audit_digest(rows),
        "discrepancy_type_counts": dict(Counter(row["discrepancy_type"] for row in rows)),
        "classification_counts": dict(Counter(row["classification"] for row in rows)),
        "case_counts": dict(Counter(row["case_id"] for row in rows)),
        "unresolved_count": sum(row["classification"] == "unresolved" for row in rows),
        "parser_projection_bug_count": sum(
            row["classification"] == "parser_projection_bug" for row in rows
        ),
    }


def write_audit_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: (
                        json.dumps(row.get(field), separators=(",", ":"), sort_keys=True)
                        if isinstance(row.get(field), (list, dict))
                        else row.get(field)
                    )
                    for field in AUDIT_FIELDS
                }
            )
