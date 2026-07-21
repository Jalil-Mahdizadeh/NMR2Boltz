#!/usr/bin/env python3
"""Compare nmr2boltz output with a deposited, hydrogen-complete NMR ensemble.

The script intentionally uses the conversion report as the source of truth for
restraint logic and topology resolution. It writes machine-readable CSV/JSON
evidence suitable for an independent audit of one real-data conversion.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import gemmi

from nmr2boltz.topology import TopologyLibrary, TopologyResolutionError


AtomKey = tuple[str, int, str]
Point = tuple[float, float, float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--cif", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    return parser.parse_args()


def _atom_priority(atom: gemmi.Atom) -> int:
    altloc = str(atom.altloc)
    return 0 if altloc in {"\x00", " ", ""} else (1 if altloc == "A" else 2)


def load_structure(
    path: Path,
    *,
    heavy_atoms_only: bool = False,
) -> tuple[list[str], list[dict[AtomKey, Point]], list[dict[tuple[str, int], str]]]:
    structure = gemmi.read_structure(str(path))
    model_ids: list[str] = []
    models: list[dict[AtomKey, Point]] = []
    residues: list[dict[tuple[str, int], str]] = []
    for position, model in enumerate(structure, start=1):
        model_ids.append(str(model.num or position))
        atoms: dict[AtomKey, Point] = {}
        priorities: dict[AtomKey, int] = {}
        names: dict[tuple[str, int], str] = {}
        for chain in model:
            for residue in chain:
                sequence_number = int(residue.seqid.num)
                residue_key = (chain.name, sequence_number)
                names.setdefault(residue_key, residue.name.upper())
                for atom in residue:
                    if heavy_atoms_only and atom.element.is_hydrogen:
                        continue
                    key = (chain.name, sequence_number, atom.name.strip())
                    priority = _atom_priority(atom)
                    if key in atoms and priority >= priorities[key]:
                        continue
                    atoms[key] = (float(atom.pos.x), float(atom.pos.y), float(atom.pos.z))
                    priorities[key] = priority
        models.append(atoms)
        residues.append(names)
    return model_ids, models, residues


def _sequence_value(record: Any, key: str) -> Any:
    return record[key] if isinstance(record, dict) else getattr(record, key)


def _component_score(target: str, observed: str) -> int:
    if target.upper() == observed.upper():
        return 3
    target_info = gemmi.find_tabulated_residue(target)
    observed_info = gemmi.find_tabulated_residue(observed)
    target_symbol = str(target_info.one_letter_code).strip().upper()
    observed_symbol = str(observed_info.one_letter_code).strip().upper()
    if (
        target_symbol
        and observed_symbol
        and target_symbol == observed_symbol
        and target_info.is_amino_acid() == observed_info.is_amino_acid()
        and target_info.is_nucleic_acid() == observed_info.is_nucleic_acid()
    ):
        return 2
    return -2


def _align_components(target: list[str], observed: list[str]) -> list[tuple[int, int]]:
    gap = -2
    scores = [[0] * (len(observed) + 1) for _ in range(len(target) + 1)]
    trace: list[list[str | None]] = [
        [None] * (len(observed) + 1) for _ in range(len(target) + 1)
    ]
    for index in range(1, len(target) + 1):
        scores[index][0] = index * gap
        trace[index][0] = "up"
    for index in range(1, len(observed) + 1):
        scores[0][index] = index * gap
        trace[0][index] = "left"
    for row in range(1, len(target) + 1):
        for column in range(1, len(observed) + 1):
            choices = [
                (
                    scores[row - 1][column - 1]
                    + _component_score(target[row - 1], observed[column - 1]),
                    "diagonal",
                ),
                (scores[row - 1][column] + gap, "up"),
                (scores[row][column - 1] + gap, "left"),
            ]
            scores[row][column], trace[row][column] = max(
                choices, key=lambda item: (item[0], item[1] == "diagonal")
            )

    aligned: list[tuple[int, int]] = []
    row = len(target)
    column = len(observed)
    while row or column:
        move = trace[row][column]
        if move == "diagonal":
            aligned.append((row - 1, column - 1))
            row -= 1
            column -= 1
        elif move == "up":
            row -= 1
        else:
            column -= 1
    aligned.reverse()
    return aligned


def align_structure_to_sequence_map(
    sequence_map: list[Any],
    models: list[dict[AtomKey, Point]],
    residue_models: list[dict[tuple[str, int], str]],
) -> tuple[list[dict[AtomKey, Point]], list[dict[str, Any]]]:
    """Re-key author-numbered coordinates onto Boltz one-based sequence indices."""
    if not residue_models:
        return models, []
    target_by_chain: dict[str, list[Any]] = defaultdict(list)
    for record in sequence_map:
        target_by_chain[str(_sequence_value(record, "boltz_chain"))].append(record)
    for records in target_by_chain.values():
        records.sort(key=lambda item: int(_sequence_value(item, "boltz_residue_index")))

    observed_by_chain: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for (chain, sequence_number), residue_name in residue_models[0].items():
        if residue_name.upper() not in {"HOH", "WAT", "DOD"}:
            observed_by_chain[chain].append((sequence_number, residue_name))

    residue_mapping: dict[tuple[str, int], tuple[str, int]] = {}
    alignment: list[dict[str, Any]] = []
    used_observed_chains: set[str] = set()
    for target_chain, records in target_by_chain.items():
        candidates = [chain for chain in observed_by_chain if chain not in used_observed_chains]
        if not candidates:
            raise ValueError(f"No unused coordinate chain is available for Boltz chain {target_chain}.")
        target_names = [str(_sequence_value(record, "residue_name")) for record in records]
        ranked: list[tuple[tuple[int, int, int, int], str, list[tuple[int, int]]]] = []
        for observed_chain in candidates:
            observed_names = [name for _number, name in observed_by_chain[observed_chain]]
            pairs = _align_components(target_names, observed_names)
            exact = sum(
                target_names[left].upper() == observed_names[right].upper()
                for left, right in pairs
            )
            compatible = sum(
                _component_score(target_names[left], observed_names[right]) > 0
                for left, right in pairs
            )
            ranked.append(
                (
                    (
                        compatible,
                        exact,
                        int(observed_chain == target_chain),
                        -abs(len(target_names) - len(observed_names)),
                    ),
                    observed_chain,
                    pairs,
                )
            )
        _score, observed_chain, pairs = max(ranked, key=lambda item: item[0])
        used_observed_chains.add(observed_chain)
        observed = observed_by_chain[observed_chain]
        exact = 0
        compatible = 0
        for target_index, observed_index in pairs:
            record = records[target_index]
            observed_number, observed_name = observed[observed_index]
            residue_mapping[(observed_chain, observed_number)] = (
                target_chain,
                int(_sequence_value(record, "boltz_residue_index")),
            )
            target_name = str(_sequence_value(record, "residue_name"))
            exact += target_name.upper() == observed_name.upper()
            compatible += _component_score(target_name, observed_name) > 0
        alignment.append(
            {
                "boltz_chain": target_chain,
                "coordinate_chain": observed_chain,
                "target_residues": len(records),
                "coordinate_residues": len(observed),
                "aligned_residues": len(pairs),
                "exact_component_matches": exact,
                "compatible_component_matches": compatible,
                "unmapped_target_residues": len(records) - len(pairs),
                "unmapped_coordinate_residues": len(observed) - len(pairs),
            }
        )

    aligned_models: list[dict[AtomKey, Point]] = []
    for atoms in models:
        aligned_models.append(
            {
                (*residue_mapping[(chain, sequence_number)], atom_name): point
                for (chain, sequence_number, atom_name), point in atoms.items()
                if (chain, sequence_number) in residue_mapping
            }
        )
    return aligned_models, alignment


def distance(point1: Point, point2: Point) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(point1, point2)))


def quantile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = fraction * (len(ordered) - 1)
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - index) + ordered[upper] * (index - lower)


def descriptive(values: Iterable[float]) -> dict[str, float | int | None]:
    materialized = list(values)
    return {
        "count": len(materialized),
        "minimum": min(materialized) if materialized else None,
        "median": statistics.median(materialized) if materialized else None,
        "mean": statistics.fmean(materialized) if materialized else None,
        "p95": quantile(materialized, 0.95),
        "maximum": max(materialized) if materialized else None,
    }


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    mean_x = statistics.fmean(xs)
    mean_y = statistics.fmean(ys)
    numerator = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    denominator = math.sqrt(
        sum((x - mean_x) ** 2 for x in xs) * sum((y - mean_y) ** 2 for y in ys)
    )
    return numerator / denominator if denominator else None


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def atom_key(payload: dict[str, Any]) -> AtomKey:
    return (
        str(payload["chain"]),
        int(payload["residue_index"]),
        str(payload["atom_name"]),
    )


class SequenceLookup:
    def __init__(self, sequence_map: list[dict[str, Any]]) -> None:
        self.exact: dict[tuple[str, str, str], tuple[str, int]] = {}
        self.loose: dict[tuple[str, str], set[tuple[str, int]]] = defaultdict(set)
        for record in sequence_map:
            destination = (str(record["boltz_chain"]), int(record["boltz_residue_index"]))
            keys = [
                (
                    str(record.get("source_chain") or ""),
                    str(record.get("source_sequence_code") or ""),
                    str(record.get("residue_name") or "").upper(),
                )
            ]
            keys.extend(tuple(str(part) for part in alias) for alias in record.get("aliases", []))
            for chain, sequence, residue in keys:
                normalized = (chain, sequence, residue.upper())
                self.exact[normalized] = destination
                self.loose[(chain, sequence)].add(destination)

    def resolve(self, endpoint: dict[str, Any]) -> tuple[str, int]:
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
            if candidate in self.exact:
                return self.exact[candidate]
        for chain, sequence, _residue in candidates:
            destinations = self.loose.get((chain, sequence), set())
            if len(destinations) == 1:
                return next(iter(destinations))
        raise KeyError(f"unmapped endpoint: {endpoint}")


def effective_distance(distances: list[float], policy: str) -> float:
    if not distances:
        raise ValueError("an atom-set branch must contain at least one distance")
    if min(distances) <= 1e-12:
        return 0.0
    if policy == "hard-or":
        return min(distances)
    inverse_sixth_sum = sum(value**-6 for value in distances)
    if policy == "mean-r6":
        inverse_sixth_sum /= len(distances)
    if policy not in {"sum-r6", "mean-r6"}:
        raise ValueError(f"unsupported averaging policy: {policy}")
    return inverse_sixth_sum ** (-1.0 / 6.0)


def endpoint_choices(
    endpoint: dict[str, Any],
    sequence_lookup: SequenceLookup,
    topology: TopologyLibrary,
    pseudoatom_policy: str,
) -> tuple[tuple[str, int], list[Any]]:
    destination = sequence_lookup.resolve(endpoint)
    residue_name = str(
        endpoint.get("canonical_residue_name") or endpoint.get("residue_name") or ""
    ).upper()
    expression = endpoint.get("atom_expression")
    if not expression:
        raise ValueError(f"endpoint has no atom expression: {endpoint}")
    choices = topology.resolve_expression(
        residue_name,
        str(expression),
        canonical_hint=endpoint.get("canonical_atom_hint"),
        pseudoatom_policy=pseudoatom_policy,
    )
    return destination, choices


def evaluate_source_groups(
    report: dict[str, Any],
    model_ids: list[str],
    models: list[dict[AtomKey, Point]],
    tolerance: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, str]]]:
    topology = TopologyLibrary()
    sequence_lookup = SequenceLookup(report["sequence_map"])
    policy = str(report["settings"]["averaging_policy"])
    pseudoatom_policy = str(report["settings"]["pseudoatom_policy"])
    detail_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    states: dict[str, dict[str, str]] = defaultdict(dict)

    for group in report["source_restraint_groups"]:
        group_id = f'{group["list_name"]}:{group["restraint_id"]}'
        prepared: list[dict[str, Any]] = []
        preparation_errors: list[str] = []
        for alternative_index, alternative in enumerate(group["alternatives"], start=1):
            try:
                destination1, choices1 = endpoint_choices(
                    alternative["endpoint1"], sequence_lookup, topology, pseudoatom_policy
                )
                destination2, choices2 = endpoint_choices(
                    alternative["endpoint2"], sequence_lookup, topology, pseudoatom_policy
                )
                upper_bound = alternative.get("upper_bound")
                if upper_bound is None:
                    raise ValueError("missing upper bound")
                for choice1 in choices1:
                    for choice2 in choices2:
                        prepared.append(
                            {
                                "alternative_index": alternative_index,
                                "destination1": destination1,
                                "destination2": destination2,
                                "atoms1": [atom.atom_name for atom in choice1.atoms],
                                "atoms2": [atom.atom_name for atom in choice2.atoms],
                                "set_semantics1": "atom-set" in choice1.semantics,
                                "set_semantics2": "atom-set" in choice2.semantics,
                                "upper_bound": float(upper_bound),
                                "label": (
                                    f'{alternative["endpoint1"].get("atom_expression")}:'
                                    f'{choice1.assignment_key or "set"}--'
                                    f'{alternative["endpoint2"].get("atom_expression")}:'
                                    f'{choice2.assignment_key or "set"}'
                                ),
                            }
                        )
            except (KeyError, ValueError, TopologyResolutionError) as exc:
                preparation_errors.append(f"alternative {alternative_index}: {exc}")

        group_model_rows: list[dict[str, Any]] = []
        for model_id, atoms in zip(model_ids, models):
            evaluated: list[dict[str, Any]] = []
            coordinate_errors: list[str] = []
            absent_set_atoms: list[str] = []
            for branch in prepared:
                pair_distances: list[float] = []
                missing: list[str] = []
                atoms1 = [
                    atom1
                    for atom1 in branch["atoms1"]
                    if (*branch["destination1"], atom1) in atoms
                ]
                atoms2 = [
                    atom2
                    for atom2 in branch["atoms2"]
                    if (*branch["destination2"], atom2) in atoms
                ]
                missing1 = sorted(set(branch["atoms1"]) - set(atoms1))
                missing2 = sorted(set(branch["atoms2"]) - set(atoms2))
                if missing1 and branch["set_semantics1"] and atoms1:
                    absent_set_atoms.extend(
                        ":".join(map(str, (*branch["destination1"], name)))
                        for name in missing1
                    )
                if missing2 and branch["set_semantics2"] and atoms2:
                    absent_set_atoms.extend(
                        ":".join(map(str, (*branch["destination2"], name)))
                        for name in missing2
                    )
                if missing1 and (not branch["set_semantics1"] or not atoms1):
                    missing.extend(
                        ":".join(map(str, (*branch["destination1"], name)))
                        for name in missing1
                    )
                if missing2 and (not branch["set_semantics2"] or not atoms2):
                    missing.extend(
                        ":".join(map(str, (*branch["destination2"], name)))
                        for name in missing2
                    )
                if missing:
                    coordinate_errors.append(
                        f'{branch["label"]}: missing {",".join(sorted(set(missing)))}'
                    )
                    continue
                for atom1 in atoms1:
                    key1 = (*branch["destination1"], atom1)
                    for atom2 in atoms2:
                        key2 = (*branch["destination2"], atom2)
                        pair_distances.append(distance(atoms[key1], atoms[key2]))
                value = effective_distance(pair_distances, policy)
                evaluated.append(
                    {
                        **branch,
                        "effective_distance": value,
                        "margin": value - branch["upper_bound"],
                    }
                )

            winner = min(evaluated, key=lambda row: row["margin"]) if evaluated else None
            complete = not preparation_errors and len(evaluated) == len(prepared)
            if winner is not None and winner["margin"] <= tolerance:
                state = "satisfied"
            elif winner is not None and complete:
                state = "violated"
            else:
                state = "indeterminate"
            states[group_id][model_id] = state
            row = {
                "group_id": group_id,
                "list_name": group["list_name"],
                "restraint_id": group["restraint_id"],
                "origin": group.get("origin"),
                "model": model_id,
                "state": state,
                "best_effective_distance": winner["effective_distance"] if winner else None,
                "winning_upper_bound": winner["upper_bound"] if winner else None,
                "margin": winner["margin"] if winner else None,
                "winning_branch": winner["label"] if winner else None,
                "branches_prepared": len(prepared),
                "branches_evaluated": len(evaluated),
                "preparation_error_count": len(preparation_errors),
                "coordinate_error_count": len(coordinate_errors),
                "absent_set_atom_count": len(set(absent_set_atoms)),
                "errors": " | ".join(preparation_errors + coordinate_errors),
            }
            detail_rows.append(row)
            group_model_rows.append(row)

        state_counts = Counter(row["state"] for row in group_model_rows)
        resolved_margins = [
            float(row["margin"]) for row in group_model_rows if row["margin"] is not None
        ]
        summary_rows.append(
            {
                "group_id": group_id,
                "origin": group.get("origin"),
                "alternative_count": len(group["alternatives"]),
                "models_satisfied": state_counts["satisfied"],
                "models_violated": state_counts["violated"],
                "models_indeterminate": state_counts["indeterminate"],
                "minimum_margin": min(resolved_margins) if resolved_margins else None,
                "mean_margin": statistics.fmean(resolved_margins) if resolved_margins else None,
                "maximum_margin": max(resolved_margins) if resolved_margins else None,
            }
        )
    return detail_rows, summary_rows, states


def separation_class(atom1: AtomKey, atom2: AtomKey) -> str:
    if atom1[0] != atom2[0]:
        return "interchain"
    separation = abs(atom1[1] - atom2[1])
    if separation == 0:
        return "intraresidue"
    if separation <= 4:
        return "sequential"
    if separation <= 23:
        return "medium_range"
    return "long_range"


def evaluate_heavy_constraints(
    report: dict[str, Any],
    model_ids: list[str],
    models: list[dict[AtomKey, Point]],
    tolerance: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    detail_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for index, constraint in enumerate(report["emitted_constraints"], start=1):
        key1 = atom_key(constraint["atom1"])
        key2 = atom_key(constraint["atom2"])
        upper_bound = float(constraint["max_distance"])
        constraint_id = f"C{index:05d}"
        model_rows: list[dict[str, Any]] = []
        for model_id, atoms in zip(model_ids, models):
            if key1 not in atoms or key2 not in atoms:
                row = {
                    "constraint_id": constraint_id,
                    "model": model_id,
                    "atom1": ":".join(map(str, key1)),
                    "atom2": ":".join(map(str, key2)),
                    "separation_class": separation_class(key1, key2),
                    "upper_bound": upper_bound,
                    "distance": None,
                    "excess": None,
                    "satisfied": None,
                    "source_groups": ";".join(constraint["source_groups"]),
                }
            else:
                value = distance(atoms[key1], atoms[key2])
                excess = value - upper_bound
                row = {
                    "constraint_id": constraint_id,
                    "model": model_id,
                    "atom1": ":".join(map(str, key1)),
                    "atom2": ":".join(map(str, key2)),
                    "separation_class": separation_class(key1, key2),
                    "upper_bound": upper_bound,
                    "distance": value,
                    "excess": excess,
                    "satisfied": excess <= tolerance,
                    "source_groups": ";".join(constraint["source_groups"]),
                }
            detail_rows.append(row)
            model_rows.append(row)
        resolved = [row for row in model_rows if row["distance"] is not None]
        distances = [float(row["distance"]) for row in resolved]
        excesses = [float(row["excess"]) for row in resolved]
        summary_rows.append(
            {
                "constraint_id": constraint_id,
                "atom1": ":".join(map(str, key1)),
                "atom2": ":".join(map(str, key2)),
                "separation_class": separation_class(key1, key2),
                "upper_bound": upper_bound,
                "models_resolved": len(resolved),
                "models_satisfied": sum(row["satisfied"] is True for row in resolved),
                "models_violated": sum(row["satisfied"] is False for row in resolved),
                "minimum_distance": min(distances) if distances else None,
                "mean_distance": statistics.fmean(distances) if distances else None,
                "maximum_distance": max(distances) if distances else None,
                "maximum_excess": max(excesses) if excesses else None,
                "source_groups": ";".join(constraint["source_groups"]),
            }
        )
    return detail_rows, summary_rows


def coordinate_crosscheck(
    pdb_ids: list[str],
    pdb_models: list[dict[AtomKey, Point]],
    cif_ids: list[str],
    cif_models: list[dict[AtomKey, Point]],
) -> dict[str, Any]:
    maximum_delta = 0.0
    mismatch_examples: list[str] = []
    common_atom_count = 0
    pdb_only_count = 0
    cif_only_count = 0
    if len(pdb_models) != len(cif_models):
        mismatch_examples.append(f"model count PDB={len(pdb_models)}, CIF={len(cif_models)}")
    for index, (pdb_atoms, cif_atoms) in enumerate(zip(pdb_models, cif_models), start=1):
        common = set(pdb_atoms) & set(cif_atoms)
        common_atom_count += len(common)
        pdb_only = set(pdb_atoms) - set(cif_atoms)
        cif_only = set(cif_atoms) - set(pdb_atoms)
        pdb_only_count += len(pdb_only)
        cif_only_count += len(cif_only)
        if (pdb_only or cif_only) and len(mismatch_examples) < 10:
            mismatch_examples.append(
                f"model {index}: PDB-only={len(pdb_only)}, CIF-only={len(cif_only)}"
            )
        for key in common:
            delta = distance(pdb_atoms[key], cif_atoms[key])
            maximum_delta = max(maximum_delta, delta)
            if delta > 1e-6 and len(mismatch_examples) < 10:
                mismatch_examples.append(f"model {index} {key}: coordinate delta {delta:.9g}")
    return {
        "pdb_model_ids": pdb_ids,
        "cif_model_ids": cif_ids,
        "pdb_model_count": len(pdb_models),
        "cif_model_count": len(cif_models),
        "common_atom_instances": common_atom_count,
        "pdb_only_atom_instances": pdb_only_count,
        "cif_only_atom_instances": cif_only_count,
        "maximum_coordinate_delta_angstrom": maximum_delta,
        "mismatch_examples": mismatch_examples,
    }


def implication_audit(
    report: dict[str, Any],
    heavy_rows: list[dict[str, Any]],
    source_states: dict[str, dict[str, str]],
    tolerance: float,
) -> dict[str, Any]:
    constraints = {
        f"C{index:05d}": constraint
        for index, constraint in enumerate(report["emitted_constraints"], start=1)
    }
    antecedent_cases = 0
    failures: list[dict[str, Any]] = []
    for row in heavy_rows:
        constraint = constraints[row["constraint_id"]]
        model_id = row["model"]
        groups = constraint["source_groups"]
        if groups and all(
            source_states.get(group, {}).get(model_id) == "satisfied" for group in groups
        ):
            antecedent_cases += 1
            if row["distance"] is None or float(row["excess"]) > tolerance:
                failures.append(
                    {
                        "constraint_id": row["constraint_id"],
                        "model": model_id,
                        "distance": row["distance"],
                        "upper_bound": row["upper_bound"],
                        "excess": row["excess"],
                        "source_groups": groups,
                    }
                )
    return {
        "antecedent_cases": antecedent_cases,
        "failure_count": len(failures),
        "failures": failures[:50],
    }


def main() -> int:
    args = parse_args()
    if args.tolerance < 0 or not math.isfinite(args.tolerance):
        raise ValueError("--tolerance must be finite and non-negative")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = json.loads(args.report.read_text(encoding="utf-8"))
    pdb_ids, pdb_models, pdb_residues = load_structure(args.pdb)
    cif_ids, cif_models, cif_residues = load_structure(args.cif)
    if not pdb_models:
        raise ValueError("PDB file contains no models")
    if any(names != pdb_residues[0] for names in pdb_residues[1:]):
        raise ValueError("PDB models do not have identical residue identities")
    if any(names != cif_residues[0] for names in cif_residues[1:]):
        raise ValueError("mmCIF models do not have identical residue identities")

    pdb_models, pdb_alignment = align_structure_to_sequence_map(
        report["sequence_map"], pdb_models, pdb_residues
    )
    cif_models, cif_alignment = align_structure_to_sequence_map(
        report["sequence_map"], cif_models, cif_residues
    )
    crosscheck = coordinate_crosscheck(pdb_ids, pdb_models, cif_ids, cif_models)
    source_rows, source_summary, source_states = evaluate_source_groups(
        report, pdb_ids, pdb_models, args.tolerance
    )
    heavy_rows, heavy_summary = evaluate_heavy_constraints(
        report, pdb_ids, pdb_models, args.tolerance
    )
    implication = implication_audit(report, heavy_rows, source_states, args.tolerance)

    write_csv(
        args.output_dir / "source_restraint_distances.csv",
        source_rows,
        [
            "group_id", "list_name", "restraint_id", "origin", "model", "state",
            "best_effective_distance", "winning_upper_bound", "margin", "winning_branch",
            "branches_prepared", "branches_evaluated", "preparation_error_count",
            "coordinate_error_count", "absent_set_atom_count", "errors",
        ],
    )
    write_csv(
        args.output_dir / "source_restraint_summary.csv",
        source_summary,
        [
            "group_id", "origin", "alternative_count", "models_satisfied", "models_violated",
            "models_indeterminate", "minimum_margin", "mean_margin", "maximum_margin",
        ],
    )
    write_csv(
        args.output_dir / "heavy_atom_distances.csv",
        heavy_rows,
        [
            "constraint_id", "model", "atom1", "atom2", "separation_class", "upper_bound",
            "distance", "excess", "satisfied", "source_groups",
        ],
    )
    write_csv(
        args.output_dir / "heavy_atom_summary.csv",
        heavy_summary,
        [
            "constraint_id", "atom1", "atom2", "separation_class", "upper_bound",
            "models_resolved", "models_satisfied", "models_violated", "minimum_distance",
            "mean_distance", "maximum_distance", "maximum_excess", "source_groups",
        ],
    )

    source_state_counts = Counter(row["state"] for row in source_rows)
    heavy_resolved = [row for row in heavy_rows if row["distance"] is not None]
    heavy_excesses = [float(row["excess"]) for row in heavy_resolved]
    heavy_violations = [row for row in heavy_resolved if row["satisfied"] is False]
    per_constraint_violated = [row for row in heavy_summary if row["models_violated"] > 0]
    top_heavy_violations = sorted(
        per_constraint_violated,
        key=lambda row: float(row["maximum_excess"]),
        reverse=True,
    )[:20]
    source_group_violations = [row for row in source_summary if row["models_violated"] > 0]
    top_source_violations = sorted(
        source_group_violations,
        key=lambda row: float(row["maximum_margin"]),
        reverse=True,
    )[:20]
    bounds = [float(row["upper_bound"]) for row in heavy_resolved]
    observed = [float(row["distance"]) for row in heavy_resolved]
    summary = {
        "inputs": {
            "conversion_report": str(args.report),
            "pdb": str(args.pdb),
            "cif": str(args.cif),
            "tolerance_angstrom": args.tolerance,
        },
        "sequence_alignment": {
            "pdb": pdb_alignment,
            "cif": cif_alignment,
        },
        "coordinate_crosscheck": crosscheck,
        "source_restraints": {
            "group_count": len(source_summary),
            "group_model_cases": len(source_rows),
            "state_counts": dict(source_state_counts),
            "groups_with_at_least_one_violated_model": len(source_group_violations),
            "groups_satisfied_in_every_model": sum(
                row["models_satisfied"] == len(pdb_models) for row in source_summary
            ),
            "groups_indeterminate_in_any_model": sum(
                row["models_indeterminate"] > 0 for row in source_summary
            ),
            "top_violations": top_source_violations,
        },
        "heavy_atom_constraints": {
            "constraint_count": len(heavy_summary),
            "constraint_model_cases": len(heavy_rows),
            "resolved_cases": len(heavy_resolved),
            "satisfied_cases": sum(row["satisfied"] is True for row in heavy_resolved),
            "violated_cases": len(heavy_violations),
            "constraints_satisfied_in_every_model": sum(
                row["models_satisfied"] == len(pdb_models) for row in heavy_summary
            ),
            "constraints_with_at_least_one_violated_model": len(per_constraint_violated),
            "constraints_missing_in_any_model": sum(
                row["models_resolved"] != len(pdb_models) for row in heavy_summary
            ),
            "separation_class_counts": dict(
                Counter(row["separation_class"] for row in heavy_summary)
            ),
            "distance_minus_bound_angstrom": descriptive(heavy_excesses),
            "bound_observed_distance_pearson": pearson(bounds, observed),
            "top_violations": top_heavy_violations,
        },
        "projection_implication_audit": implication,
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if implication["failure_count"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
