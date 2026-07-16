from __future__ import annotations

import csv
from decimal import Decimal, ROUND_CEILING
import itertools
import json
import math
import random
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable

import yaml

from .model import (
    AmbiguousGroup,
    BoltzAtom,
    ConversionReport,
    EmittedConstraint,
    ProjectedAlternative,
)


class FlowList(list):
    pass


class BoltzYamlDumper(yaml.SafeDumper):
    pass


def _flow_list_representer(dumper: yaml.SafeDumper, data: FlowList) -> yaml.Node:
    return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)


BoltzYamlDumper.add_representer(FlowList, _flow_list_representer)


def _atom_list(atom: BoltzAtom) -> FlowList:
    return FlowList([atom.chain, atom.residue_index, atom.atom_name])


def _rounded(value: float) -> float:
    return float(Decimal(str(value)).quantize(Decimal("0.000001"), rounding=ROUND_CEILING))


def _constraint_payload(constraints: Iterable[EmittedConstraint]) -> dict[str, Any]:
    return {
        "constraints": [
            {
                "atom_contact": {
                    "atom1": _atom_list(item.atom1),
                    "atom2": _atom_list(item.atom2),
                    "max_distance": _rounded(item.max_distance),
                    "force": True,
                }
            }
            for item in constraints
        ]
    }


def _dump_yaml(payload: Any) -> str:
    return yaml.dump(
        payload,
        Dumper=BoltzYamlDumper,
        sort_keys=False,
        allow_unicode=True,
        default_flow_style=False,
        width=120,
    )


def write_outputs(
    report: ConversionReport,
    output_dir: str | Path,
    *,
    hypothesis_count: int = 0,
    random_seed: int = 0,
) -> list[Path]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    report_path = output / "conversion_report.json"
    report_path.write_text(
        json.dumps(report.to_dict(), indent=2, sort_keys=False, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    written.append(report_path)

    yaml_path = output / "boltz_constraints.yaml"
    yaml_header = (
        "# Conservative heavy-atom constraints safe to use as simultaneous (AND) Boltz contacts.\n"
        "# Merge the `constraints` list into the Boltz input YAML and run Boltz-2 with potentials enabled.\n"
    )
    yaml_path.write_text(
        yaml_header + _dump_yaml(_constraint_payload(report.emitted_constraints)),
        encoding="utf-8",
    )
    written.append(yaml_path)

    text_path = output / "heavy_atom_constraints.txt"
    with text_path.open("w", encoding="utf-8") as handle:
        handle.write("# [heavy-atom X -- heavy-atom Y : conservative upper bound in angstrom]\n")
        for item in report.emitted_constraints:
            handle.write(
                f"[{item.atom1.display()} -- {item.atom2.display()} : {_rounded(item.max_distance):.6f}]\n"
            )
    written.append(text_path)

    tsv_path = output / "heavy_atom_constraints.tsv"
    _write_tsv(
        tsv_path,
        [
            "chain_1",
            "residue_index_1",
            "atom_1",
            "chain_2",
            "residue_index_2",
            "atom_2",
            "max_distance_A",
            "raw_projected_distance_A",
            "source_groups",
            "boltz_adjustment",
        ],
        [
            [
                item.atom1.chain,
                item.atom1.residue_index,
                item.atom1.atom_name,
                item.atom2.chain,
                item.atom2.residue_index,
                item.atom2.atom_name,
                f"{_rounded(item.max_distance):.6f}",
                f"{item.raw_projected_distance:.6f}",
                ";".join(item.source_groups),
                item.boltz_adjustment or "",
            ]
            for item in report.emitted_constraints
        ],
    )
    written.append(tsv_path)

    sequence_path = output / "sequence_map.tsv"
    _write_tsv(
        sequence_path,
        [
            "source_chain",
            "source_sequence_code",
            "source_residue_name",
            "boltz_chain",
            "boltz_residue_index",
            "mapping_source",
            "aliases",
            "warnings",
        ],
        [
            [
                record.source_chain,
                record.source_sequence_code,
                record.residue_name,
                record.boltz_chain,
                record.boltz_residue_index,
                record.source,
                json.dumps(record.aliases, ensure_ascii=False),
                "; ".join(record.warnings),
            ]
            for record in report.sequence_map
        ],
    )
    written.append(sequence_path)

    ambiguity_path = output / "ambiguous_groups.tsv"
    ambiguity_rows: list[list[Any]] = []
    for group in report.ambiguous_groups:
        for index, alternative in enumerate(group.alternatives, start=1):
            ambiguity_rows.append(
                [
                    group.group_id,
                    group.restraint_id,
                    index,
                    alternative.atom1.chain,
                    alternative.atom1.residue_index,
                    alternative.atom1.atom_name,
                    alternative.atom2.chain,
                    alternative.atom2.residue_index,
                    alternative.atom2.atom_name,
                    f"{_rounded(alternative.max_distance):.6f}",
                    f"{alternative.averaging_factor:.8f}",
                    alternative.explicit_pair_count,
                    ";".join(alternative.source_rows),
                    " | ".join(alternative.source_endpoints),
                    group.reason,
                ]
            )
    _write_tsv(
        ambiguity_path,
        [
            "group_id",
            "restraint_id",
            "alternative_number",
            "chain_1",
            "residue_index_1",
            "atom_1",
            "chain_2",
            "residue_index_2",
            "atom_2",
            "max_distance_A",
            "averaging_factor",
            "explicit_pair_count",
            "source_rows",
            "source_endpoints",
            "reason_not_emitted",
        ],
        ambiguity_rows,
    )
    written.append(ambiguity_path)

    rejection_path = output / "rejections.tsv"
    _write_tsv(
        rejection_path,
        ["group_id", "reason", "details", "row_ids", "endpoint"],
        [
            [
                item.group_id,
                item.reason,
                item.details,
                ";".join(item.row_ids),
                item.endpoint or "",
            ]
            for item in report.rejections
        ],
    )
    written.append(rejection_path)

    union_path = output / "proposed_atom_contact_unions.yaml"
    union_payload = {
        "x_nmr2boltz_schema": "proposed-atom-contact-union-v1",
        "compatible_with_current_boltzui_atom_contact_parser": False,
        "note": (
            "Each item is one OR group. Implement as a shared union index in the Boltz distance "
            "potential; do not mark every alternative as a simultaneous token contact."
        ),
        "constraints": [_union_payload(group) for group in report.ambiguous_groups],
    }
    union_path.write_text(
        "# PROPOSED schema for preserving assignment ambiguity; not accepted by the current BoltzUI parser.\n"
        + _dump_yaml(union_payload),
        encoding="utf-8",
    )
    written.append(union_path)

    summary_path = output / "summary.txt"
    summary_path.write_text(_summary_text(report), encoding="utf-8")
    written.append(summary_path)

    hypothesis_directory = output / "hypotheses"
    if hypothesis_count > 0 and report.ambiguous_groups:
        hypothesis_paths = write_hypotheses(
            report,
            hypothesis_directory,
            count=hypothesis_count,
            seed=random_seed,
        )
        written.extend(hypothesis_paths)
    else:
        _clear_hypothesis_outputs(hypothesis_directory)

    return written


def _union_payload(group: AmbiguousGroup) -> dict[str, Any]:
    return {
        "atom_contact_union": {
            "group_id": group.group_id,
            "restraint_id": group.restraint_id,
            "alternatives": [
                {
                    "atom1": _atom_list(alternative.atom1),
                    "atom2": _atom_list(alternative.atom2),
                    "max_distance": _rounded(alternative.max_distance),
                    "source_rows": alternative.source_rows,
                }
                for alternative in group.alternatives
            ],
            "force": True,
        }
    }


def write_hypotheses(
    report: ConversionReport,
    directory: str | Path,
    *,
    count: int,
    seed: int,
) -> list[Path]:
    destination = Path(directory)
    destination.mkdir(parents=True, exist_ok=True)
    _clear_hypothesis_outputs(destination)
    minimum = float(report.settings.get("boltz_min_distance_angstrom", 2.0))
    maximum = float(report.settings.get("boltz_max_distance_angstrom", 20.0))
    choice_sets: list[list[ProjectedAlternative]] = []
    used_groups: list[AmbiguousGroup] = []
    skipped_groups: list[str] = []
    for group in report.ambiguous_groups:
        candidates = [item for item in group.alternatives if item.max_distance <= maximum]
        if not candidates:
            skipped_groups.append(group.group_id)
            continue
        choice_sets.append(candidates)
        used_groups.append(group)
    if not choice_sets:
        return []

    total_combinations = math.prod(len(choices) for choices in choice_sets)
    requested = min(count, total_combinations)
    rng = random.Random(seed)
    selections: list[tuple[int, ...]] = []
    if total_combinations <= 100_000:
        all_indices = list(itertools.product(*(range(len(items)) for items in choice_sets)))
        rng.shuffle(all_indices)
        selections = all_indices[:requested]
    else:
        seen: set[tuple[int, ...]] = set()
        attempts = 0
        max_attempts = max(10_000, requested * 200)
        while len(selections) < requested and attempts < max_attempts:
            candidate = tuple(rng.randrange(len(items)) for items in choice_sets)
            if candidate not in seen:
                seen.add(candidate)
                selections.append(candidate)
            attempts += 1

    written: list[Path] = []
    manifest: dict[str, Any] = {
        "warning": (
            "Each file makes one explicit assignment hypothesis for every included ambiguous OR group. "
            "A hypothesis is not a logically conservative conversion; run multiple hypotheses and validate "
            "against the original NMR restraints."
        ),
        "seed": seed,
        "requested": count,
        "generated": len(selections),
        "total_possible_combinations": total_combinations,
        "skipped_groups_without_boltz_compatible_alternatives": skipped_groups,
        "hypotheses": [],
    }
    for number, selection in enumerate(selections, start=1):
        selected = [choice_sets[group_index][choice_index] for group_index, choice_index in enumerate(selection)]
        constraints = _combine_hypothesis_constraints(
            report.emitted_constraints,
            selected,
            minimum=minimum,
            maximum=maximum,
        )
        name = f"hypothesis_{number:04d}.yaml"
        path = destination / name
        path.write_text(
            "# Assignment hypothesis: one alternative selected from each ambiguous NMR restraint group.\n"
            + _dump_yaml(_constraint_payload(constraints)),
            encoding="utf-8",
        )
        written.append(path)
        manifest["hypotheses"].append(
            {
                "file": name,
                "selections": [
                    {
                        "group_id": used_groups[index].group_id,
                        "alternative_index": selection[index] + 1,
                        "atom1": asdict_atom(selected[index].atom1),
                        "atom2": asdict_atom(selected[index].atom2),
                        "max_distance": selected[index].max_distance,
                    }
                    for index in range(len(selection))
                ],
            }
        )
    manifest_path = destination / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    written.append(manifest_path)
    return written


def _clear_hypothesis_outputs(directory: Path) -> None:
    if not directory.is_dir():
        return
    for path in directory.glob("hypothesis_*.yaml"):
        path.unlink()
    (directory / "manifest.json").unlink(missing_ok=True)


def _combine_hypothesis_constraints(
    safe: list[EmittedConstraint],
    selected: list[ProjectedAlternative],
    *,
    minimum: float,
    maximum: float,
) -> list[EmittedConstraint]:
    merged: dict[tuple[BoltzAtom, BoltzAtom], EmittedConstraint] = {
        item.pair_key: EmittedConstraint(
            atom1=item.atom1,
            atom2=item.atom2,
            max_distance=item.max_distance,
            source_groups=list(item.source_groups),
            raw_projected_distance=item.raw_projected_distance,
            boltz_adjustment=item.boltz_adjustment,
            provenance=list(item.provenance),
        )
        for item in safe
    }
    for alternative in selected:
        raw = alternative.max_distance
        if raw > maximum:
            continue
        adjusted = max(minimum, raw)
        pair = alternative.pair_key
        if pair in merged:
            current = merged[pair]
            if adjusted < current.max_distance:
                current.max_distance = adjusted
                current.raw_projected_distance = min(current.raw_projected_distance, raw)
            current.source_groups = list(dict.fromkeys(current.source_groups + [alternative.group_id]))
        else:
            merged[pair] = EmittedConstraint(
                atom1=pair[0],
                atom2=pair[1],
                max_distance=adjusted,
                source_groups=[alternative.group_id],
                raw_projected_distance=raw,
                boltz_adjustment=(
                    f"raised to Boltz minimum {minimum:g} A" if raw < minimum else None
                ),
                provenance=[],
            )
    return [merged[key] for key in sorted(merged)]


def asdict_atom(atom: BoltzAtom) -> dict[str, Any]:
    return {"chain": atom.chain, "residue_index": atom.residue_index, "atom_name": atom.atom_name}


def _write_tsv(path: Path, header: list[str], rows: Iterable[Iterable[Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _summary_text(report: ConversionReport) -> str:
    stats = report.statistics
    lines = [
        "NMR-to-Boltz heavy-atom projection summary",
        "=" * 46,
        f"Input: {report.input_file}",
        f"Detected format: {report.detected_format}",
        f"Restraint groups read: {stats.get('restraint_groups_read', 0)}",
        f"Source alternatives read: {stats.get('source_alternatives_read', 0)}",
        f"Safe groups before pair deduplication: {stats.get('safe_groups_before_pair_deduplication', 0)}",
        f"Unique Boltz constraints emitted: {stats.get('emitted_unique_heavy_atom_constraints', 0)}",
        f"Ambiguous OR groups quarantined: {stats.get('ambiguous_or_groups_not_emitted', 0)}",
        f"Rejection records: {stats.get('rejection_records', 0)}",
        "",
        "Important:",
        "- boltz_constraints.yaml contains only constraints that may be imposed simultaneously.",
        "- ambiguous_groups.tsv preserves true OR alternatives; adding all of them would overconstrain the model.",
        "- proposed_atom_contact_unions.yaml is a proposed schema, not accepted by the current BoltzUI parser.",
        "- conversion_report.json is the lossless audit/provenance record.",
        "",
    ]
    if report.warnings:
        lines.append("Global warnings:")
        lines.extend(f"- {warning}" for warning in report.warnings)
        lines.append("")
    return "\n".join(lines)
