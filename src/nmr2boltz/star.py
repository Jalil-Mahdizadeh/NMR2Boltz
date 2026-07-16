from __future__ import annotations

import bz2
import csv
import gzip
import lzma
import json
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import pynmrstar

from .model import (
    Endpoint,
    RawAlternative,
    RestraintGroup,
    SequenceRecord,
    as_float,
    clean,
)
from .topology import ComponentTopology


class StarDataError(ValueError):
    pass


def _parse_entry_file(path: Path) -> Any:
    def from_file(filename: str) -> Any:
        try:
            return pynmrstar.Entry.from_file(filename, raise_parse_warnings=False)
        except TypeError as exc:
            if "raise_parse_warnings" not in str(exc):
                raise
            return pynmrstar.Entry.from_file(filename)

    suffix = path.suffix.lower()
    if suffix not in {".gz", ".bz2", ".xz", ".lzma"}:
        return from_file(str(path))
    opener = {".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open, ".lzma": lzma.open}[suffix]
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    try:
        try:
            return pynmrstar.Entry.from_string(text, raise_parse_warnings=False)
        except TypeError as exc:
            if "raise_parse_warnings" not in str(exc):
                raise
            return pynmrstar.Entry.from_string(text)
    except AttributeError:
        import tempfile

        with tempfile.NamedTemporaryFile("w", suffix=".str", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            return from_file(handle.name)


def _category(loop: Any) -> str:
    return str(getattr(loop, "category", "")).strip().lstrip("_").lower()


def _tag_leaf(tag: Any) -> str:
    text = str(tag).strip().lstrip("_")
    return text.split(".")[-1].lower()


def loop_rows(loop: Any) -> list[dict[str, Any]]:
    tags = [_tag_leaf(tag) for tag in getattr(loop, "tags", [])]
    rows: list[dict[str, Any]] = []
    for raw in getattr(loop, "data", []):
        rows.append({tag: raw[index] if index < len(raw) else None for index, tag in enumerate(tags)})
    return rows


def saveframe_tags(saveframe: Any) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for tag in getattr(saveframe, "tags", []):
        if isinstance(tag, (list, tuple)) and len(tag) >= 2:
            result[_tag_leaf(tag[0])] = tag[1]
    return result


def pick(row: dict[str, Any], *names: str) -> Any:
    for name in names:
        value = clean(row.get(name.lower()))
        if value is not None:
            return value
    return None


def _natural_sequence_key(value: str) -> tuple[int, int, str]:
    match = re.fullmatch(r"([+-]?\d+)(.*)", value)
    if match:
        return (0, int(match.group(1)), match.group(2))
    return (1, 0, value)


@dataclass
class SequenceResolver:
    records: list[SequenceRecord] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    _exact: dict[tuple[str, str, str], SequenceRecord] = field(default_factory=dict)
    _loose: dict[tuple[str, str], list[SequenceRecord]] = field(default_factory=lambda: defaultdict(list))

    def add(self, record: SequenceRecord, override: bool = False) -> None:
        aliases = list(record.aliases)
        aliases.append((record.source_chain, record.source_sequence_code, record.residue_name))
        dedup_aliases: list[tuple[str, str, str]] = []
        for chain, seq, residue in aliases:
            key = (str(chain), str(seq), str(residue).upper())
            if key not in dedup_aliases:
                dedup_aliases.append(key)

        if override:
            replaced: list[SequenceRecord] = []
            for key in dedup_aliases:
                existing = self._exact.get(key)
                if existing is not None and all(existing is not item for item in replaced):
                    replaced.append(existing)
            for existing in replaced:
                self.records = [item for item in self.records if item is not existing]
                for key, value in list(self._exact.items()):
                    if value is existing:
                        del self._exact[key]
                for key, values in list(self._loose.items()):
                    retained = [item for item in values if item is not existing]
                    if retained:
                        self._loose[key] = retained
                    else:
                        del self._loose[key]

        for key in dedup_aliases:
            chain, seq, _ = key
            if override or key not in self._exact:
                self._exact[key] = record
            loose_key = (chain, seq)
            if record not in self._loose[loose_key]:
                self._loose[loose_key].append(record)
        record.aliases = dedup_aliases
        if record not in self.records:
            self.records.append(record)

    def resolve(self, endpoint: Endpoint) -> tuple[SequenceRecord, list[str]]:
        warnings: list[str] = []
        residue_candidates = [
            endpoint.residue_name,
            endpoint.canonical_residue_name,
        ]
        locators = [
            (endpoint.chain_code, endpoint.sequence_code),
            (endpoint.canonical_chain_code, endpoint.canonical_sequence_code),
            (endpoint.chain_code, endpoint.canonical_sequence_code),
            (endpoint.canonical_chain_code, endpoint.sequence_code),
        ]
        for chain, sequence in locators:
            if chain is None or sequence is None:
                continue
            for residue in residue_candidates:
                if residue is None:
                    continue
                key = (chain, sequence, residue.upper())
                if key in self._exact:
                    return self._exact[key], warnings
        for chain, sequence in locators:
            if chain is None or sequence is None:
                continue
            candidates = self._loose.get((chain, sequence), [])
            if len(candidates) == 1:
                record = candidates[0]
                if endpoint.residue_name and endpoint.residue_name.upper() != record.residue_name.upper():
                    warnings.append(
                        f"residue-name mismatch at {chain}:{sequence}: restraint has "
                        f"{endpoint.residue_name}, sequence map has {record.residue_name}"
                    )
                return record, warnings
            if len(candidates) > 1:
                raise StarDataError(
                    f"Ambiguous residue mapping for endpoint {endpoint.display()}: "
                    f"{len(candidates)} sequence records match {chain}:{sequence}."
                )
        # A missing chain can be resolved only if the sequence code is globally unique.
        sequences = [value for value in {endpoint.sequence_code, endpoint.canonical_sequence_code} if value]
        for sequence in sequences:
            candidates = [record for record in self.records if record.source_sequence_code == sequence]
            if len(candidates) == 1:
                warnings.append(
                    f"missing/unresolved chain for {endpoint.display()} inferred as {candidates[0].source_chain}"
                )
                return candidates[0], warnings
        raise StarDataError(f"No Boltz residue mapping for endpoint {endpoint.display()}.")


@dataclass
class ParsedStarDocument:
    entry: Any
    detected_format: str
    sequence_resolver: SequenceResolver
    restraint_groups: list[RestraintGroup]
    embedded_topologies: list[ComponentTopology]
    warnings: list[str]


def parse_star_document(
    path: str | Path,
    *,
    format_hint: str = "auto",
    missing_upper_policy: str = "reject",
    residue_map_path: str | Path | None = None,
    chain_map: dict[str, str] | None = None,
    allow_inferred_sequence_map: bool = False,
    origins: set[str] | None = None,
) -> ParsedStarDocument:
    allowed_missing_upper_policies = {
        "reject",
        "upper-linear",
        "target-plus-uncertainty",
        "target",
    }
    if missing_upper_policy not in allowed_missing_upper_policies:
        raise StarDataError(
            f"Unsupported missing upper-bound policy: {missing_upper_policy!r}"
        )
    input_path = Path(path)
    try:
        entry = _parse_entry_file(input_path)
    except Exception as exc:
        raise StarDataError(f"Unable to parse STAR/NEF file {input_path}: {exc}") from exc

    categories = {_category(loop) for saveframe in entry for loop in saveframe}
    if format_hint == "auto":
        if "nef_distance_restraint" in categories:
            detected = "nef"
        elif "gen_dist_constraint" in categories:
            detected = "nmr-star"
        else:
            raise StarDataError(
                "No _nef_distance_restraint or _Gen_dist_constraint loop was found."
            )
    elif format_hint in {"nef", "nmr-star", "nmrstar"}:
        detected = "nmr-star" if format_hint in {"nmr-star", "nmrstar"} else "nef"
    else:
        raise StarDataError(f"Unsupported format hint: {format_hint}")

    warnings: list[str] = []
    groups = extract_restraint_groups(
        entry,
        detected,
        missing_upper_policy=missing_upper_policy,
        origins=origins,
        warnings=warnings,
    )
    resolver = extract_sequence_map(entry, detected, chain_map=chain_map or {})
    warnings.extend(resolver.warnings)
    if not resolver.records and allow_inferred_sequence_map:
        resolver = infer_sequence_map(groups, chain_map=chain_map or {})
        warnings.extend(resolver.warnings)
    if residue_map_path is not None:
        overlay_residue_map(resolver, residue_map_path, chain_map=chain_map or {})
    if not resolver.records:
        warnings.append(
            "No sequence mapping was found. Supply --residue-map, or use "
            "--allow-inferred-sequence-map after verifying numbering manually."
        )

    embedded_topologies = extract_embedded_topologies(entry)
    return ParsedStarDocument(
        entry=entry,
        detected_format=detected,
        sequence_resolver=resolver,
        restraint_groups=groups,
        embedded_topologies=embedded_topologies,
        warnings=warnings,
    )


def extract_restraint_groups(
    entry: Any,
    detected_format: str,
    *,
    missing_upper_policy: str,
    origins: set[str] | None,
    warnings: list[str],
) -> list[RestraintGroup]:
    grouped: dict[tuple[str, str], list[RawAlternative]] = defaultdict(list)
    group_origins: dict[tuple[str, str], str | None] = {}
    group_warnings: dict[tuple[str, str], list[str]] = defaultdict(list)

    target_category = "nef_distance_restraint" if detected_format == "nef" else "gen_dist_constraint"
    for sf_index, saveframe in enumerate(entry):
        sf_name = str(getattr(saveframe, "name", f"saveframe_{sf_index + 1}"))
        sf_tags = saveframe_tags(saveframe)
        if detected_format == "nef":
            origin = clean(sf_tags.get("restraint_origin"))
        else:
            origin = clean(sf_tags.get("constraint_type")) or clean(sf_tags.get("restraint_origin"))
        normalized_origin = (origin or "unknown").lower()
        if origins is not None and normalized_origin not in origins and "all" not in origins:
            continue

        for loop in saveframe:
            if _category(loop) != target_category:
                continue
            for row_number, row in enumerate(loop_rows(loop), start=1):
                if detected_format == "nef":
                    alternative = _parse_nef_row(
                        row,
                        sf_name,
                        row_number,
                        origin,
                        missing_upper_policy,
                    )
                else:
                    alternative = _parse_nmrstar_row(
                        row,
                        sf_name,
                        row_number,
                        origin,
                        missing_upper_policy,
                    )
                key = (sf_name, alternative.restraint_id)
                grouped[key].append(alternative)
                group_origins[key] = origin

    result: list[RestraintGroup] = []
    for (list_name, restraint_id), alternatives in grouped.items():
        deduplicated: dict[tuple[Any, ...], RawAlternative] = {}
        for alternative in alternatives:
            key = alternative.dedup_key()
            if key in deduplicated:
                existing = deduplicated[key]
                existing.row_ids.extend(row for row in alternative.row_ids if row not in existing.row_ids)
                existing.warnings.append(
                    "duplicate canonical expansion row collapsed to its author-level atom expression"
                )
            else:
                deduplicated[key] = alternative
        unique = list(deduplicated.values())
        complex_ids = {alt.combination_id for alt in unique if alt.combination_id is not None}
        non_or_member_logic = {
            alt.member_logic_code.strip().upper()
            for alt in unique
            if alt.member_logic_code is not None
            and alt.member_logic_code.strip().upper() != "OR"
        }
        complex_logic = bool(complex_ids or non_or_member_logic)
        if complex_ids:
            group_warnings[(list_name, restraint_id)].append(
                "non-null restraint combination identifier detected; complex AND/OR logic is not flattened"
            )
        if non_or_member_logic:
            codes = ", ".join(sorted(non_or_member_logic))
            group_warnings[(list_name, restraint_id)].append(
                f"NMR-STAR member logic code(s) {codes} detected; only explicit OR is flattened"
            )
        result.append(
            RestraintGroup(
                source_format=detected_format,
                list_name=list_name,
                restraint_id=restraint_id,
                alternatives=unique,
                origin=group_origins[(list_name, restraint_id)],
                complex_logic=complex_logic,
                warnings=group_warnings[(list_name, restraint_id)],
            )
        )
    if not result:
        warnings.append(f"No {target_category} rows matched the selected origin filter.")
    return result


def _parse_nef_row(
    row: dict[str, Any],
    list_name: str,
    row_number: int,
    origin: str | None,
    missing_upper_policy: str,
) -> RawAlternative:
    row_id = pick(row, "index") or str(row_number)
    restraint_id = pick(row, "restraint_id") or f"row-{row_id}"
    endpoint1 = Endpoint(
        chain_code=pick(row, "chain_code_1"),
        sequence_code=pick(row, "sequence_code_1"),
        residue_name=pick(row, "residue_name_1"),
        atom_expression=pick(row, "atom_name_1"),
    )
    endpoint2 = Endpoint(
        chain_code=pick(row, "chain_code_2"),
        sequence_code=pick(row, "sequence_code_2"),
        residue_name=pick(row, "residue_name_2"),
        atom_expression=pick(row, "atom_name_2"),
    )
    upper, bound_source, bound_warnings = _select_upper_bound(
        explicit=as_float(row.get("upper_limit")),
        upper_linear=as_float(row.get("upper_linear_limit")),
        target=as_float(row.get("target_value")),
        uncertainty=as_float(row.get("target_value_uncertainty")),
        policy=missing_upper_policy,
    )
    return RawAlternative(
        source_format="nef",
        list_name=list_name,
        restraint_id=str(restraint_id),
        endpoint1=endpoint1,
        endpoint2=endpoint2,
        upper_bound=upper,
        lower_bound=as_float(row.get("lower_limit")),
        target_value=as_float(row.get("target_value")),
        target_uncertainty=as_float(row.get("target_value_uncertainty")),
        upper_linear_limit=as_float(row.get("upper_linear_limit")),
        weight=as_float(row.get("weight")),
        origin=origin,
        combination_id=pick(row, "restraint_combination_id"),
        row_ids=[str(row_id)],
        bound_source=bound_source,
        warnings=bound_warnings,
    )


def _parse_nmrstar_row(
    row: dict[str, Any],
    list_name: str,
    row_number: int,
    origin: str | None,
    missing_upper_policy: str,
) -> RawAlternative:
    row_id = pick(row, "index_id", "index") or str(row_number)
    restraint_id = pick(row, "id", "constraint_id") or f"row-{row_id}"

    def endpoint(side: int) -> Endpoint:
        suffix = str(side)
        author_atom_expression = pick(
            row,
            f"auth_atom_name_{suffix}",
            f"auth_atom_id_{suffix}",
        )
        canonical_atom = pick(row, f"atom_id_{suffix}")
        return Endpoint(
            chain_code=pick(row, f"auth_asym_id_{suffix}", f"pdb_strand_id_{suffix}"),
            sequence_code=pick(row, f"auth_seq_id_{suffix}"),
            residue_name=pick(row, f"auth_comp_id_{suffix}"),
            atom_expression=author_atom_expression or canonical_atom,
            canonical_chain_code=pick(row, f"entity_assembly_id_{suffix}"),
            canonical_sequence_code=pick(
                row,
                f"comp_index_id_{suffix}",
                f"seq_id_{suffix}",
            ),
            canonical_residue_name=pick(row, f"comp_id_{suffix}"),
            canonical_atom_hint=canonical_atom,
        )

    upper, bound_source, bound_warnings = _select_upper_bound(
        explicit=as_float(row.get("distance_upper_bound_val")),
        upper_linear=as_float(row.get("upper_linear_limit")),
        target=as_float(row.get("target_val")),
        uncertainty=as_float(row.get("target_val_uncertainty")),
        policy=missing_upper_policy,
    )
    return RawAlternative(
        source_format="nmr-star",
        list_name=list_name,
        restraint_id=str(restraint_id),
        endpoint1=endpoint(1),
        endpoint2=endpoint(2),
        upper_bound=upper,
        lower_bound=as_float(row.get("distance_lower_bound_val")),
        target_value=as_float(row.get("target_val")),
        target_uncertainty=as_float(row.get("target_val_uncertainty")),
        upper_linear_limit=as_float(row.get("upper_linear_limit")),
        weight=as_float(row.get("weight")),
        origin=origin,
        combination_id=pick(row, "combination_id"),
        member_id=pick(row, "member_id"),
        member_logic_code=pick(row, "member_logic_code"),
        row_ids=[str(row_id)],
        bound_source=bound_source,
        warnings=bound_warnings,
    )


def _select_upper_bound(
    *,
    explicit: float | None,
    upper_linear: float | None,
    target: float | None,
    uncertainty: float | None,
    policy: str,
) -> tuple[float | None, str, list[str]]:
    if explicit is not None:
        return explicit, "explicit_upper_bound", []
    if policy == "reject":
        return None, "missing", ["no explicit upper bound was supplied"]
    if policy == "upper-linear" and upper_linear is not None:
        return upper_linear, "upper_linear_limit", [
            "upper linear limit used as a heuristic upper bound by explicit policy"
        ]
    if policy == "target-plus-uncertainty" and target is not None:
        if uncertainty is not None and uncertainty < 0:
            return None, "missing", [
                "target uncertainty must be non-negative before it can define an upper bound"
            ]
        derived = target + (uncertainty or 0.0)
        return derived, "target_plus_uncertainty", [
            "target value plus uncertainty used as a heuristic upper bound by explicit policy"
        ]
    if policy == "target" and target is not None:
        return target, "target_value", [
            "target value used as a heuristic upper bound by explicit policy"
        ]
    return None, "missing", [
        f"missing upper-bound policy {policy!r} could not derive a value from available fields"
    ]


def extract_sequence_map(
    entry: Any,
    detected_format: str,
    *,
    chain_map: dict[str, str],
) -> SequenceResolver:
    if detected_format == "nef":
        return _extract_nef_sequence(entry, chain_map)
    return _extract_nmrstar_sequence(entry, chain_map)


def _extract_nef_sequence(entry: Any, chain_map: dict[str, str]) -> SequenceResolver:
    resolver = SequenceResolver()
    rows: list[dict[str, Any]] = []
    for saveframe in entry:
        for loop in saveframe:
            if _category(loop) == "nef_sequence":
                rows.extend(loop_rows(loop))
    if not rows:
        return resolver
    indexed = list(enumerate(rows))
    indexed.sort(
        key=lambda item: (
            int(pick(item[1], "index")) if str(pick(item[1], "index") or "").isdigit() else item[0]
        )
    )
    per_chain: dict[str, int] = defaultdict(int)
    for _, row in indexed:
        if (pick(row, "linking") or "").lower() == "dummy":
            resolver.warnings.append(
                f"NEF dummy residue {pick(row, 'chain_code')}:{pick(row, 'sequence_code')} was not mapped to Boltz"
            )
            continue
        chain = pick(row, "chain_code") or "A"
        sequence = pick(row, "sequence_code")
        residue = pick(row, "residue_name")
        if sequence is None or residue is None:
            resolver.warnings.append("incomplete _nef_sequence row skipped")
            continue
        per_chain[chain] += 1
        record = SequenceRecord(
            source_chain=chain,
            source_sequence_code=sequence,
            residue_name=residue,
            boltz_chain=chain_map.get(chain, chain),
            boltz_residue_index=per_chain[chain],
            source="_nef_sequence",
        )
        resolver.add(record)
    return resolver


def _extract_nmrstar_sequence(entry: Any, chain_map: dict[str, str]) -> SequenceResolver:
    resolver = SequenceResolver()
    assembly_by_id: dict[str, dict[str, str | None]] = {}
    assembly_ids_by_entity: dict[str, list[str]] = defaultdict(list)
    chem_rows: list[dict[str, Any]] = []
    poly_rows: list[dict[str, Any]] = []
    for saveframe in entry:
        for loop in saveframe:
            category = _category(loop)
            rows = loop_rows(loop)
            if category == "entity_assembly":
                for row in rows:
                    assembly_id = pick(row, "id")
                    if assembly_id is None:
                        continue
                    entity_id = pick(row, "entity_id")
                    assembly_by_id[assembly_id] = {
                        "asym": pick(row, "asym_id", "pdb_chain_id"),
                        "pdb": pick(row, "pdb_chain_id"),
                        "entity": entity_id,
                    }
                    if entity_id:
                        assembly_ids_by_entity[entity_id].append(assembly_id)
            elif category == "chem_comp_assembly":
                chem_rows.extend(rows)
            elif category == "entity_poly_seq":
                poly_rows.extend(rows)

    staged: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for order, row in enumerate(chem_rows):
        assembly_id = pick(row, "entity_assembly_id")
        assembly = assembly_by_id.get(assembly_id or "", {})
        author_chain = pick(row, "auth_asym_id")
        chain = author_chain or clean(assembly.get("asym")) or assembly_id or "A"
        comp_index = pick(row, "comp_index_id", "seq_id")
        author_seq = pick(row, "auth_seq_id") or comp_index
        comp_id = pick(row, "comp_id")
        author_comp = pick(row, "auth_comp_id") or comp_id
        if comp_index is None or author_seq is None or comp_id is None or author_comp is None:
            continue
        staged[chain].append(
            {
                "order": order,
                "assembly_id": assembly_id,
                "asym": clean(assembly.get("asym")),
                "author_chain": author_chain,
                "author_seq": author_seq,
                "comp_index": comp_index,
                "comp_id": comp_id,
                "author_comp": author_comp,
            }
        )

    if not staged and poly_rows:
        for order, row in enumerate(poly_rows):
            entity_id = pick(row, "entity_id")
            assemblies = assembly_ids_by_entity.get(entity_id or "", []) or [entity_id or "A"]
            for assembly_id in assemblies:
                assembly = assembly_by_id.get(assembly_id, {})
                chain = clean(assembly.get("asym")) or assembly_id
                num = pick(row, "num", "comp_index_id")
                comp_id = pick(row, "mon_id", "comp_id")
                if num is None or comp_id is None:
                    continue
                staged[chain].append(
                    {
                        "order": order,
                        "assembly_id": assembly_id,
                        "asym": clean(assembly.get("asym")),
                        "author_chain": chain,
                        "author_seq": num,
                        "comp_index": num,
                        "comp_id": comp_id,
                        "author_comp": comp_id,
                    }
                )

    for chain, chain_rows in staged.items():
        chain_rows.sort(key=lambda row: row["order"])
        integer_indices: list[int] = []
        all_integer = True
        for row in chain_rows:
            try:
                integer_indices.append(int(str(row["comp_index"])))
            except ValueError:
                all_integer = False
                break
        contiguous = all_integer and sorted(integer_indices) == list(range(1, len(integer_indices) + 1))
        if not contiguous:
            resolver.warnings.append(
                f"NMR-STAR component indices for chain {chain} were not contiguous 1..N; Boltz indices use sequence-loop order"
            )
        for position, row in enumerate(chain_rows, start=1):
            boltz_index = int(str(row["comp_index"])) if contiguous else position
            source_chain = str(row["author_chain"] or chain)
            source_seq = str(row["author_seq"])
            residue = str(row["author_comp"])
            aliases: list[tuple[str, str, str]] = []
            for alias_chain in sorted(
                {
                    source_chain,
                    str(row["assembly_id"] or ""),
                    str(row["asym"] or ""),
                }
            ):
                if alias_chain:
                    aliases.append((alias_chain, str(row["comp_index"]), str(row["comp_id"])))
                    aliases.append((alias_chain, source_seq, residue))
            resolver.add(
                SequenceRecord(
                    source_chain=source_chain,
                    source_sequence_code=source_seq,
                    residue_name=residue,
                    boltz_chain=chain_map.get(source_chain, chain_map.get(chain, chain)),
                    boltz_residue_index=boltz_index,
                    source="_Chem_comp_assembly" if chem_rows else "_Entity_poly_seq",
                    aliases=aliases,
                )
            )
    return resolver


def infer_sequence_map(
    groups: Iterable[RestraintGroup],
    *,
    chain_map: dict[str, str],
) -> SequenceResolver:
    resolver = SequenceResolver()
    residues: dict[str, dict[tuple[str, str], None]] = defaultdict(dict)
    for group in groups:
        for alternative in group.alternatives:
            for endpoint in (alternative.endpoint1, alternative.endpoint2):
                chain = endpoint.chain_code or endpoint.canonical_chain_code
                sequence = endpoint.sequence_code or endpoint.canonical_sequence_code
                residue = endpoint.residue_name or endpoint.canonical_residue_name
                if chain and sequence and residue:
                    residues[chain][(sequence, residue)] = None
    for chain, values in residues.items():
        ordered = sorted(values, key=lambda item: _natural_sequence_key(item[0]))
        for index, (sequence, residue) in enumerate(ordered, start=1):
            resolver.add(
                SequenceRecord(
                    source_chain=chain,
                    source_sequence_code=sequence,
                    residue_name=residue,
                    boltz_chain=chain_map.get(chain, chain),
                    boltz_residue_index=index,
                    source="inferred-from-restraint-identifiers",
                    warnings=["sequence position inferred; verify against the exact Boltz input sequence"],
                )
            )
    if resolver.records:
        resolver.warnings.append(
            "Sequence mapping was inferred from restraint identifiers, not a molecular-system loop. Verify every Boltz index."
        )
    return resolver


def overlay_residue_map(
    resolver: SequenceResolver,
    path: str | Path,
    *,
    chain_map: dict[str, str],
) -> None:
    map_path = Path(path)
    if map_path.suffix.lower() == ".json":
        payload = json.loads(map_path.read_text(encoding="utf-8"))
        rows = payload["residues"] if isinstance(payload, dict) and "residues" in payload else payload
    else:
        text = map_path.read_text(encoding="utf-8")
        try:
            dialect = csv.Sniffer().sniff(text[:4096], delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel_tab if "\t" in text.partition("\n")[0] else csv.excel
        rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    for row in rows:
        normalized = {str(key).strip().lower(): value for key, value in row.items()}
        source_chain = clean(
            normalized.get("source_chain")
            or normalized.get("chain_code")
            or normalized.get("chain")
        )
        source_seq = clean(
            normalized.get("source_sequence_code")
            or normalized.get("sequence_code")
            or normalized.get("seq")
        )
        residue = clean(
            normalized.get("source_residue_name")
            or normalized.get("residue_name")
            or normalized.get("comp_id")
        )
        boltz_chain = clean(normalized.get("boltz_chain"))
        boltz_index = clean(
            normalized.get("boltz_residue_index")
            or normalized.get("boltz_index")
        )
        if not all((source_chain, source_seq, residue, boltz_index)):
            raise StarDataError(
                "Residue-map rows require source_chain, source_sequence_code, "
                "source_residue_name, and boltz_residue_index."
            )
        try:
            index = int(str(boltz_index))
        except ValueError as exc:
            raise StarDataError(f"Invalid Boltz residue index {boltz_index!r} in {map_path}") from exc
        if index < 1:
            raise StarDataError(
                f"Invalid Boltz residue index {boltz_index!r} in {map_path}; index must be positive."
            )
        resolver.add(
            SequenceRecord(
                source_chain=source_chain,
                source_sequence_code=source_seq,
                residue_name=residue,
                boltz_chain=boltz_chain or chain_map.get(source_chain, source_chain),
                boltz_residue_index=index,
                source=f"user-map:{map_path}",
            ),
            override=True,
        )


def extract_embedded_topologies(entry: Any) -> list[ComponentTopology]:
    atom_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    bond_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for saveframe in entry:
        sf_values = saveframe_tags(saveframe)
        sf_comp = clean(sf_values.get("id")) if str(getattr(saveframe, "category", "")).lstrip("_").lower() == "chem_comp" else None
        for loop in saveframe:
            category = _category(loop)
            if category == "chem_comp_atom":
                for row in loop_rows(loop):
                    comp = pick(row, "comp_id") or sf_comp
                    if comp:
                        atom_rows[comp.upper()].append(row)
            elif category == "chem_comp_bond":
                for row in loop_rows(loop):
                    comp = pick(row, "comp_id") or sf_comp
                    if comp:
                        bond_rows[comp.upper()].append(row)
    topologies: list[ComponentTopology] = []
    for comp in sorted(set(atom_rows) | set(bond_rows)):
        topology = ComponentTopology(comp_id=comp, source="embedded-NMR-STAR-chem-comp")
        for row in atom_rows.get(comp, []):
            topology.add_atom(
                pick(row, "atom_id", "pdb_atom_id"),
                pick(row, "type_symbol", "atom_type"),
            )
        for row in bond_rows.get(comp, []):
            topology.add_bond(
                pick(row, "atom_id_1"),
                pick(row, "atom_id_2"),
                as_float(row.get("value_dist") or row.get("bond_length")),
                as_float(row.get("value_dist_esd") or row.get("bond_length_error")),
            )
        if topology.atom_elements:
            topologies.append(topology)
    return topologies
