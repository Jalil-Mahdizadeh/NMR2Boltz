from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from itertools import chain as itertools_chain
from pathlib import Path
from typing import Any

import yaml

from .model import BoltzAtom, ConversionReport


class TargetValidationError(ValueError):
    """Raised when a Boltz target is malformed or incompatible with a conversion."""


PROTEIN_CODES = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
    "ASX": "B",
    "GLX": "Z",
    "UNK": "X",
}
RNA_CODES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "RA": "A",
    "RC": "C",
    "RG": "G",
    "RU": "U",
    "ADE": "A",
    "CYT": "C",
    "GUA": "G",
    "URA": "U",
}
DNA_CODES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "DA": "A",
    "DC": "C",
    "DG": "G",
    "DT": "T",
    "ADE": "A",
    "CYT": "C",
    "GUA": "G",
    "THY": "T",
}


@dataclass(frozen=True)
class TargetEntity:
    chain: str
    entity_type: str
    sequence: str | None
    modifications: dict[int, str] = field(default_factory=dict)
    ligand_ccd: str | None = None
    ligand_smiles: str | None = None

    @property
    def length(self) -> int:
        return len(self.sequence) if self.sequence is not None else 1


@dataclass(frozen=True)
class ValidationIssue:
    severity: str
    code: str
    message: str
    chain: str | None = None
    residue_index: int | None = None


@dataclass
class TargetValidationResult:
    target_file: str
    target_version: Any
    target_chains: list[str]
    checked_sequence_records: int
    checked_constraints: int
    checked_exact_constraints: int
    checked_union_groups: int
    checked_union_alternatives: int
    mapped_positions: int
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def errors(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.severity == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.severity == "warning"]

    def to_dict(self) -> dict[str, Any]:
        return {
            "target_file": self.target_file,
            "target_version": self.target_version,
            "target_chains": self.target_chains,
            "checked_sequence_records": self.checked_sequence_records,
            "checked_constraints": self.checked_constraints,
            "checked_exact_constraints": self.checked_exact_constraints,
            "checked_union_groups": self.checked_union_groups,
            "checked_union_alternatives": self.checked_union_alternatives,
            "mapped_positions": self.mapped_positions,
            "error_count": len(self.errors),
            "warning_count": len(self.warnings),
            "issues": [asdict(issue) for issue in self.issues],
        }

    def require_valid(self) -> None:
        if not self.errors:
            return
        details = "\n".join(f"- [{issue.code}] {issue.message}" for issue in self.errors[:20])
        remainder = len(self.errors) - 20
        if remainder > 0:
            details += f"\n- ... and {remainder} additional error(s)"
        raise TargetValidationError(f"Boltz target validation failed:\n{details}")


def _nonempty_text(value: Any, description: str) -> str:
    text = str(value).strip() if value is not None else ""
    if not text:
        raise TargetValidationError(f"{description} must be a non-empty value.")
    return text


def _parse_modifications(payload: dict[str, Any], length: int, chain: str) -> dict[int, str]:
    raw = payload.get("modifications", [])
    if raw is None:
        return {}
    if not isinstance(raw, list):
        raise TargetValidationError(f"Target chain {chain} modifications must be a list.")
    result: dict[int, str] = {}
    for item in raw:
        if not isinstance(item, dict) or "position" not in item or "ccd" not in item:
            raise TargetValidationError(
                f"Target chain {chain} modification rows require position and ccd."
            )
        try:
            position = int(item["position"])
        except (TypeError, ValueError) as exc:
            raise TargetValidationError(
                f"Target chain {chain} has invalid modification position {item['position']!r}."
            ) from exc
        if position < 1 or position > length:
            raise TargetValidationError(
                f"Target chain {chain} modification position {position} is outside 1..{length}."
            )
        if position in result:
            raise TargetValidationError(
                f"Target chain {chain} repeats modification position {position}."
            )
        result[position] = _nonempty_text(item["ccd"], f"Modification CCD at {chain}:{position}").upper()
    return result


def load_boltz_target(path: str | Path) -> tuple[Any, dict[str, TargetEntity]]:
    target_path = Path(path)
    try:
        payload = yaml.safe_load(target_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise TargetValidationError(f"Unable to read Boltz target {target_path}: {exc}") from exc
    except yaml.YAMLError as exc:
        raise TargetValidationError(f"Unable to parse Boltz target {target_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise TargetValidationError("Boltz target must be a YAML mapping.")
    version = payload.get("version")
    if version not in {None, 1, "1"}:
        raise TargetValidationError(f"Unsupported Boltz target version {version!r}; expected version 1.")
    sequences = payload.get("sequences")
    if not isinstance(sequences, list) or not sequences:
        raise TargetValidationError("Boltz target requires a non-empty sequences list.")

    entities: dict[str, TargetEntity] = {}
    allowed_types = {"protein", "dna", "rna", "ligand"}
    for position, entry in enumerate(sequences, start=1):
        if not isinstance(entry, dict):
            raise TargetValidationError(f"Target sequence entry {position} must be a mapping.")
        kinds = [key for key in entry if key in allowed_types]
        if len(kinds) != 1:
            raise TargetValidationError(
                f"Target sequence entry {position} must contain exactly one of {sorted(allowed_types)}."
            )
        entity_type = kinds[0]
        entity_payload = entry[entity_type]
        if not isinstance(entity_payload, dict):
            raise TargetValidationError(f"Target {entity_type} entry {position} must be a mapping.")
        raw_ids = entity_payload.get("id")
        ids = raw_ids if isinstance(raw_ids, list) else [raw_ids]
        if not ids:
            raise TargetValidationError(f"Target {entity_type} entry {position} has no chain IDs.")

        sequence: str | None = None
        ligand_ccd: str | None = None
        ligand_smiles: str | None = None
        if entity_type == "ligand":
            has_ccd = entity_payload.get("ccd") is not None
            has_smiles = entity_payload.get("smiles") is not None
            if has_ccd == has_smiles:
                raise TargetValidationError(
                    f"Target ligand entry {position} requires exactly one of ccd or smiles."
                )
            ligand_ccd = (
                _nonempty_text(entity_payload.get("ccd"), "Ligand CCD").upper() if has_ccd else None
            )
            ligand_smiles = (
                _nonempty_text(entity_payload.get("smiles"), "Ligand SMILES") if has_smiles else None
            )
            modifications: dict[int, str] = {}
        else:
            sequence = "".join(_nonempty_text(entity_payload.get("sequence"), "Polymer sequence").split()).upper()
            if not re.fullmatch(r"[A-Z]+", sequence):
                raise TargetValidationError(
                    f"Target {entity_type} sequence entry {position} must contain letters only."
                )
            modifications = _parse_modifications(entity_payload, len(sequence), str(raw_ids))

        for raw_id in ids:
            chain = _nonempty_text(raw_id, f"Target {entity_type} chain ID")
            if chain in entities:
                raise TargetValidationError(f"Duplicate target chain ID {chain!r}.")
            entities[chain] = TargetEntity(
                chain=chain,
                entity_type=entity_type,
                sequence=sequence,
                modifications=dict(modifications),
                ligand_ccd=ligand_ccd,
                ligand_smiles=ligand_smiles,
            )
    return version, entities


def _residue_symbol(residue_name: str, entity_type: str) -> str | None:
    name = residue_name.strip().upper()
    if entity_type == "protein":
        if len(name) == 1 and name in set(PROTEIN_CODES.values()):
            return name
        return PROTEIN_CODES.get(name)
    if entity_type == "rna":
        return RNA_CODES.get(name)
    if entity_type == "dna":
        return DNA_CODES.get(name)
    return None


def _validate_position(
    *,
    chain: str,
    residue_index: int,
    residue_name: str | None,
    entities: dict[str, TargetEntity],
    issues: list[ValidationIssue],
    check_identity: bool,
) -> bool:
    entity = entities.get(chain)
    if entity is None:
        issues.append(
            ValidationIssue(
                "error",
                "target_chain_missing",
                f"Mapped chain {chain!r} is absent from the Boltz target.",
                chain,
                residue_index,
            )
        )
        return False
    if residue_index < 1 or residue_index > entity.length:
        issues.append(
            ValidationIssue(
                "error",
                "target_index_out_of_range",
                f"Mapped position {chain}:{residue_index} is outside target length {entity.length}.",
                chain,
                residue_index,
            )
        )
        return False
    if not check_identity or residue_name is None:
        return True

    actual = residue_name.strip().upper()
    modification = entity.modifications.get(residue_index)
    if modification is not None:
        if actual != modification:
            issues.append(
                ValidationIssue(
                    "error",
                    "target_modification_mismatch",
                    f"Source residue {actual} at {chain}:{residue_index} does not match target modification {modification}.",
                    chain,
                    residue_index,
                )
            )
        return True
    if entity.entity_type == "ligand":
        if entity.ligand_ccd is not None and actual != entity.ligand_ccd:
            issues.append(
                ValidationIssue(
                    "error",
                    "target_ligand_mismatch",
                    f"Source residue {actual} at {chain}:1 does not match target ligand CCD {entity.ligand_ccd}.",
                    chain,
                    residue_index,
                )
            )
        elif entity.ligand_smiles is not None:
            issues.append(
                ValidationIssue(
                    "warning",
                    "smiles_ligand_identity_unverified",
                    f"Cannot verify source residue {actual} against a SMILES-only target ligand at {chain}:1.",
                    chain,
                    residue_index,
                )
            )
        return True

    expected = entity.sequence[residue_index - 1] if entity.sequence is not None else None
    observed = _residue_symbol(actual, entity.entity_type)
    if expected == "X":
        issues.append(
            ValidationIssue(
                "warning",
                "ambiguous_target_residue",
                f"Target residue {chain}:{residue_index} is X; identity of source residue {actual} was not proven.",
                chain,
                residue_index,
            )
        )
    elif observed is None:
        issues.append(
            ValidationIssue(
                "error",
                "unsupported_source_residue",
                f"Source residue {actual} at {chain}:{residue_index} is not canonical for target type {entity.entity_type} and is not declared as a modification.",
                chain,
                residue_index,
            )
        )
    elif observed != expected:
        issues.append(
            ValidationIssue(
                "error",
                "target_residue_mismatch",
                f"Source residue {actual} ({observed}) at {chain}:{residue_index} does not match target residue {expected}.",
                chain,
                residue_index,
            )
        )
    return True


def validate_report_against_target(
    report: ConversionReport,
    target_path: str | Path,
) -> TargetValidationResult:
    version, entities = load_boltz_target(target_path)
    issues: list[ValidationIssue] = []
    mapped: dict[tuple[str, int], str] = {}
    inferred_mapping = False
    for record in report.sequence_map:
        key = (record.boltz_chain, record.boltz_residue_index)
        residue = record.residue_name.strip().upper()
        prior = mapped.get(key)
        if prior is not None and prior != residue:
            issues.append(
                ValidationIssue(
                    "error",
                    "target_mapping_collision",
                    f"Multiple source residues map to {key[0]}:{key[1]}: {prior} and {residue}.",
                    key[0],
                    key[1],
                )
            )
        mapped[key] = residue
        inferred_mapping = inferred_mapping or record.source == "inferred-from-restraint-identifiers"
        _validate_position(
            chain=key[0],
            residue_index=key[1],
            residue_name=residue,
            entities=entities,
            issues=issues,
            check_identity=True,
        )

    if inferred_mapping:
        issues.append(
            ValidationIssue(
                "warning",
                "inferred_sequence_mapping",
                "At least one sequence position was inferred from restraint identifiers; target matching reduces but does not eliminate numbering risk.",
            )
        )

    checked_atoms: set[BoltzAtom] = set()
    constraint_atoms = (
        atom
        for constraint in report.emitted_constraints
        for atom in (constraint.atom1, constraint.atom2)
    )
    union_atoms = (
        atom
        for group in report.ambiguous_groups
        for alternative in group.alternatives
        for atom in (alternative.atom1, alternative.atom2)
    )
    for atom in itertools_chain(constraint_atoms, union_atoms):
        if atom in checked_atoms:
            continue
        checked_atoms.add(atom)
        _validate_position(
            chain=atom.chain,
            residue_index=atom.residue_index,
            residue_name=None,
            entities=entities,
            issues=issues,
            check_identity=False,
        )

    exact_count = len(report.emitted_constraints)
    union_group_count = len(report.ambiguous_groups)
    union_alternative_count = sum(
        len(group.alternatives) for group in report.ambiguous_groups
    )

    for chain, entity in entities.items():
        mapped_count = sum(mapped_chain == chain for mapped_chain, _index in mapped)
        if 0 < mapped_count < entity.length:
            issues.append(
                ValidationIssue(
                    "warning",
                    "partial_target_coverage",
                    f"Sequence map covers {mapped_count} of {entity.length} positions in target chain {chain}.",
                    chain,
                )
            )

    return TargetValidationResult(
        target_file=str(Path(target_path)),
        target_version=version,
        target_chains=sorted(entities),
        checked_sequence_records=len(report.sequence_map),
        checked_constraints=exact_count + union_group_count,
        checked_exact_constraints=exact_count,
        checked_union_groups=union_group_count,
        checked_union_alternatives=union_alternative_count,
        mapped_positions=len(mapped),
        issues=issues,
    )
