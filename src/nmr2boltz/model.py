from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

NULL_VALUES = {None, "", ".", "?"}


def clean(value: Any) -> str | None:
    """Return a stripped STAR value or ``None`` for STAR null tokens."""
    if value is None:
        return None
    text = str(value).strip()
    return None if text in {"", ".", "?"} else text


def as_float(value: Any) -> float | None:
    text = clean(value)
    if text is None:
        return None
    try:
        result = float(text)
    except (TypeError, ValueError):
        return None
    if result != result or result in {float("inf"), float("-inf")}:
        return None
    return result


@dataclass(frozen=True, order=True)
class BoltzAtom:
    chain: str
    residue_index: int
    atom_name: str

    def display(self) -> str:
        return f"{self.chain}:{self.residue_index}:{self.atom_name}"


@dataclass(frozen=True, order=True)
class BoltzToken:
    chain: str
    residue_index: int

    def display(self) -> str:
        return f"{self.chain}:{self.residue_index}"


@dataclass
class SequenceRecord:
    source_chain: str
    source_sequence_code: str
    residue_name: str
    boltz_chain: str
    boltz_residue_index: int
    source: str
    aliases: list[tuple[str, str, str]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


@dataclass
class Endpoint:
    chain_code: str | None
    sequence_code: str | None
    residue_name: str | None
    atom_expression: str | None
    canonical_chain_code: str | None = None
    canonical_sequence_code: str | None = None
    canonical_residue_name: str | None = None
    canonical_atom_hint: str | None = None
    canonical_atom_set: list[str] = field(default_factory=list)

    def display(self) -> str:
        return ":".join(
            [
                self.chain_code or "?",
                self.sequence_code or "?",
                self.residue_name or "?",
                self.atom_expression or "?",
            ]
        )


@dataclass
class RawAlternative:
    source_format: str
    list_name: str
    restraint_id: str
    endpoint1: Endpoint
    endpoint2: Endpoint
    upper_bound: float | None
    lower_bound: float | None = None
    target_value: float | None = None
    target_uncertainty: float | None = None
    upper_linear_limit: float | None = None
    weight: float | None = None
    origin: str | None = None
    combination_id: str | None = None
    member_id: str | None = None
    member_logic_code: str | None = None
    row_ids: list[str] = field(default_factory=list)
    canonical_expansions: list[dict[str, str | None]] = field(default_factory=list)
    bound_source: str = "explicit_upper_bound"
    warnings: list[str] = field(default_factory=list)

    def dedup_key(self) -> tuple[Any, ...]:
        def endpoint_key(endpoint: Endpoint) -> tuple[Any, ...]:
            expression = endpoint.atom_expression or ""
            # NMR-STAR conversions of XPLOR/CNS restraints commonly retain
            # "#" as the author-level digit wildcard. Topology resolution
            # normalizes it to NEF "%"; recognize it here as well so repeated
            # canonical expansion rows are collapsed before projection.
            is_set_expression = any(
                symbol in expression for symbol in ("%", "*", "#", "x", "y")
            )
            # Canonical expansion rows for an author atom set are deduplicated.
            # For a plain author atom name, differing canonical atoms remain
            # separate alternatives rather than being silently discarded.
            canonical_hint = None if is_set_expression else endpoint.canonical_atom_hint
            return (
                endpoint.chain_code,
                endpoint.sequence_code,
                endpoint.residue_name,
                endpoint.atom_expression,
                canonical_hint,
                tuple(endpoint.canonical_atom_set),
            )

        e1 = endpoint_key(self.endpoint1)
        e2 = endpoint_key(self.endpoint2)
        ends = tuple(sorted((e1, e2), key=lambda item: tuple(x or "" for x in item)))
        return (
            ends,
            self.upper_bound,
            self.lower_bound,
            self.target_value,
            self.target_uncertainty,
            self.upper_linear_limit,
            self.weight,
            self.combination_id,
            self.bound_source,
        )


@dataclass
class RestraintGroup:
    source_format: str
    list_name: str
    restraint_id: str
    alternatives: list[RawAlternative]
    origin: str | None = None
    complex_logic: bool = False
    warnings: list[str] = field(default_factory=list)

    @property
    def group_id(self) -> str:
        return f"{self.list_name}:{self.restraint_id}"


@dataclass
class AtomChoice:
    atom_name: str
    element: str
    parent_atom: str
    bond_length_upper: float
    resolution_source: str
    warnings: list[str] = field(default_factory=list)


@dataclass
class AtomSetChoice:
    atoms: list[AtomChoice]
    expression: str
    assignment_key: str | None = None
    semantics: str = "explicit"
    warnings: list[str] = field(default_factory=list)


@dataclass
class ProjectedAlternative:
    atom1: BoltzAtom
    atom2: BoltzAtom
    max_distance: float
    source_upper_bound: float
    averaging_policy: str
    averaging_factor: float
    explicit_pair_count: int
    bond_offset: float
    group_id: str
    source_observation: dict[str, Any] = field(default_factory=dict)
    source_rows: list[str] = field(default_factory=list)
    source_endpoints: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    raw_projected_distance: float | None = None
    boltz_adjustment: str | None = None

    @property
    def pair_key(self) -> tuple[BoltzAtom, BoltzAtom]:
        return tuple(sorted((self.atom1, self.atom2)))  # type: ignore[return-value]


@dataclass
class EmittedConstraint:
    atom1: BoltzAtom
    atom2: BoltzAtom
    max_distance: float
    source_groups: list[str]
    raw_projected_distance: float
    boltz_adjustment: str | None = None
    provenance: list[dict[str, Any]] = field(default_factory=list)

    @property
    def pair_key(self) -> tuple[BoltzAtom, BoltzAtom]:
        return tuple(sorted((self.atom1, self.atom2)))  # type: ignore[return-value]


@dataclass
class AmbiguousGroup:
    group_id: str
    restraint_id: str
    list_name: str
    alternatives: list[ProjectedAlternative]
    reason: str
    source_format: str
    origin: str | None = None
    warnings: list[str] = field(default_factory=list)


@dataclass
class Rejection:
    group_id: str
    reason: str
    details: str
    row_ids: list[str] = field(default_factory=list)
    endpoint: str | None = None
    provenance: dict[str, Any] = field(default_factory=dict)


@dataclass
class TokenContribution:
    source_kind: str
    source_groups: list[str]
    raw_bound: float
    candidate_bound: float
    adjustments: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)


@dataclass
class TokenConstraint:
    token1: BoltzToken
    token2: BoltzToken
    max_distance: float
    raw_candidate_bounds: list[float]
    source_groups: list[str]
    source_kinds: list[str]
    adjustments: list[str] = field(default_factory=list)
    contributions: list[TokenContribution] = field(default_factory=list)

    @property
    def pair_key(self) -> tuple[BoltzToken, BoltzToken]:
        return tuple(sorted((self.token1, self.token2)))  # type: ignore[return-value]


@dataclass
class TokenProjectionOmission:
    source_kind: str
    source_groups: list[str]
    reason: str
    details: str
    raw_bound: float | None = None
    provenance: dict[str, Any] = field(default_factory=dict)


@dataclass
class ConversionReport:
    input_file: str
    detected_format: str
    settings: dict[str, Any]
    sequence_map: list[SequenceRecord]
    emitted_constraints: list[EmittedConstraint]
    ambiguous_groups: list[AmbiguousGroup]
    rejections: list[Rejection]
    warnings: list[str]
    statistics: dict[str, Any]
    source_restraint_groups: list[RestraintGroup] = field(default_factory=list)
    target_component_topologies: dict[str, dict[str, Any]] = field(default_factory=dict)
    target_validation: dict[str, Any] | None = None
    token_constraints: list[TokenConstraint] = field(default_factory=list)
    token_projection_omissions: list[TokenProjectionOmission] = field(default_factory=list)
    token_projection_statistics: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        # Canonical expansion provenance is optional and is populated only
        # when multiple canonical rows collapse to one author-level set.
        for group in payload["source_restraint_groups"]:
            for alternative in group["alternatives"]:
                if not alternative["canonical_expansions"]:
                    del alternative["canonical_expansions"]
        # Union-bound adjustment provenance is present only when the producer
        # had to weaken a projected alternative to the Boltz minimum.
        for group in payload["ambiguous_groups"]:
            for alternative in group["alternatives"]:
                if alternative["raw_projected_distance"] is None:
                    del alternative["raw_projected_distance"]
                if alternative["boltz_adjustment"] is None:
                    del alternative["boltz_adjustment"]
        return payload
