from __future__ import annotations

import itertools
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from .model import AtomChoice, AtomSetChoice, BoltzAtom, ConversionReport, SequenceRecord, clean

# Conservative upper envelopes for ordinary covalent X-H bond lengths (angstrom).
# They are deliberately a little longer than common idealized geometry values.
DEFAULT_XH_UPPER = {
    "C": 1.12,
    "N": 1.08,
    "O": 1.02,
    "S": 1.36,
    "SE": 1.49,
    "P": 1.45,
}
FALLBACK_XH_UPPER = 1.55
HYDROGEN_ELEMENTS = {"H", "D", "T"}


class TopologyResolutionError(ValueError):
    pass


class AtomTopologyValidationError(ValueError):
    """Raised when an executable constraint lacks exact component-topology evidence."""


@dataclass(frozen=True)
class AtomTopologyViolation:
    atom: BoltzAtom
    residue_name: str | None
    topology_source: str | None
    reason: str

    def to_dict(self) -> dict[str, str | int | None]:
        return {
            "chain": self.atom.chain,
            "residue_number": self.atom.residue_index,
            "residue_name": self.residue_name,
            "atom_name": self.atom.atom_name,
            "topology_source": self.topology_source,
            "reason": self.reason,
        }


@dataclass
class ComponentTopology:
    comp_id: str
    atom_elements: dict[str, str] = field(default_factory=dict)
    hydrogen_parent: dict[str, str] = field(default_factory=dict)
    hydrogen_bond_upper: dict[str, float] = field(default_factory=dict)
    atom_aliases: dict[str, str] = field(default_factory=dict)
    source: str = "unknown"

    def add_atom(self, atom_name: str, element: str | None) -> None:
        atom = clean(atom_name)
        if atom is None:
            return
        elem = (clean(element) or infer_element(atom)).upper()
        self.atom_elements[atom] = elem

    def add_bond(
        self,
        atom1: str,
        atom2: str,
        distance: float | None = None,
        distance_error: float | None = None,
    ) -> None:
        a1 = clean(atom1)
        a2 = clean(atom2)
        if a1 is None or a2 is None:
            return
        e1 = self.atom_elements.get(a1, infer_element(a1)).upper()
        e2 = self.atom_elements.get(a2, infer_element(a2)).upper()
        if e1 in HYDROGEN_ELEMENTS and e2 not in HYDROGEN_ELEMENTS:
            hydrogen, parent, parent_element = a1, a2, e2
        elif e2 in HYDROGEN_ELEMENTS and e1 not in HYDROGEN_ELEMENTS:
            hydrogen, parent, parent_element = a2, a1, e1
        else:
            return
        default = DEFAULT_XH_UPPER.get(parent_element, FALLBACK_XH_UPPER)
        if distance is not None and math.isfinite(distance) and distance > 0:
            error = distance_error if distance_error is not None and math.isfinite(distance_error) else 0.0
            candidate = distance + max(0.03, 3.0 * max(0.0, error))
            upper = max(default, candidate)
        else:
            upper = default
        self.hydrogen_parent[hydrogen] = parent
        self.hydrogen_bond_upper[hydrogen] = upper
        self.atom_elements.setdefault(hydrogen, "H")
        self.atom_elements.setdefault(parent, parent_element)

    def add_hydrogen(self, hydrogen: str, parent: str, parent_element: str | None = None) -> None:
        elem = (parent_element or infer_element(parent)).upper()
        self.atom_elements[hydrogen] = "H"
        self.atom_elements.setdefault(parent, elem)
        self.hydrogen_parent[hydrogen] = parent
        self.hydrogen_bond_upper[hydrogen] = DEFAULT_XH_UPPER.get(elem, FALLBACK_XH_UPPER)

    def add_alias(self, alias: str, canonical_atom: str) -> None:
        self.atom_aliases[alias] = canonical_atom

    def available_atoms(self) -> set[str]:
        # Aliases are intentionally excluded from pattern expansion so wildcard
        # cardinality counts physical atoms, not duplicate naming conventions.
        return set(self.atom_elements) | set(self.hydrogen_parent) | set(self.hydrogen_parent.values())

    def atom_choice(self, atom_name: str, resolution_source: str) -> AtomChoice:
        element = self.atom_elements.get(atom_name, infer_element(atom_name)).upper()
        if element in HYDROGEN_ELEMENTS or atom_name in self.hydrogen_parent:
            parent = self.hydrogen_parent.get(atom_name)
            if parent is None:
                raise TopologyResolutionError(
                    f"Hydrogen {self.comp_id}:{atom_name} has no uniquely resolved heavy-atom parent."
                )
            parent_element = self.atom_elements.get(parent, infer_element(parent)).upper()
            return AtomChoice(
                atom_name=atom_name,
                element="H",
                parent_atom=parent,
                bond_length_upper=self.hydrogen_bond_upper.get(
                    atom_name, DEFAULT_XH_UPPER.get(parent_element, FALLBACK_XH_UPPER)
                ),
                resolution_source=resolution_source,
            )
        return AtomChoice(
            atom_name=atom_name,
            element=element,
            parent_atom=atom_name,
            bond_length_upper=0.0,
            resolution_source=resolution_source,
        )


def infer_element(atom_name: str) -> str:
    name = atom_name.strip().upper()
    if not name:
        return "X"
    if name[0].isdigit():
        name = name.lstrip("0123456789")
    if not name:
        return "X"
    # Atom names in biomolecular dictionaries are normally element-first.  Check
    # common two-letter elements before falling back to the first character.
    two = name[:2]
    if two in {"SE", "CL", "BR", "NA", "MG", "ZN", "FE", "MN", "CO", "NI", "CU"}:
        return two
    if name[0] in {"D", "T"}:
        return "H"
    return name[0]


def _leading_digit_alias(name: str) -> str | None:
    match = re.fullmatch(r"([123])([HDT].+)", name)
    if not match:
        return None
    return f"{match.group(2)}{match.group(1)}"


def _normalize_legacy_expression(expression: str) -> tuple[str, list[str]]:
    warnings: list[str] = []
    result = expression.strip()
    if "#" in result:
        result = result.replace("#", "%")
        warnings.append("legacy '#' atom wildcard interpreted as NEF '%' wildcard")
    if '"' in result:
        result = result.replace('"', "''")
        warnings.append("double-quote prime notation normalized to two apostrophes")
    return result, warnings


def _compile_atom_pattern(expression: str) -> tuple[re.Pattern[str], list[str]]:
    pieces: list[str] = ["^"]
    xy_names: list[str] = []
    xy_count = 0
    for char in expression:
        if char == "%":
            pieces.append(r"[0-9]+")
        elif char == "*":
            pieces.append(r"\S*")
        elif char in {"x", "y"}:
            group_name = f"xy{xy_count}"
            xy_count += 1
            xy_names.append(group_name)
            # Non-greedy is important for HGx%: HG11 should resolve branch x=1,
            # while the '%' consumes the final proton number.
            # NEF x/y replaces a stereochemical numeric suffix. Apostrophes
            # are literal parts of nucleotide atom names (for example H4')
            # and must never be consumed as an assignment branch.
            pieces.append(fr"(?P<{group_name}>[0-9]+?)")
        else:
            pieces.append(re.escape(char))
    pieces.append("$")
    return re.compile("".join(pieces)), xy_names


def _pseudoatom_to_expression(name: str) -> str | None:
    # Optional approximation only. NEF defines Q/M pseudoatoms as geometric
    # pseudoatoms, not atom sets, so the default policy rejects these names.
    match = re.fullmatch(r"[QM]([A-Z][A-Z0-9']*)", name.upper())
    if not match:
        return None
    suffix = match.group(1)
    return f"H{suffix}%"


class TopologyLibrary:
    def __init__(
        self,
        external_ccd_paths: Iterable[str | Path] = (),
        bond_length_config: str | Path | None = None,
    ) -> None:
        self.components: dict[str, ComponentTopology] = build_builtin_topologies()
        self.external_ccd_paths = [Path(path) for path in external_ccd_paths]
        self._loaded_external_files: set[Path] = set()
        self._loaded_external_components: set[tuple[Path, str]] = set()
        self.element_bond_overrides: dict[str, float] = {}
        self.component_bond_overrides: dict[str, dict[str, float]] = {}
        if bond_length_config is not None:
            self._load_bond_length_config(Path(bond_length_config))
        for topology in self.components.values():
            self._apply_bond_overrides(topology)

    def _load_bond_length_config(self, path: Path) -> None:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            raise TopologyResolutionError(f"Unable to read bond-length configuration {path}: {exc}") from exc
        for element, value in payload.get("element_upper", {}).items():
            try:
                numeric = float(value)
            except (TypeError, ValueError) as exc:
                raise TopologyResolutionError(
                    f"Bond-length upper value for element {element} must be numeric."
                ) from exc
            if not math.isfinite(numeric) or numeric <= 0:
                raise TopologyResolutionError(
                    f"Bond-length upper value for element {element} must be positive and finite."
                )
            self.element_bond_overrides[str(element).upper()] = numeric
        for comp_id, mapping in payload.get("components", {}).items():
            comp_values: dict[str, float] = {}
            for hydrogen, value in mapping.items():
                try:
                    numeric = float(value)
                except (TypeError, ValueError) as exc:
                    raise TopologyResolutionError(
                        f"Bond-length upper value for {comp_id}:{hydrogen} must be numeric."
                    ) from exc
                if not math.isfinite(numeric) or numeric <= 0:
                    raise TopologyResolutionError(
                        f"Bond-length upper value for {comp_id}:{hydrogen} must be positive and finite."
                    )
                comp_values[str(hydrogen)] = numeric
            self.component_bond_overrides[str(comp_id).upper()] = comp_values

    def _apply_bond_overrides(self, topology: ComponentTopology) -> None:
        for hydrogen, parent in topology.hydrogen_parent.items():
            parent_element = topology.atom_elements.get(parent, infer_element(parent)).upper()
            if parent_element in self.element_bond_overrides:
                topology.hydrogen_bond_upper[hydrogen] = self.element_bond_overrides[parent_element]
        for hydrogen, value in self.component_bond_overrides.get(topology.comp_id.upper(), {}).items():
            if hydrogen not in topology.hydrogen_parent:
                raise TopologyResolutionError(
                    f"Bond-length configuration names unknown hydrogen {topology.comp_id}:{hydrogen}."
                )
            topology.hydrogen_bond_upper[hydrogen] = value

    def register(self, topology: ComponentTopology, replace: bool = True) -> None:
        key = topology.comp_id.upper()
        self._apply_bond_overrides(topology)
        if replace or key not in self.components:
            self.components[key] = topology

    def get(self, comp_id: str) -> ComponentTopology | None:
        key = comp_id.upper()
        topology = self.components.get(key)
        if topology is not None:
            return topology
        self._load_external_for_component(key)
        return self.components.get(key)

    def _load_external_for_component(self, comp_id: str) -> None:
        for path in self.external_ccd_paths:
            if path.is_dir():
                candidates = [
                    path / f"{comp_id}.cif",
                    path / f"{comp_id.lower()}.cif",
                    path / comp_id[0].lower() / f"{comp_id}.cif",
                ]
                for candidate in candidates:
                    if candidate.exists():
                        self._load_ccd_file(candidate, only_component=comp_id)
                        if comp_id in self.components:
                            return
            elif path.exists() and (path, comp_id) not in self._loaded_external_components:
                self._load_ccd_file(path, only_component=comp_id)
                self._loaded_external_components.add((path, comp_id))
                if comp_id in self.components:
                    return

    def _load_ccd_file(self, path: Path, only_component: str | None = None) -> None:
        try:
            import gemmi
        except ImportError as exc:  # pragma: no cover - dependency is in the distribution
            raise TopologyResolutionError(
                "Reading external CCD mmCIF files requires the 'gemmi' package."
            ) from exc
        try:
            document = gemmi.cif.read_file(str(path))
        except Exception as exc:  # pragma: no cover - gemmi error details vary
            raise TopologyResolutionError(f"Unable to read CCD file {path}: {exc}") from exc
        self._loaded_external_files.add(path)
        for block in document:
            comp_id = clean(block.find_value("_chem_comp.id"))
            if comp_id is None:
                name = str(block.name)
                comp_id = name[5:] if name.lower().startswith("comp_") else name
            key = comp_id.upper()
            if only_component and key != only_component.upper():
                continue
            topology = ComponentTopology(comp_id=key, source=f"ccd:{path}")
            try:
                atom_table = block.find(
                    ["_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol"]
                )
                for row in atom_table:
                    topology.add_atom(str(row[0]), str(row[1]))
            except Exception:
                pass
            bond_tags = [
                "_chem_comp_bond.atom_id_1",
                "_chem_comp_bond.atom_id_2",
                "_chem_comp_bond.value_dist",
                "_chem_comp_bond.value_dist_esd",
            ]
            try:
                bond_table = block.find(bond_tags)
                for row in bond_table:
                    distance = _optional_float(str(row[2]))
                    error = _optional_float(str(row[3]))
                    topology.add_bond(str(row[0]), str(row[1]), distance, error)
            except Exception:
                try:
                    bond_table = block.find(bond_tags[:2])
                    for row in bond_table:
                        topology.add_bond(str(row[0]), str(row[1]))
                except Exception:
                    pass
            if topology.atom_elements:
                self.register(topology)

    def resolve_expression(
        self,
        comp_id: str,
        expression: str,
        canonical_hint: str | None = None,
        pseudoatom_policy: str = "reject",
    ) -> list[AtomSetChoice]:
        comp = comp_id.upper()
        topology = self.get(comp)
        original_expression = expression
        expression, normalization_warnings = _normalize_legacy_expression(expression)

        if topology is None:
            # A named heavy atom can still be passed through without component
            # topology. Unknown hydrogens cannot be projected safely.
            hint = canonical_hint or expression
            if not _looks_like_hydrogen(hint):
                warning = f"component {comp} topology unavailable; heavy atom name passed through"
                choice = AtomChoice(
                    atom_name=hint,
                    element=infer_element(hint),
                    parent_atom=hint,
                    bond_length_upper=0.0,
                    resolution_source="unverified-heavy-name",
                    warnings=[warning],
                )
                return [
                    AtomSetChoice(
                        atoms=[choice],
                        expression=original_expression,
                        semantics="explicit-heavy-unverified",
                        warnings=normalization_warnings + [warning],
                    )
                ]
            raise TopologyResolutionError(
                f"No topology is available for component {comp}; cannot map hydrogen {expression}."
            )

        available = topology.available_atoms()
        pseudo_expression_early = _pseudoatom_to_expression(expression)
        if expression not in available and pseudo_expression_early is not None and pseudoatom_policy == "reject":
            canonical_note = (
                f" Canonical Atom_ID {canonical_hint!r} was preserved for audit but is not "
                "trusted as a complete pseudoatom expansion."
                if canonical_hint
                else ""
            )
            raise TopologyResolutionError(
                f"{comp}:{expression} appears to be a geometric pseudoatom. "
                "Geometric pseudoatoms are not equivalent to wildcard atom sets; use "
                "--pseudoatom-policy atomset only after expert review."
                + canonical_note
            )
        has_wildcard = any(symbol in expression for symbol in ("%", "*"))
        has_xy = "x" in expression or "y" in expression
        aliases = [expression]
        leading_alias = _leading_digit_alias(expression)
        if leading_alias:
            aliases.append(leading_alias)
        if expression == "HN":
            aliases.append("H")
        # A canonical Atom_ID is a fallback for an explicit author name. It must
        # not override an author wildcard/x-y expression, because translated
        # NMR-STAR commonly repeats one wildcard as several canonical rows.
        if canonical_hint and not (has_wildcard or has_xy):
            aliases.append(canonical_hint)

        # Exact resolution precedes wildcard processing. This protects unusual
        # components that legitimately contain characters with wildcard meaning.
        for candidate in aliases:
            resolved_candidate = topology.atom_aliases.get(candidate, candidate)
            if resolved_candidate in available:
                choice = topology.atom_choice(resolved_candidate, f"{topology.source}:exact")
                warnings = list(normalization_warnings)
                if resolved_candidate != expression:
                    warnings.append(
                        f"atom name {expression!r} resolved through naming alias {resolved_candidate!r}"
                    )
                return [
                    AtomSetChoice(
                        atoms=[choice],
                        expression=original_expression,
                        semantics="explicit",
                        warnings=warnings,
                    )
                ]

        # Case-insensitive fallback is allowed only when unique, and is always
        # reported because NEF identifiers are formally case-sensitive.
        lower_lookup: dict[str, set[str]] = {}
        for name in available:
            lower_lookup.setdefault(name.lower(), set()).add(name)
        for alias, canonical in topology.atom_aliases.items():
            lower_lookup.setdefault(alias.lower(), set()).add(canonical)
        lower_matches = sorted(lower_lookup.get(expression.lower(), set()))
        if len(lower_matches) == 1:
            candidate = lower_matches[0]
            choice = topology.atom_choice(candidate, f"{topology.source}:unique-case-fold")
            return [
                AtomSetChoice(
                    atoms=[choice],
                    expression=original_expression,
                    semantics="explicit",
                    warnings=normalization_warnings
                    + [f"case-sensitive atom name {expression!r} resolved uniquely as {candidate!r}"],
                )
            ]

        if has_wildcard or has_xy:
            pattern, xy_names = _compile_atom_pattern(expression)
            grouped: dict[tuple[str, ...], list[str]] = {}
            for atom_name in sorted(available):
                match = pattern.fullmatch(atom_name)
                if not match:
                    continue
                # Expressions beginning with H/D/T denote nuclei, not a heavy
                # element whose atom name merely happens to begin with H.
                if expression[:1].upper() in {"H", "D", "T"}:
                    element = topology.atom_elements.get(atom_name, infer_element(atom_name)).upper()
                    if element not in HYDROGEN_ELEMENTS and atom_name not in topology.hydrogen_parent:
                        continue
                key = tuple(match.group(name) for name in xy_names)
                grouped.setdefault(key, []).append(atom_name)
            if grouped:
                results: list[AtomSetChoice] = []
                for key, names in sorted(grouped.items()):
                    atoms = [
                        topology.atom_choice(name, f"{topology.source}:pattern") for name in names
                    ]
                    if has_xy and has_wildcard:
                        semantics = "atom-set-with-stereo-branch"
                    elif has_xy:
                        semantics = "stereo-assignment-alternative"
                    else:
                        semantics = "atom-set"
                    results.append(
                        AtomSetChoice(
                            atoms=atoms,
                            expression=original_expression,
                            assignment_key="/".join(key) if key else None,
                            semantics=semantics,
                            warnings=list(normalization_warnings),
                        )
                    )
                return results

        if canonical_hint and canonical_hint in available:
            choice = topology.atom_choice(canonical_hint, f"{topology.source}:canonical-fallback")
            return [
                AtomSetChoice(
                    atoms=[choice],
                    expression=original_expression,
                    semantics="canonical-fallback",
                    warnings=normalization_warnings
                    + [
                        f"author atom expression {expression!r} did not match topology; "
                        f"canonical Atom_ID {canonical_hint!r} used as a single fallback"
                    ],
                )
            ]

        pseudo_expression = _pseudoatom_to_expression(expression)
        if pseudo_expression is not None:
            if pseudoatom_policy == "reject":
                raise TopologyResolutionError(
                    f"{comp}:{expression} appears to be a geometric pseudoatom. "
                    "Geometric pseudoatoms are not equivalent to wildcard atom sets; use "
                    "--pseudoatom-policy atomset only after expert review."
                )
            if pseudoatom_policy == "atomset":
                results = self.resolve_expression(
                    comp,
                    pseudo_expression,
                    canonical_hint=None,
                    pseudoatom_policy="reject",
                )
                for result in results:
                    result.expression = original_expression
                    result.semantics = "pseudoatom-approximated-as-atom-set"
                    result.warnings.append(
                        f"geometric pseudoatom {expression!r} approximated as {pseudo_expression!r} atom set"
                    )
                return results

        if not _looks_like_hydrogen(expression):
            # Heavy atoms can be valid even when a compact built-in topology did
            # not enumerate them. Passing one through does not require a bond map.
            choice = AtomChoice(
                atom_name=expression,
                element=topology.atom_elements.get(expression, infer_element(expression)),
                parent_atom=expression,
                bond_length_upper=0.0,
                resolution_source=f"{topology.source}:unlisted-heavy",
                warnings=["heavy atom was not enumerated in the selected topology"],
            )
            return [
                AtomSetChoice(
                    atoms=[choice],
                    expression=original_expression,
                    semantics="explicit-heavy-unlisted",
                    warnings=normalization_warnings + choice.warnings,
                )
            ]

        available_h = sorted(topology.hydrogen_parent)
        preview = ", ".join(available_h[:16]) + (" ..." if len(available_h) > 16 else "")
        raise TopologyResolutionError(
            f"Unable to resolve {comp}:{expression}. Available mapped hydrogens include: {preview or 'none'}."
        )

    def resolve_canonical_atom_set(
        self,
        comp_id: str,
        atom_names: Iterable[str],
        *,
        author_expression: str,
    ) -> AtomSetChoice:
        """Resolve a topology-verified canonical expansion as one physical set."""
        topology = self.get(comp_id.upper())
        if topology is None:
            raise TopologyResolutionError(
                f"No topology is available for component {comp_id.upper()}; "
                "cannot verify a canonical atom-set expansion."
            )
        names = list(dict.fromkeys(str(name) for name in atom_names))
        if len(names) < 2:
            raise TopologyResolutionError(
                "A reconstructed canonical atom set must contain at least two atoms."
            )
        atoms = [
            topology.atom_choice(name, f"{topology.source}:canonical-expansion")
            for name in names
        ]
        parents = {atom.parent_atom for atom in atoms}
        if len(parents) != 1 or any(atom.element != "H" for atom in atoms):
            raise TopologyResolutionError(
                f"Canonical atoms {names!r} in {comp_id.upper()} are not one physical "
                "proton set on a single heavy parent."
            )
        return AtomSetChoice(
            atoms=atoms,
            expression=author_expression,
            semantics="canonical-expansion-atom-set",
            warnings=[
                "NMR-STAR canonical OR rows reconstructed as one topology-verified "
                "author-level proton set"
            ],
        )


def component_topology_snapshot(
    records: Iterable[SequenceRecord], library: TopologyLibrary
) -> dict[str, dict[str, object]]:
    """Freeze exact component atom membership used by one conversion.

    Unknown components are represented explicitly with an empty atom set. This
    keeps the output validator fail-closed and makes the evidence auditable.
    """
    snapshot: dict[str, dict[str, object]] = {}
    for comp_id in sorted({record.residue_name.strip().upper() for record in records}):
        topology = library.get(comp_id)
        snapshot[comp_id] = {
            "available": topology is not None,
            "source": topology.source if topology is not None else None,
            "atoms": sorted(topology.available_atoms()) if topology is not None else [],
        }
    return snapshot


def mapped_residue_index(
    records: Iterable[SequenceRecord],
) -> dict[tuple[str, int], str]:
    """Return an unambiguous Boltz position-to-component mapping."""
    mapped: dict[tuple[str, int], str] = {}
    for record in records:
        key = (record.boltz_chain, record.boltz_residue_index)
        comp_id = record.residue_name.strip().upper()
        prior = mapped.get(key)
        if prior is not None and prior != comp_id:
            raise AtomTopologyValidationError(
                f"Conflicting mapped residue identities at {key[0]}:{key[1]}: "
                f"{prior} and {comp_id}."
            )
        mapped[key] = comp_id
    return mapped


def atom_topology_violations(
    atoms: Iterable[BoltzAtom],
    *,
    mapped_residues: dict[tuple[str, int], str],
    component_topologies: dict[str, dict[str, object]],
) -> list[AtomTopologyViolation]:
    """Check atoms against component dictionaries, never coordinate observations."""
    violations: list[AtomTopologyViolation] = []
    for atom in atoms:
        comp_id = mapped_residues.get((atom.chain, atom.residue_index))
        if comp_id is None:
            violations.append(
                AtomTopologyViolation(atom, None, None, "mapped_residue_not_found")
            )
            continue
        evidence = component_topologies.get(comp_id)
        source = str(evidence.get("source")) if evidence and evidence.get("source") else None
        if not evidence or evidence.get("available") is not True:
            violations.append(
                AtomTopologyViolation(atom, comp_id, source, "component_topology_unavailable")
            )
            continue
        available = evidence.get("atoms")
        if not isinstance(available, list) or atom.atom_name not in available:
            violations.append(
                AtomTopologyViolation(atom, comp_id, source, "atom_absent_from_component_topology")
            )
    return violations


def emitted_atom_topology_violations(report: ConversionReport) -> list[AtomTopologyViolation]:
    mapped = mapped_residue_index(report.sequence_map)
    return atom_topology_violations(
        (
            atom
            for constraint in report.emitted_constraints
            for atom in (constraint.atom1, constraint.atom2)
        ),
        mapped_residues=mapped,
        component_topologies=report.target_component_topologies,
    )


def output_atom_topology_violations(report: ConversionReport) -> list[AtomTopologyViolation]:
    """Validate every atom that can be written to an executable constraint file."""
    mapped = mapped_residue_index(report.sequence_map)
    return atom_topology_violations(
        itertools.chain(
            (
                atom
                for constraint in report.emitted_constraints
                for atom in (constraint.atom1, constraint.atom2)
            ),
            (
                atom
                for group in report.ambiguous_groups
                for alternative in group.alternatives
                for atom in (alternative.atom1, alternative.atom2)
            ),
        ),
        mapped_residues=mapped,
        component_topologies=report.target_component_topologies,
    )


def require_valid_emitted_atom_topology(report: ConversionReport) -> None:
    """Fail before executable YAML can contain an unproven atom."""
    violations = emitted_atom_topology_violations(report)
    if not violations:
        return
    details = "; ".join(
        f"{item.atom.display()} mapped to {item.residue_name or '?'} ({item.reason})"
        for item in violations[:20]
    )
    remainder = len(violations) - 20
    if remainder > 0:
        details += f"; and {remainder} additional violation(s)"
    raise AtomTopologyValidationError(
        "Executable atom-contact topology validation failed: " + details
    )


def require_valid_output_atom_topology(report: ConversionReport) -> None:
    """Fail before either exact or union output can contain an unproven atom."""
    violations = output_atom_topology_violations(report)
    if not violations:
        return
    details = "; ".join(
        f"{item.atom.display()} mapped to {item.residue_name or '?'} ({item.reason})"
        for item in violations[:20]
    )
    remainder = len(violations) - 20
    if remainder > 0:
        details += f"; and {remainder} additional violation(s)"
    raise AtomTopologyValidationError(
        "Executable atom-contact topology validation failed: " + details
    )


def _looks_like_hydrogen(atom_name: str) -> bool:
    stripped = atom_name.lstrip("123456789").upper()
    return stripped.startswith(("H", "D", "T"))


def _optional_float(value: str) -> float | None:
    text = clean(value)
    if text is None:
        return None
    try:
        result = float(text)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def build_builtin_topologies() -> dict[str, ComponentTopology]:
    topologies: dict[str, ComponentTopology] = {}

    def component(comp_id: str, source: str = "builtin-standard") -> ComponentTopology:
        topo = ComponentTopology(comp_id=comp_id, source=source)
        topologies[comp_id] = topo
        return topo

    def add(topo: ComponentTopology, parent: str, *hydrogens: str, parent_element: str | None = None) -> None:
        for hydrogen in hydrogens:
            topo.add_hydrogen(hydrogen, parent, parent_element)

    def add_heavy(topo: ComponentTopology, *atom_names: str) -> None:
        for atom_name in atom_names:
            topo.add_atom(atom_name, infer_element(atom_name))

    protein_ids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    ]
    proteins = {comp_id: component(comp_id) for comp_id in protein_ids}
    protein_sidechains = {
        "ALA": ("CB",),
        "ARG": ("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"),
        "ASN": ("CB", "CG", "OD1", "ND2"),
        "ASP": ("CB", "CG", "OD1", "OD2"),
        "CYS": ("CB", "SG"),
        "GLN": ("CB", "CG", "CD", "OE1", "NE2"),
        "GLU": ("CB", "CG", "CD", "OE1", "OE2"),
        "GLY": (),
        "HIS": ("CB", "CG", "ND1", "CD2", "CE1", "NE2"),
        "ILE": ("CB", "CG1", "CG2", "CD1"),
        "LEU": ("CB", "CG", "CD1", "CD2"),
        "LYS": ("CB", "CG", "CD", "CE", "NZ"),
        "MET": ("CB", "CG", "SD", "CE"),
        "PHE": ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
        "PRO": ("CB", "CG", "CD"),
        "SER": ("CB", "OG"),
        "THR": ("CB", "OG1", "CG2"),
        "TRP": ("CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
        "TYR": ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"),
        "VAL": ("CB", "CG1", "CG2"),
    }
    for comp_id, topo in proteins.items():
        add_heavy(topo, "N", "CA", "C", "O", "OXT", *protein_sidechains[comp_id])
        if comp_id != "PRO":
            add(topo, "N", "H")
            topo.add_alias("HN", "H")
        add(topo, "N", "H1", "H2", "H3")
        topo.add_alias("HT1", "H1")
        topo.add_alias("HT2", "H2")
        topo.add_alias("HT3", "H3")
        if comp_id != "GLY":
            add(topo, "CA", "HA")
            topo.add_alias("HA1", "HA")
        add(topo, "OXT", "HXT", parent_element="O")

    add(proteins["ALA"], "CB", "HB1", "HB2", "HB3")
    add(proteins["ARG"], "CB", "HB2", "HB3")
    add(proteins["ARG"], "CG", "HG2", "HG3")
    add(proteins["ARG"], "CD", "HD2", "HD3")
    add(proteins["ARG"], "NE", "HE", parent_element="N")
    add(proteins["ARG"], "NH1", "HH11", "HH12", parent_element="N")
    add(proteins["ARG"], "NH2", "HH21", "HH22", parent_element="N")
    add(proteins["ASN"], "CB", "HB2", "HB3")
    add(proteins["ASN"], "ND2", "HD21", "HD22", parent_element="N")
    add(proteins["ASP"], "CB", "HB2", "HB3")
    add(proteins["ASP"], "OD1", "HD1", parent_element="O")
    add(proteins["ASP"], "OD2", "HD2", parent_element="O")
    add(proteins["CYS"], "CB", "HB2", "HB3")
    add(proteins["CYS"], "SG", "HG", parent_element="S")
    add(proteins["GLN"], "CB", "HB2", "HB3")
    add(proteins["GLN"], "CG", "HG2", "HG3")
    add(proteins["GLN"], "NE2", "HE21", "HE22", parent_element="N")
    add(proteins["GLU"], "CB", "HB2", "HB3")
    add(proteins["GLU"], "CG", "HG2", "HG3")
    add(proteins["GLU"], "OE1", "HE1", parent_element="O")
    add(proteins["GLU"], "OE2", "HE2", parent_element="O")
    add(proteins["GLY"], "CA", "HA2", "HA3")
    add(proteins["HIS"], "CB", "HB2", "HB3")
    add(proteins["HIS"], "ND1", "HD1", parent_element="N")
    add(proteins["HIS"], "CD2", "HD2")
    add(proteins["HIS"], "CE1", "HE1")
    add(proteins["HIS"], "NE2", "HE2", parent_element="N")
    add(proteins["ILE"], "CB", "HB")
    add(proteins["ILE"], "CG1", "HG12", "HG13")
    add(proteins["ILE"], "CG2", "HG21", "HG22", "HG23")
    add(proteins["ILE"], "CD1", "HD11", "HD12", "HD13")
    add(proteins["LEU"], "CB", "HB2", "HB3")
    add(proteins["LEU"], "CG", "HG")
    add(proteins["LEU"], "CD1", "HD11", "HD12", "HD13")
    add(proteins["LEU"], "CD2", "HD21", "HD22", "HD23")
    add(proteins["LYS"], "CB", "HB2", "HB3")
    add(proteins["LYS"], "CG", "HG2", "HG3")
    add(proteins["LYS"], "CD", "HD2", "HD3")
    add(proteins["LYS"], "CE", "HE2", "HE3")
    add(proteins["LYS"], "NZ", "HZ1", "HZ2", "HZ3", parent_element="N")
    add(proteins["MET"], "CB", "HB2", "HB3")
    add(proteins["MET"], "CG", "HG2", "HG3")
    add(proteins["MET"], "CE", "HE1", "HE2", "HE3")
    add(proteins["PHE"], "CB", "HB2", "HB3")
    add(proteins["PHE"], "CD1", "HD1")
    add(proteins["PHE"], "CD2", "HD2")
    add(proteins["PHE"], "CE1", "HE1")
    add(proteins["PHE"], "CE2", "HE2")
    add(proteins["PHE"], "CZ", "HZ")
    add(proteins["PRO"], "N", "H", "H2", "H3", parent_element="N")
    add(proteins["PRO"], "CB", "HB2", "HB3")
    add(proteins["PRO"], "CG", "HG2", "HG3")
    add(proteins["PRO"], "CD", "HD2", "HD3")
    add(proteins["SER"], "CB", "HB2", "HB3")
    add(proteins["SER"], "OG", "HG", parent_element="O")
    add(proteins["THR"], "CB", "HB")
    add(proteins["THR"], "OG1", "HG1", parent_element="O")
    add(proteins["THR"], "CG2", "HG21", "HG22", "HG23")
    add(proteins["TRP"], "CB", "HB2", "HB3")
    add(proteins["TRP"], "CD1", "HD1")
    add(proteins["TRP"], "NE1", "HE1", parent_element="N")
    add(proteins["TRP"], "CE3", "HE3")
    add(proteins["TRP"], "CZ2", "HZ2")
    add(proteins["TRP"], "CZ3", "HZ3")
    add(proteins["TRP"], "CH2", "HH2")
    add(proteins["TYR"], "CB", "HB2", "HB3")
    add(proteins["TYR"], "CD1", "HD1")
    add(proteins["TYR"], "CD2", "HD2")
    add(proteins["TYR"], "CE1", "HE1")
    add(proteins["TYR"], "CE2", "HE2")
    add(proteins["TYR"], "OH", "HH", parent_element="O")
    add(proteins["VAL"], "CB", "HB")
    add(proteins["VAL"], "CG1", "HG11", "HG12", "HG13")
    add(proteins["VAL"], "CG2", "HG21", "HG22", "HG23")

    # Common PDB protonation/variant aliases.
    variant_map = {
        "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
        "ASH": "ASP", "GLH": "GLU", "LYN": "LYS",
        "CYX": "CYS", "CYM": "CYS",
    }
    for variant, parent in variant_map.items():
        base = topologies[parent]
        topologies[variant] = ComponentTopology(
            comp_id=variant,
            atom_elements=dict(base.atom_elements),
            hydrogen_parent=dict(base.hydrogen_parent),
            hydrogen_bond_upper=dict(base.hydrogen_bond_upper),
            atom_aliases=dict(base.atom_aliases),
            source=f"builtin-variant-of-{parent}",
        )

    mse = ComponentTopology(
        comp_id="MSE",
        atom_elements=dict(topologies["MET"].atom_elements),
        hydrogen_parent=dict(topologies["MET"].hydrogen_parent),
        hydrogen_bond_upper=dict(topologies["MET"].hydrogen_bond_upper),
        atom_aliases=dict(topologies["MET"].atom_aliases),
        source="builtin-MSE",
    )
    mse.atom_elements.pop("SD", None)
    mse.atom_elements["SE"] = "SE"
    topologies["MSE"] = mse
    sec = ComponentTopology(
        comp_id="SEC",
        atom_elements=dict(topologies["CYS"].atom_elements),
        hydrogen_parent=dict(topologies["CYS"].hydrogen_parent),
        hydrogen_bond_upper=dict(topologies["CYS"].hydrogen_bond_upper),
        atom_aliases=dict(topologies["CYS"].atom_aliases),
        source="builtin-SEC",
    )
    sec.atom_elements.pop("SG", None)
    sec.atom_elements["SE"] = "SE"
    sec.hydrogen_parent["HG"] = "SE"
    sec.hydrogen_bond_upper["HG"] = DEFAULT_XH_UPPER["SE"]
    topologies["SEC"] = sec

    def nucleic(comp_id: str, *, deoxy: bool = False) -> ComponentTopology:
        topo = component(comp_id, source="builtin-nucleic-acid")
        add_heavy(
            topo,
            "P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "O4'",
            "C3'", "O3'", "C2'", "C1'",
        )
        add(topo, "C1'", "H1'")
        add(topo, "C2'", "H2'", "H2''", "H2'1", "H2'2")
        if not deoxy:
            add_heavy(topo, "O2'")
            add(topo, "O2'", "HO2'", parent_element="O")
        add(topo, "C3'", "H3'")
        add(topo, "O3'", "HO3'", parent_element="O")
        add(topo, "C4'", "H4'")
        add(topo, "C5'", "H5'", "H5''", "H5'1", "H5'2")
        add(topo, "O5'", "HO5'", parent_element="O")
        add(topo, "OP1", "HOP1", parent_element="O")
        add(topo, "OP2", "HOP2", parent_element="O")
        add(topo, "OP3", "HOP3", parent_element="O")
        return topo

    adenines = ["A", "DA", "ADE", "AMP"]
    cytosines = ["C", "DC", "CYT", "CMP"]
    guanines = ["G", "DG", "GUA", "GMP"]
    uracils = ["U", "URA", "UMP"]
    thymines = ["T", "DT", "THY", "TMP"]
    inosines = ["I", "DI", "INO", "IMP"]
    for comp_id in adenines:
        topo = nucleic(comp_id, deoxy=comp_id == "DA")
        add_heavy(topo, "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4")
        add(topo, "C2", "H2")
        add(topo, "C8", "H8")
        add(topo, "N6", "H61", "H62", parent_element="N")
    for comp_id in cytosines:
        topo = nucleic(comp_id, deoxy=comp_id == "DC")
        add_heavy(topo, "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6")
        add(topo, "C5", "H5")
        add(topo, "C6", "H6")
        add(topo, "N4", "H41", "H42", parent_element="N")
    for comp_id in guanines:
        topo = nucleic(comp_id, deoxy=comp_id == "DG")
        add_heavy(topo, "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4")
        add(topo, "N1", "H1", parent_element="N")
        add(topo, "C8", "H8")
        add(topo, "N2", "H21", "H22", parent_element="N")
    for comp_id in uracils:
        topo = nucleic(comp_id)
        add_heavy(topo, "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6")
        add(topo, "N3", "H3", parent_element="N")
        add(topo, "C5", "H5")
        add(topo, "C6", "H6")
    for comp_id in thymines:
        topo = nucleic(comp_id, deoxy=comp_id == "DT")
        add_heavy(topo, "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "C7", "C5M")
        add(topo, "N3", "H3", parent_element="N")
        add(topo, "C6", "H6")
        add(topo, "C7", "H71", "H72", "H73")
        add(topo, "C5M", "H51", "H52", "H53")
    for comp_id in inosines:
        topo = nucleic(comp_id, deoxy=comp_id == "DI")
        add_heavy(topo, "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N3", "C4")
        add(topo, "N1", "H1", parent_element="N")
        add(topo, "C2", "H2")
        add(topo, "C8", "H8")

    return topologies
