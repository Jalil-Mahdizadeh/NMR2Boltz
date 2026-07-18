"""Deterministic mathematical and converter stress validation for nmr2boltz.

This script intentionally uses only runtime dependencies so it can be executed
inside the production Docker image against a read-only bind mount of the repo.
"""

from __future__ import annotations

import argparse
import copy
import gzip
import json
import math
import random
import tempfile
from pathlib import Path

import yaml

from nmr2boltz.model import BoltzAtom, ProjectedAlternative
from nmr2boltz.output import _rounded, write_outputs
from nmr2boltz.project import (
    ProjectionSettings,
    _merge_independent_constraints,
    _merge_or_alternatives,
    project_document,
)
from nmr2boltz.star import parse_star_document
from nmr2boltz.topology import TopologyLibrary


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


def _unit_vector(rng: random.Random) -> tuple[float, float, float]:
    while True:
        vector = tuple(rng.uniform(-1.0, 1.0) for _ in range(3))
        norm = math.sqrt(sum(value * value for value in vector))
        if norm > 1e-12:
            return tuple(value / norm for value in vector)  # type: ignore[return-value]


def _add(
    point: tuple[float, float, float],
    direction: tuple[float, float, float],
    length: float,
) -> tuple[float, float, float]:
    return tuple(point[index] + direction[index] * length for index in range(3))  # type: ignore[return-value]


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[index] - b[index]) ** 2 for index in range(3)))


def validate_sum_r6(rng: random.Random, iterations: int) -> None:
    for _ in range(iterations):
        count = rng.randint(1, 128)
        distances = [10 ** rng.uniform(math.log10(1.2), math.log10(40.0)) for _ in range(count)]
        effective = sum(distance ** -6 for distance in distances) ** (-1.0 / 6.0)
        upper = effective * rng.uniform(1.0, 3.0)
        projected = count ** (1.0 / 6.0) * upper
        assert effective <= upper * (1.0 + 1e-12)
        assert min(distances) <= projected * (1.0 + 1e-12)


def validate_triangle_projection(rng: random.Random, iterations: int) -> None:
    for _ in range(iterations):
        parent1 = tuple(rng.uniform(-100.0, 100.0) for _ in range(3))
        length1 = rng.uniform(0.0, 1.55)
        length2 = rng.uniform(0.0, 1.55)
        proton_upper = rng.uniform(1.5, 20.0)
        proton1 = _add(parent1, _unit_vector(rng), length1)
        proton2 = _add(proton1, _unit_vector(rng), rng.uniform(0.0, proton_upper))
        parent2 = _add(proton2, _unit_vector(rng), length2)
        assert _distance(proton1, proton2) <= proton_upper * (1.0 + 1e-12)
        assert _distance(parent1, parent2) <= proton_upper + length1 + length2 + 1e-10


def validate_outward_rounding(rng: random.Random, iterations: int) -> None:
    for _ in range(iterations):
        value = rng.uniform(2.0, 20.0)
        rounded = _rounded(value)
        assert rounded >= value
        assert rounded - value < 1.000001e-6


def _projected(group: str, distance: float) -> ProjectedAlternative:
    return ProjectedAlternative(
        atom1=BoltzAtom("A", 1, "CA"),
        atom2=BoltzAtom("A", 9, "CB"),
        max_distance=distance,
        source_upper_bound=distance,
        averaging_policy="hard-or",
        averaging_factor=1.0,
        explicit_pair_count=1,
        bond_offset=0.0,
        group_id=group,
    )


def validate_boolean_merge_order(rng: random.Random, iterations: int) -> None:
    settings = ProjectionSettings()
    for case in range(iterations):
        values = [rng.uniform(2.0, 20.0) for _ in range(rng.randint(2, 16))]
        order = list(range(len(values)))
        rng.shuffle(order)
        or_items = [_projected("or-group", values[index]) for index in order]
        merged_or = _merge_or_alternatives(or_items)
        assert len(merged_or) == 1
        assert math.isclose(merged_or[0].max_distance, max(values), rel_tol=0.0, abs_tol=1e-12)

        and_items = [_projected(f"and-{case}-{index}", values[index]) for index in order]
        merged_and, rejected = _merge_independent_constraints(and_items, settings)
        assert not rejected and len(merged_and) == 1
        assert math.isclose(merged_and[0].max_distance, min(values), rel_tol=0.0, abs_tol=1e-12)


def validate_builtin_topologies() -> int:
    library = TopologyLibrary()
    resolved = 0
    for comp_id, topology in sorted(library.components.items()):
        for hydrogen, expected_parent in sorted(topology.hydrogen_parent.items()):
            choices = library.resolve_expression(comp_id, hydrogen)
            atoms = [atom for choice in choices for atom in choice.atoms]
            assert len(atoms) == 1
            assert atoms[0].parent_atom == expected_parent
            assert math.isfinite(atoms[0].bond_length_upper)
            assert atoms[0].bond_length_upper > 0.0
            resolved += 1
    assert resolved > 100
    return resolved


def _convert(path: Path, settings: ProjectionSettings | None = None):
    parsed = parse_star_document(path)
    report = project_document(
        parsed,
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=settings or ProjectionSettings(),
    )
    assert report.source_restraint_groups
    pairs = [constraint.pair_key for constraint in report.emitted_constraints]
    assert len(pairs) == len(set(pairs))
    for constraint in report.emitted_constraints:
        assert math.isfinite(constraint.max_distance)
        assert 2.0 <= constraint.max_distance <= 20.0
        assert constraint.provenance
        assert "source_observation" in constraint.provenance[0]
    return report


def validate_end_to_end() -> dict[str, int]:
    reports = {
        "nef": _convert(FIXTURES / "example.nef"),
        "nmrstar": _convert(FIXTURES / "example.str"),
        "custom_component": _convert(FIXTURES / "custom_component.nef"),
    }
    for policy in ("sum-r6", "mean-r6", "hard-or"):
        _convert(FIXTURES / "example.nef", ProjectionSettings(averaging_policy=policy))

    with tempfile.TemporaryDirectory(prefix="nmr2boltz-stress-") as temporary:
        root = Path(temporary)
        compressed = root / "example.nef.gz"
        with (FIXTURES / "example.nef").open("rb") as source, gzip.open(compressed, "wb") as target:
            target.write(source.read())
        _convert(compressed)

        first = root / "first"
        second = root / "second"
        write_outputs(copy.deepcopy(reports["nef"]), first, hypothesis_count=32, random_seed=20260716)
        write_outputs(copy.deepcopy(reports["nef"]), second, hypothesis_count=32, random_seed=20260716)
        relative_files = sorted(path.relative_to(first) for path in first.rglob("*") if path.is_file())
        assert relative_files
        for relative in relative_files:
            assert (first / relative).read_bytes() == (second / relative).read_bytes()
        payload = yaml.safe_load(
            (first / "atom_constraints_exact.yaml").read_text(encoding="utf-8")
        )
        for wrapped in payload["constraints"]:
            contact = wrapped["atom_contact"]
            assert contact["force"] is True
            assert 2.0 <= contact["max_distance"] <= 20.0
            assert len(contact["atom1"]) == len(contact["atom2"]) == 3

    return {
        name: len(report.emitted_constraints)
        for name, report in reports.items()
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20260716)
    args = parser.parse_args()
    if args.iterations < 1:
        parser.error("--iterations must be positive")
    rng = random.Random(args.seed)

    validate_sum_r6(rng, args.iterations)
    validate_triangle_projection(rng, args.iterations)
    validate_outward_rounding(rng, max(1_000, args.iterations // 4))
    validate_boolean_merge_order(rng, max(1_000, args.iterations // 10))
    topology_atoms = validate_builtin_topologies()
    emitted = validate_end_to_end()

    print(json.dumps({
        "status": "passed",
        "seed": args.seed,
        "sum_r6_cases": args.iterations,
        "triangle_cases": args.iterations,
        "rounding_cases": max(1_000, args.iterations // 4),
        "boolean_merge_cases": max(1_000, args.iterations // 10),
        "builtin_hydrogens_resolved": topology_atoms,
        "emitted_constraints_by_fixture": emitted,
        "deterministic_hypotheses": 32,
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
