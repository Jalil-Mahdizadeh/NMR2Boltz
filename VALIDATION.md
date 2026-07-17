# Validation status

Status: PASS WITH EXACT REVIEWED CORPUS LIMITATIONS
Date: 2026-07-17

## Regression and stress validation

- 55 Pytest tests passed.
- 100,000 randomized sum-r6 implication cases passed.
- 100,000 constructive triangle-inequality cases passed.
- 25,000 outward-rounding cases passed.
- 10,000 randomized OR-max/AND-min merge-order cases passed.
- All 850 built-in hydrogen-parent mappings resolved.
- Deterministic NEF, NMR-STAR, compressed-input, custom-component, averaging-policy, and 32-hypothesis paths passed.
- Target-YAML validation checks chain existence, one-based index bounds, canonical residue identity, declared modifications, mapping collisions, and emitted-contact positions.
- Sequence-only NEF/NMR-STAR entries produce auditable empty distance conversions.
- Sequence/residue conflicts are rejected explicitly before atom-topology resolution.
- PDB coordinates are aligned to one-based Boltz sequence indices before distance evaluation.
- NMR-STAR sequence-alias provenance is sorted deterministically rather than depending on Python set iteration.
- Sequence-map ordering remains deterministic across Python hash seeds.
- Nucleotide x/y expressions consume only stereochemical numeric suffixes;
  `H4y` no longer matches the literal prime in `H4'`.
- Every discrepancy allowlist predicate has positive and adversarial tests, and
  equivalent source evidence with a changed projection reaches the explicit
  `parser_projection_bug` gate path.

## Paired-format deposited-data validation

Twelve deposited NMR ensembles were converted from both NEF and NMR-STAR.
All 24 conversions completed; the two 8S8O inputs are valid empty distance
conversions containing sequence, shifts, and torsion data but no distance loop.
The default safe projection emitted 13,010 NEF and 11,841 NMR-STAR contacts.
Resolved contact/model satisfaction was 99.88% and 99.86%, respectively.
Across 390,930 model/contact cases with satisfied source antecedents, the
projected heavy-atom implication had zero failures.

The paired audit contains 4,195 discrepancies: 4,136 match an explicit tested
semantic predicate and 59 are deposition inconsistencies. There are zero
unresolved discrepancies and zero remaining parser/projection bugs. The audit
digest is
`f6b32f8fd521da9371ce08fdc4828550c078161deda41ae108033b39609aa244`.

The exact reviewed missing-coordinate set contains 106 format/contact records:
45 each for NEF and STAR in partial-coordinate 8R1X, and 8 each for NEF and
STAR in the 9CCH `ZN`/`ZN*` deposition inconsistency. Every identity, bound,
and source group is stored in `benchmark/reviewed_baseline.json`; the set digest
is `60ed1b785e328fb2d972de6f5d691765012eb622afa489b7a47f41909c8eb141`.
The CI gate passes only for this exact set and fails on additions, removals, or
changes.

Only three positive-distance cases have exact NEF/STAR pair-and-bound parity
(21CC, 9D99, and 9KG4); 8S8O also has exact empty-output parity. Large differences
caused by pseudoatoms, atom-set multiplicity, and source inconsistencies remain
visible and are not silently approximated. The complete current report is under
`benchmark/output`. This benchmark validates conversion behavior, not the
predictive accuracy of a GPU Boltz folding campaign.

## Container status

- Image: `nmr2boltz:0.1.0-validated`
- Image ID / repository digest:
  `sha256:366529d372944a278b74fa99c7e499c6f5aa73836dd8e60ed75513e91c0bb000`
- Reported size: 243 MB; default user: `65532:65532`; other nmr2boltz images: 0.
- An offline (`--network none`), read-only-root smoke test ran as UID/GID
  `65532:65532` and resolved nucleotide `H4y` only to `H41/H42` on parent N4.
- An offline, read-only-root 6M6O NMR-STAR conversion emitted 2,877 contacts,
  7 quarantined OR groups, and 37 rejections. Its
  `heavy_atom_constraints.tsv` and `rejections.tsv` were byte-identical to the
  committed benchmark artifacts (SHA-256
  `ea1b6f17490457f860de36ae4a3854604c9ff8e7ad974da88c6004b3589b989a`
  and `6e0a5981c748369e1fd492a399510e7614e7698ee97e588719ecc8027bb6e76c`).
