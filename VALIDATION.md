# Validation status

Status: PASS WITH EXACT REVIEWED CORPUS LIMITATIONS
Date: 2026-07-18

## Regression and stress validation

- 79 Pytest tests passed.
- 100,000 randomized sum-r6 implication cases passed.
- 100,000 constructive triangle-inequality cases passed.
- 25,000 outward-rounding cases passed.
- 10,000 randomized OR-max/AND-min merge-order cases passed.
- All 845 built-in hydrogen-parent mappings resolved.
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
- Both projected endpoints are checked against the mapped standard or declared
  CCD component before deduplication, and the output writer independently
  validates every exact and union endpoint against the frozen target-topology
  snapshot.
- Exact and ambiguous constraints are serialized into separate metadata-free
  files, with deterministic ordering and conservative six-decimal formatting.
- Invalid contacts are quarantined deterministically with source rows, mapped
  residue/component identity, atom, restraint group, and original bounds;
  coordinate absence is not used as topology evidence.

## Paired-format deposited-data validation

Twelve deposited NMR ensembles were converted from both NEF and NMR-STAR.
All 24 conversions completed; the two 8S8O inputs are valid empty distance
conversions containing sequence, shifts, and torsion data but no distance loop.
The default safe projection emitted 12,998 NEF and 11,829 NMR-STAR contacts.
Resolved contact/model satisfaction was 99.88% and 99.86%, respectively.
Across 390,930 model/contact cases with satisfied source antecedents, the
projected heavy-atom implication had zero failures.

The paired audit contains 4,179 discrepancies: 4,136 match an explicit tested
semantic predicate and 43 are deposition inconsistencies. There are zero
unresolved discrepancies and zero remaining parser/projection bugs. The audit
digest is
`3b4e4f3d54377d7f80854a1525db0f5351895fbb5c3b6bf2e63f2e3a649c86c8`.

The exact reviewed missing-coordinate set contains 82 format/contact records:
41 each for NEF and STAR in partial-coordinate 8R1X. Four malformed 8R1X
contacts and all eight 9CCH GLN/ZN contacts per format are now quarantined
before coordinate evaluation. Every remaining identity, bound, and source group
is stored in `benchmark/reviewed_baseline.json`; the set digest is
`ce8bc2710b5068d3d1cec33794f77bb4e806019a2116a9f4127102b25cfb63e2`.
The CI gate passes only for this exact set and fails on additions, removals, or
changes. It also fails on any emitted atom-topology violation. The current
corpus has zero such violations.

Only three positive-distance cases have exact NEF/STAR pair-and-bound parity
(21CC, 9D99, and 9KG4); 8S8O also has exact empty-output parity. Large differences
caused by pseudoatoms, atom-set multiplicity, and source inconsistencies remain
visible and are not silently approximated. The complete current report is under
`benchmark/output`. This benchmark validates conversion behavior, not the
predictive accuracy of a GPU Boltz folding campaign.

## Container status

- Image: `nmr2boltz:0.1.0-validated`
- Image ID / repository digest:
  `sha256:e1b8da2544111093b21f952077c9a601233d46d618c9b390b7f9ae4b835d7dc2`
- Reported size: 243 MB; default user: `65532:65532`; other nmr2boltz images: 0.
- An offline (`--network none`), read-only-root smoke test ran as UID/GID
  `65532:65532` and proved the standard GLN/LEU/TRP invalid-atom exclusions
  alongside valid adenine N6 and guanine O6 membership.
- An offline, read-only-root 6M6O NMR-STAR conversion emitted 2,877 exact
  contacts, 7 atom-contact union groups, and 37 rejections. Its
  `heavy_atom_constraints.tsv` and `rejections.tsv` were byte-identical to the
  committed benchmark artifacts (SHA-256
  `ea1b6f17490457f860de36ae4a3854604c9ff8e7ad974da88c6004b3589b989a`
  and `16a0616e27ee91eab290c00f61f2aa084980492ca6411f7aa7aaae44eb809f2d`).
- An offline, read-only-root 9CCH NEF conversion emitted 3,350 contacts,
  quarantined exactly 8 GLN B48/ZN rows, and wrote zero matching invalid atoms
  to `atom_constraints_exact.yaml` or `atom_constraints_union.yaml`.
