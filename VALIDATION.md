# Validation status

Status: PASS WITH EXACT REVIEWED CORPUS LIMITATIONS
Date: 2026-07-20

## Regression and stress validation

- 138 Pytest tests passed.
- 100,000 randomized sum-r6 implication cases passed.
- 100,000 constructive triangle-inequality cases passed.
- 25,000 outward-rounding cases passed.
- 10,000 randomized OR-max/AND-min merge-order cases passed.
- All 845 built-in hydrogen-parent mappings resolved.
- Deterministic NEF, NMR-STAR, compressed-input, custom-component, averaging-policy, and 32-hypothesis paths passed.
- Target-YAML validation checks chain existence, one-based index bounds,
  canonical residue identity, declared modifications, mapping collisions, and
  every exact and union-alternative position. Its report counts exact
  constraints, union groups, and union alternatives explicitly.
- Sequence-only NEF/NMR-STAR entries produce auditable empty distance conversions.
- Sequence/residue conflicts are rejected explicitly before atom-topology resolution.
- Distinct source residues cannot share one mapped Boltz chain/index, including
  same-component collisions; both projection and target validation reject them.
- Missing-chain author/canonical sequence identifiers are evaluated in stable
  order and conflicting unique resolutions fail closed.
- PDB coordinates are aligned to one-based Boltz sequence indices before distance evaluation.
- NMR-STAR sequence-alias provenance is sorted deterministically rather than depending on Python set iteration.
- Distinct NMR-STAR canonical atom spellings that topology proves are aliases
  of the same physical proton are exercised on both endpoints; the rows retain
  their source identities and deduplicate only after projection.
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
- Author and canonical component identities are frozen together, so a
  topology-verified canonical alias validates correctly without accepting an
  atom absent from every declared component.
- CCD hydrogens with multiple heavy parents are rejected deterministically,
  independent of bond-row order; malformed CCD atom/bond tables raise
  contextual errors instead of being silently skipped.
- Exact and ambiguous constraints are serialized into separate metadata-free
  files, with deterministic ordering and conservative six-decimal formatting.
- Exact and union bounds share the executable interval policy: sub-minimum
  values are weakened to the configured minimum, while any over-maximum union
  alternative quarantines its complete OR group. The writer independently
  rejects non-finite or out-of-interval exact and union bounds.
- `--exclude-intrachain` was exercised end to end for NEF and NMR-STAR
  protein, DNA, and RNA chains. The tests cover exact contacts, all-inter-chain
  unions, mixed-scope OR quarantine, mapped-chain identity, empty constraint
  files, strict status, benchmark manifests, and the final writer invariant.
- `--exclude-intraresidue` is enforced before the exact/union split. Regression
  tests cover exact contacts, all-intraresidue unions, mixed
  intraresidue/inter-residue unions, unaffected inter-residue unions, CLI
  output, deterministic provenance and counters, and the final writer
  invariant.
- The version-controlled CLI reference is checked against the live argparse
  tree so every global, `convert`, and `benchmark` flag has a documented row
  and default value.
- Invalid contacts are quarantined deterministically with source rows, mapped
  residue/component identity, atom, restraint group, and original bounds;
  coordinate absence is not used as topology evidence.
- Output artifacts are written as one staged bundle. Simulated staging-write
  and commit failures preserve the complete prior directory and leave no
  staging/backup residue.

## Paired-format deposited-data validation

Twelve deposited NMR ensembles were converted from both NEF and NMR-STAR.
All 24 conversions completed; the two 8S8O inputs are valid empty distance
conversions containing sequence, shifts, and torsion data but no distance loop.
The default safe projection emitted 12,998 NEF and 11,829 NMR-STAR contacts.
Resolved contact/model satisfaction was 99.88% and 99.86%, respectively.

NMR-STAR canonical OR expansions are now reconstructed before projection when
author endpoint identity, bounds/weights/logic, component topology, and the
complete canonical Cartesian product agree. In 43JX this changed 591 emitted
STAR bounds without changing the 1,716 emitted heavy-pair identities. Restraint
468 now preserves STAR rows 796--798 as one three-proton set (`N=3`) and matches
the NEF projected bound at 7.428047646 A. Observed reconstructed multiplicities
cover `N=2`, `3`, `4`, `6`, and `9`; incomplete or unverifiable candidates are
rejected rather than guessed.

Across 379,449 model/contact cases with satisfied source antecedents, the
projected heavy-atom implication had zero failures.

The paired audit contains 4,177 discrepancies: 4,134 match an explicit tested
semantic predicate and 43 are deposition inconsistencies. There are zero
unresolved discrepancies and zero remaining parser/projection bugs. The audit
digest is
`e710f092a339fcf8d5d7a57d35207cc4b05bea25fe92c0c30b1f97815605c677`.

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
caused by geometric pseudoatoms, genuinely different atom-set/assignment
semantics, and source inconsistencies remain visible and are not silently
approximated. The complete current report is under
`benchmark/output`. This benchmark validates conversion behavior, not the
predictive accuracy of a GPU Boltz folding campaign.

## Exact-contact distance matrices

`benchmark/distance_check` contains one CSV for each of the 12 PDB IDs. The
files cover the union of 13,275 exact heavy-atom pairs emitted by NEF and
NMR-STAR and contain 227,920 requested pair/model distance cells. Of these,
227,100 resolve from deposited coordinates. The remaining 820 blank cells are
exactly 41 pairs absent from every model of the partial-coordinate 8R1X
ensemble; no other entry has a missing distance cell.

The generator independently aligns PDB author numbering through each format's
sequence map. A common pair is accepted only if NEF and NMR-STAR resolve it in
the same models and their computed distances agree within 1e-6 A. It also
requires the exact executable YAML pair set and conservatively rounded bound to
match `conversion_report.json`. Two consecutive generations produced
byte-identical files, and every CSV digest is pinned in
`benchmark/distance_check/distance_check_summary.json`.

The CSVs intentionally exclude ambiguous `atom_contact_union` groups because a
row identified only by an atom pair cannot preserve restraint-group OR
semantics. Their alternatives remain available in the normal union and
provenance artifacts. Descriptive satisfaction fractions in the distance-check
README compare deposited structures with source-driven exact bounds; they do
not establish predictive accuracy.

## Container status

- Image: `nmr2boltz:0.1.0-validated`
- Image ID / repository digest:
  `sha256:dd488022ae6a425812f0268fb33e2abf5003e25278c62b408a95ed9b85d01882`
- Reported size: 243 MB; default user: `65532:65532`; other nmr2boltz images: 0.
- The rebuilt image passed the 100,000-case `sum-r6` and triangle-inequality
  stress checks, 25,000 outward-rounding checks, 10,000 merge-order checks,
  all 845 built-in hydrogen-parent mappings, and deterministic 32-hypothesis
  generation while offline, read-only, and non-root.
- An offline, read-only-root, non-root conversion using the documented writable
  `/workspace` parent mount emitted 3 exact constraints and 1 union group into
  a complete 10-file staged bundle.
- Offline, read-only-root, non-root NEF and NMR-STAR smoke runs of
  `--exclude-intrachain` each emitted 2 exact contacts and 1 all-inter-chain
  union, filtered 3 restraint groups including 1 mixed-scope OR group, and
  contained no same-chain YAML alternative. The two formats produced
  byte-identical exact and union YAML.
- An offline (`--network none`), read-only-root smoke test ran as UID/GID
  `65532:65532` and proved the standard GLN/LEU/TRP invalid-atom exclusions
  alongside valid adenine N6 and guanine O6 membership.
- An offline, read-only-root 43JX NMR-STAR conversion emitted 1,716 exact
  contacts, 16 atom-contact union groups, and 38 rejections. Restraint 468
  preserved rows 796--798, used `N=3`, and projected to 7.428047646 A. Its
  `heavy_atom_constraints.tsv` was byte-identical to the committed benchmark
  artifact (SHA-256
  `017eec00dfaafc8655e2c5867a4d81951fffd077647f33dc5f11a01e72c67c84`).
- An offline, read-only-root 6M6O NMR-STAR conversion emitted 2,877 exact
  contacts, 7 atom-contact union groups, and 37 rejections. Its
  `heavy_atom_constraints.tsv` and `rejections.tsv` were byte-identical to the
  committed benchmark artifacts (SHA-256
  `ea1b6f17490457f860de36ae4a3854604c9ff8e7ad974da88c6004b3589b989a`
  and `16a0616e27ee91eab290c00f61f2aa084980492ca6411f7aa7aaae44eb809f2d`).
- An offline, read-only-root 9CCH NEF conversion emitted 3,350 contacts,
  quarantined exactly 8 GLN B48/ZN rows, and wrote zero matching invalid atoms
  to `atom_constraints_exact.yaml` or `atom_constraints_union.yaml`.
