# nmr2boltz toolkit release notes

## Unreleased

- Reject distinct source residues that collide at one Boltz chain/index even
  when their component names match. The same generic injectivity check now
  protects projection, final atom-topology validation, and `--target-yaml`.
- Resolve missing-chain endpoints deterministically. Author and canonical
  sequence identifiers that select different records now fail closed instead
  of depending on `PYTHONHASHSEED`.
- Reject CCD hydrogens bonded to multiple heavy atoms independently of bond-row
  order, and surface malformed CCD atom/bond tables as contextual topology
  errors rather than silently suppressing them.
- Preserve declared author and canonical component identities in the frozen
  topology mapping, allowing topology-verified aliases such as HSD/HIS while
  still requiring the emitted atom to exist in at least one declared
  component.
- Write conversion artifacts through a sibling staging directory and commit
  them as one rollback-safe bundle. Simulated write and commit failures retain
  the complete previous output instead of leaving mixed-generation files.
- Fixed a reachable NMR-STAR canonical-normalization crash caused by an
  undefined collection name when distinct canonical spellings on both
  endpoints resolve through topology aliases to the same physical protons.
  Alias-equivalent rows now remain duplicate OR evidence, retain deterministic
  source provenance, and pass through normal projected-pair deduplication
  without being mistaken for a multi-proton atom-set expansion.
- Fixed `--exclude-intraresidue` at the projected-group boundary. Exact
  contacts retain their established rejection behavior; all-intraresidue and
  mixed intraresidue/inter-residue OR groups are now quarantined in full with
  deterministic `intraresidue_filtered` provenance and counters. A final
  writer invariant prevents local contacts from re-entering either output.
- Audited `--exclude-intrachain` for the analogous defect. It already applies
  before the exact/union split, quarantines mixed groups in full, and validates
  final exact and union output independently; no code correction was needed.
- Applied the Boltz executable distance interval generically to ambiguous
  unions. Sub-minimum alternatives are raised with raw-bound provenance; if
  any alternative exceeds the maximum, the complete OR group is quarantined
  rather than narrowed. The final writer now fails closed on any non-finite or
  out-of-interval exact or union bound.
- Extended `--target-yaml` validation to every union-alternative endpoint and
  added explicit exact-constraint, union-group, and union-alternative counts.
  Updated BoltzUI guidance to state that atom-contact guidance remains active
  independently of optional `--use_potentials` physical steering.
- Added deterministic per-entry exact-contact distance matrices under
  `benchmark/distance_check`. The generator places NEF and NMR-STAR executable
  bounds beside one aligned Euclidean distance per deposited PDB model,
  cross-checks YAML against conversion provenance, fails closed on
  format-dependent mapping, preserves missing coordinate observations as blank
  cells, and records SHA-256 artifact digests.
- Added a complete `docs/CLI_REFERENCE.md` covering every global, `convert`,
  and `benchmark` flag, including aliases, repeatability, effective defaults,
  positional arguments, and exit statuses. Regression coverage now prevents a
  parser flag from being added without updating the reference.
- Added `--exclude-intrachain` for protein, DNA, and RNA complexes. The filter
  uses mapped Boltz chain IDs, applies to both exact and union outputs, removes
  mixed intra/inter OR groups in full, records deterministic
  `intrachain_filtered` provenance and counts, and is independently enforced
  again by the output writer. Benchmark manifests support
  `options: {exclude_intrachain: true}`.
- Added a generic pre-projection NMR-STAR normalization pass for complete,
  topology-verified canonical OR expansions of one author-level proton set.
  One- and two-sided expansions now use the correct `sum-r6` multiplicity
  (`N=3`, `6`, or `9` for the covered methyl/methylene cases), while branches
  with different heavy parents remain union alternatives. Incomplete,
  inconsistent, or unverifiable expansions are rejected in full with source
  rows and canonical expansions preserved.
- Replaced the mixed legacy constraint artifacts with two deterministic,
  metadata-free outputs: `atom_constraints_exact.yaml` for non-ambiguous
  `atom_contact` constraints and `atom_constraints_union.yaml` for ambiguous
  `atom_contact_union` groups. Every alternative retains its own conservatively
  outward-rounded six-decimal bound.
- Added a generic fail-closed atom-membership boundary after residue mapping and
  heavy-atom projection. A contact is quarantined in full as
  `atom_not_present_in_mapped_residue` unless both atoms exist in the mapped
  standard or declared CCD component; the YAML writer independently rechecks
  the frozen target topology before creating files.
- Expanded standard amino-acid and RNA/DNA heavy-atom inventories, retained
  embedded/external CCD support for modified residues, ligands, and ions, and
  added deterministic rejection provenance with source rows, mapped identity,
  restraint group, and original bounds.
- Quarantined 4 malformed contacts per format in 8R1X and 8 GLN/ZN contacts per
  format in 9CCH. No other corpus entry had this defect.
- Added clean polymer-only `sequences.fasta` output for every conversion.
- Replaced permissive paired-format discrepancy classification with tested,
  fail-closed predicates for wildcard/OR, x/y assignment, rejected Q/M
  pseudoatoms, and verified canonical aliases.
- Fixed nucleotide x/y matching so `H4y` cannot consume the literal prime in
  `H4'`; this removes a false 8R1X sugar-atom branch and restores six NEF
  contacts present in the canonical STAR representation.
- Pinned all reviewed coordinate limitations by exact contact identity, bound,
  source provenance, and SHA-256 digest. Added, removed, or changed contacts
  now fail CI while the exact reviewed 8R1X/9CCH set passes.
- Rebuilt and validated `nmr2boltz:0.1.0-validated` as non-root and offline;
  image digest is
  `sha256:dd488022ae6a425812f0268fb33e2abf5003e25278c62b408a95ed9b85d01882`.
- Added a paired NEF/NMR-STAR corpus runner with sequence-aware PDB ensemble
  alignment, exact contact/bound parity metrics, and per-case output directories.
- Sequence/residue conflicts now use the explicit `sequence_residue_mismatch`
  rejection reason before topology resolution.
- NEF/NMR-STAR files with sequence data but no distance loop now produce an
  auditable empty conversion rather than failing format detection.
- Expanded the regression suite from 39 to 138 tests, including positive and
  adversarial checks for every discrepancy predicate and fail-closed gate.
- Added fail-fast `--target-yaml` validation for chain IDs, residue indices,
  canonical residue identities, declared modifications, mapping collisions, and
  exact and union-alternative positions.
- Added the versioned, checksum-aware `nmr2boltz benchmark` command.
- Made NMR-STAR sequence-alias provenance ordering deterministic across Python
  hash seeds and container runs.
- Added protein, DNA, RNA, and modified-polymer target validation.

Release date: 2026-07-14

This research-grade reference implementation converts NEF 1.1 and NMR-STAR distance restraints into heavy-atom upper-bound constraints for the BoltzUI `atom_contact` extension.

## Included

- category-driven NEF and NMR-STAR parsing with PyNMRSTAR and Gemmi backends;
- gzip, bzip2, and xz/lzma input handling;
- standard protein, DNA, and RNA hydrogen-to-parent topology;
- embedded and external CCD/custom-component topology support;
- explicit sequence-code to Boltz residue-index mapping;
- conservative triangle-inequality projection from proton to parent-heavy-atom bounds;
- selectable atom-set semantics: unnormalized sum-r6, normalized mean-r6, or hard-OR;
- first-class handling of restraint OR groups, wildcard atom sets, x/y assignments, and pseudoatoms;
- separate exact and union constraint YAML, compact text, TSV, JSON audit, and optional hypothesis batches;
- Dockerfile, pinned dependency set, synthetic NEF/NMR-STAR fixtures, tests, and an expert-facing methods document.

## Safety defaults

The default mode does not turn unresolved OR alternatives into simultaneous restraints. It does not silently invent an upper bound from a target value. It does not clip a projected bound above the Boltz 20 A ceiling down to 20 A. It does not interpret geometric Q/M pseudoatoms as atom sets unless explicitly requested.

## Validation

See `VALIDATION.md` and `benchmark/output/BENCHMARK_REPORT.md`. The current
validation covers regression tests, deterministic randomized stress tests, and
paired NEF/NMR-STAR coordinate audits for 12 deposited NMR ensembles.
