# nmr2boltz toolkit release notes

## Unreleased

- Added clean polymer-only `sequences.fasta` output for every conversion.
- Added a paired NEF/NMR-STAR corpus runner with sequence-aware PDB ensemble
  alignment, exact contact/bound parity metrics, and per-case output directories.
- Sequence/residue conflicts now use the explicit `sequence_residue_mismatch`
  rejection reason before topology resolution.
- NEF/NMR-STAR files with sequence data but no distance loop now produce an
  auditable empty conversion rather than failing format detection.
- Expanded the regression suite from 39 to 44 tests.
- Added fail-fast `--target-yaml` validation for chain IDs, residue indices,
  canonical residue identities, declared modifications, mapping collisions, and
  emitted-contact positions.
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
- safe current-Boltz YAML, compact text, TSV, JSON audit, proposed union restraints, and optional hypothesis batches;
- Dockerfile, pinned dependency set, synthetic NEF/NMR-STAR fixtures, tests, and an expert-facing methods document.

## Safety defaults

The default mode does not turn unresolved OR alternatives into simultaneous restraints. It does not silently invent an upper bound from a target value. It does not clip a projected bound above the Boltz 20 A ceiling down to 20 A. It does not interpret geometric Q/M pseudoatoms as atom sets unless explicitly requested.

## Validation

See `VALIDATION.md` and `benchmark/output/BENCHMARK_REPORT.md`. The current
validation covers regression tests, deterministic randomized stress tests, and
paired NEF/NMR-STAR coordinate audits for 12 deposited NMR ensembles.
