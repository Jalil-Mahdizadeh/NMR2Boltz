# nmr2boltz toolkit release notes

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

See `VALIDATION.md` and `workspace/output/REAL_TEST_6M6O.md`. The current
validation covers regression tests, deterministic randomized stress tests,
container isolation, and a 10-model deposited-data benchmark for PDB 6M6O.
