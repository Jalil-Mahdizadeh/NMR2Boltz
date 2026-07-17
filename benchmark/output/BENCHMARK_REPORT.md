# nmr2boltz paired-format benchmark

All 12 deposited ensembles were converted from both NEF and NMR-STAR with conservative defaults, then audited against every PDB conformer using sequence-aware coordinate alignment.

| Case | NEF contacts | STAR contacts | NEF PDB satisfaction | STAR PDB satisfaction | Exact parity | Implication failures |
|---|---:|---:|---:|---:|---:|---:|
| 21CC | 176 | 176 | 99.72% | 99.72% | yes | 0 |
| 43JX | 1716 | 1716 | 99.96% | 99.92% | no | 0 |
| 6M6O | 2819 | 2877 | 99.72% | 99.60% | no | 0 |
| 8R1X | 1118 | 1118 | 99.74% | 99.74% | no | 0 |
| 8S8O | 0 | 0 | N/A | N/A | yes | 0 |
| 9CCH | 3358 | 3524 | 99.98% | 100.00% | no | 0 |
| 9D99 | 456 | 456 | 100.00% | 100.00% | yes | 0 |
| 9KG4 | 406 | 406 | 99.74% | 99.74% | yes | 0 |
| 9PQH | 3 | 44 | 3.33% | 100.00% | no | 0 |
| 9SGX | 1578 | 765 | 100.00% | 100.00% | no | 0 |
| 9VQ1 | 53 | 56 | 95.00% | 94.20% | no | 0 |
| 9VUY | 1327 | 703 | 99.96% | 100.00% | no | 0 |

## Current result

- 24/24 conversions completed; 2 are valid empty distance conversions for 8S8O.
- NEF: 13010 contacts, 99.88% resolved PDB satisfaction.
- NMR-STAR: 11841 contacts, 99.86% resolved PDB satisfaction.
- Conservative implication failures: 0 across 390930 satisfied-antecedent cases.
- Exact NEF/STAR pair-and-bound parity: 3/11 positive-distance cases; 8S8O also has exact empty-output parity.

## Row-level format discrepancy audit

- 4195 NEF-only, STAR-only, or different-bound contacts were audited to source rows and physical proton sets.
- Classifications: 4136 scientifically expected format differences; 59 deposition inconsistencies; 0 unresolved; 0 parser/projection bugs.
- Reviewed audit digest: `f6b32f8fd521da9371ce08fdc4828550c078161deda41ae108033b39609aa244`.
- `expected_format_difference` is allowlisted only for wildcard-set versus explicit OR, x/y assignment versus a compatible physical set, rejected Q/M pseudoatoms, or verified canonical aliases.
- 9PQH contains 43 contact discrepancies caused by its NEF sequence/residue numbering conflict; its remaining two different-bound rows are the expected explicit-OR versus wildcard atom-set distinction.
- 9CCH contains 16 deposition inconsistencies because `ZN` and `ZN*` are unverified heavy-atom names that topology cannot prove equivalent.
- 43JX, 6M6O, 9SGX, and 9VUY differences are justified by explicit-OR, wildcard atom-set, x/y assignment, or rejected geometric-pseudoatom semantics.

## Fail-closed gate

- Gate status: **PASS** (0 failure categories).
- Coordinate resolution gaps: 53 NEF and 53 STAR contacts. These are not omitted from the denominator silently.
- The current gaps are 45 contacts per format in partial-coordinate 8R1X and 8 contacts per format in 9CCH restraints that deposit `ZN`/`ZN*` on GLN 48 rather than the coordinate zinc residue.
- The exact 106-contact reviewed coordinate set is pinned by digest `60ed1b785e328fb2d972de6f5d691765012eb622afa489b7a47f41909c8eb141`; any addition, removal, identity, bound, or provenance change fails CI.
- Any implication failure, unresolved discrepancy, parser/projection bug, discrepancy-digest change, scientific-metric change, or reviewed-coordinate-set change makes the command exit nonzero.

This validates conversion safety and format behavior against structures refined with the deposited restraints; it is not an independent Boltz prediction-accuracy benchmark.
