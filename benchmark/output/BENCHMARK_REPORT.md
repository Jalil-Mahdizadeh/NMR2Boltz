# nmr2boltz paired-format benchmark

All 12 deposited ensembles were converted from both NEF and NMR-STAR with conservative defaults, then audited against every PDB conformer using sequence-aware coordinate alignment.

| Case | NEF exact | NEF unions | STAR exact | STAR unions | NEF PDB satisfaction | STAR PDB satisfaction | Pair/bound parity | Implication failures |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 21CC | 176 | 16 | 176 | 16 | 99.72% | 99.72% | yes | 0 |
| 43JX | 1716 | 16 | 1716 | 16 | 99.96% | 99.96% | no | 0 |
| 6M6O | 2819 | 71 | 2877 | 7 | 99.72% | 99.60% | no | 0 |
| 8R1X | 1114 | 165 | 1114 | 165 | 99.74% | 99.74% | no | 0 |
| 8S8O | 0 | 0 | 0 | 0 | N/A | N/A | yes | 0 |
| 9CCH | 3350 | 292 | 3516 | 75 | 99.98% | 100.00% | no | 0 |
| 9D99 | 456 | 0 | 456 | 0 | 100.00% | 100.00% | yes | 0 |
| 9KG4 | 406 | 0 | 406 | 0 | 99.74% | 99.74% | yes | 0 |
| 9PQH | 3 | 0 | 44 | 15 | 3.33% | 100.00% | no | 0 |
| 9SGX | 1578 | 766 | 765 | 0 | 100.00% | 100.00% | no | 0 |
| 9VQ1 | 53 | 3 | 56 | 0 | 95.00% | 94.20% | no | 0 |
| 9VUY | 1327 | 727 | 703 | 0 | 99.96% | 100.00% | no | 0 |

## Current result

- 24/24 conversions completed; 2 are valid empty distance conversions for 8S8O.
- NEF: 12998 exact contacts and 2056 union groups; 99.88% resolved PDB satisfaction for exact contacts.
- NMR-STAR: 11829 exact contacts and 294 union groups; 99.86% resolved PDB satisfaction for exact contacts.
- Conservative implication failures: 0 across 379449 satisfied-antecedent cases.
- Final executable-topology violations: 0; every emitted endpoint is proven by its mapped component dictionary.
- Exact NEF/STAR pair-and-bound parity: 3/11 positive-distance cases; 8S8O also has exact empty-output parity.

## Constraint output split

- `atom_constraints_exact.yaml` contains only non-ambiguous `atom_contact` constraints counted in the `exact` columns.
- `atom_constraints_union.yaml` contains only ambiguous `atom_contact_union` OR groups counted in the `unions` columns; every alternative retains its own conservatively rounded bound.
- Both files are written for every conversion, including explicit empty `constraints: []` files. Exact contacts never appear in the union file, and union groups never appear in the exact file.
- PDB satisfaction and implication denominators cover exact contacts only. Union alternatives are not treated as simultaneous contacts because that would convert OR semantics into AND.

## Row-level format discrepancy audit

- 4177 NEF-only, STAR-only, or different-bound contacts were audited to source rows and physical proton sets.
- Classifications: 4134 scientifically expected format differences; 43 deposition inconsistencies; 0 unresolved; 0 parser/projection bugs.
- Reviewed audit digest: `e710f092a339fcf8d5d7a57d35207cc4b05bea25fe92c0c30b1f97815605c677`.
- `expected_format_difference` is allowlisted only for wildcard-set versus explicit OR, x/y assignment versus a compatible physical set, rejected Q/M pseudoatoms, or verified canonical aliases.
- Complete topology-verified NMR-STAR canonical OR expansions are reconstructed before projection, so their atom-set multiplicity is applied once and their source rows remain auditable.
- 9PQH contains 43 contact discrepancies caused by its NEF sequence/residue numbering conflict; its remaining two different-bound rows are the expected explicit-OR versus wildcard atom-set distinction.
- 9CCH's 16 cross-format `ZN`/`ZN*` discrepancy rows disappeared because all eight invalid contacts per format are now quarantined before emission; its remaining format differences satisfy the semantic allowlist.
- 43JX, 6M6O, 9SGX, and 9VUY differences are justified by explicit-OR, wildcard atom-set, x/y assignment, or rejected geometric-pseudoatom semantics.

## Fail-closed gate

- Gate status: **PASS** (0 failure categories).
- Coordinate resolution gaps: 41 NEF and 41 STAR contacts. These are not omitted from the denominator silently.
- The current gaps are 41 contacts per format in the partial-coordinate 8R1X ensemble. Four malformed 8R1X contacts and all eight GLN/ZN 9CCH contacts per format are quarantined before coordinate evaluation.
- The exact 82-contact reviewed coordinate set is pinned by digest `ce8bc2710b5068d3d1cec33794f77bb4e806019a2116a9f4127102b25cfb63e2`; any addition, removal, identity, bound, or provenance change fails CI.
- Any emitted atom-topology violation, implication failure, unresolved discrepancy, parser/projection bug, discrepancy-digest change, scientific-metric change, or reviewed-coordinate-set change makes the command exit nonzero.

This validates conversion safety and format behavior against structures refined with the deposited restraints; it is not an independent Boltz prediction-accuracy benchmark.
