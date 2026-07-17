# nmr2boltz paired-format benchmark

All 12 deposited ensembles were converted from both NEF and NMR-STAR with conservative defaults, then audited against every PDB conformer using sequence-aware coordinate alignment.

| Case | NEF contacts | STAR contacts | NEF PDB satisfaction | STAR PDB satisfaction | Exact parity | Implication failures |
|---|---:|---:|---:|---:|---:|---:|
| 21CC | 176 | 176 | 99.72% | 99.72% | yes | 0 |
| 43JX | 1716 | 1716 | 99.96% | 99.92% | no | 0 |
| 6M6O | 2819 | 2877 | 99.72% | 99.60% | no | 0 |
| 8R1X | 1112 | 1118 | 99.74% | 99.74% | no | 0 |
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
- NEF: 13004 contacts, 99.88% resolved PDB satisfaction.
- NMR-STAR: 11841 contacts, 99.86% resolved PDB satisfaction.
- Conservative implication failures: 0 across 390810 satisfied-antecedent cases.
- Exact NEF/STAR pair-and-bound parity: 3/11 positive-distance cases; 8S8O also has exact empty-output parity.

## Comparison with the initial audit

- Successful conversions: 22 -> 24.
- No-distance failures: 2 -> 0.
- FASTA outputs: 0 -> 24.
- Emitted contacts: 13004 -> 13004 NEF and 11841 -> 11841 STAR under unchanged conservative projection policies.
- Sequence/residue conflicts are now explicit: 312 rejection records use `sequence_residue_mismatch`.
- Large pseudoatom and atom-set parity differences remain intentionally visible rather than being silently approximated.

This validates conversion safety and format behavior against structures refined with the deposited restraints; it is not an independent Boltz prediction-accuracy benchmark.
