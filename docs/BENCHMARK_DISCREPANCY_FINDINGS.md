# Paired-format discrepancy findings

The July 2026 corpus run audited all 4,177 heavy-contact discrepancies rather
than treating unequal contact counts as parser failures. The complete evidence
is in `benchmark/output/FORMAT_DISCREPANCY_AUDIT.tsv`; each row includes both
formats' source restraint and row IDs, atom expressions, canonical expansions,
physical proton sets, pseudoatom outcome, `N`, averaging policy, source upper
bound, projected heavy pair, projected terms, and final bound.

Expected differences are accepted only by tested predicates: 728 wildcard-set
versus explicit-OR rows, 1,866 x/y-assignment versus compatible physical-set
rows, 1,444 rejected geometric-pseudoatom rows, and 96 rows requiring both the
x/y and pseudoatom predicates. The verified-alias predicate is regression-tested
but is not needed by this corpus. Arbitrary representation and rejection
differences are not allowlisted.

| Entry | Audited rows | Classification | Scientific finding |
|---|---:|---|---|
| 9PQH | 45 | 43 deposition inconsistencies; 2 expected differences | The NEF restraint numbering is shifted relative to its NEF sequence. Forty-one contacts encounter an explicit residue-name conflict and two more resolve the same restraint ID to different target positions. The two common-pair bound differences are explicit expanded OR rows in NEF versus wildcard atom sets in STAR. |
| 9SGX | 813 | expected format differences | NEF uses resolvable `%` and x/y atom/assignment forms. STAR uses geometric Q/M pseudoatoms for the counterpart restraints; these remain rejected by default because no program-specific pseudoatom correction is deposited. |
| 9VUY | 734 | expected format differences | The 631 NEF-only, 7 STAR-only, and 96 different-bound rows trace to x/y or explicit atoms in NEF versus geometric pseudoatoms/canonical expansions in STAR. Replacing Q/M pseudoatoms with atom sets would be an unverified approximation. |
| 43JX | 353 | expected format differences | Both formats emit the same 1,716 heavy pairs. Complete topology-verified STAR canonical expansions are now reconstructed before projection: 591 emitted bounds changed, and restraint 468 now matches NEF with `N=3` and 7.428047646 A. The remaining rows are proven wildcard/explicit-OR or x/y-assignment versus physical-set differences, not lost canonical multiplicity. |
| 6M6O | 149 | expected format differences | NEF x/y rows encode assignment alternatives while STAR `#` expressions encode physical atom sets. The resulting 1 NEF-only, 59 STAR-only, and 89 different-bound contacts preserve those distinct logical statements. |

The remaining audited cases contain 160 expected differences in 8R1X, 1,920
expected differences in 9CCH, and 3 expected differences in 9VQ1. The prior 16
9CCH `ZN`/`ZN*` discrepancy rows are no longer executable discrepancies: all
eight contacts per format are quarantined because GLN B48 has no such atom.
21CC, 9D99, and 9KG4 have exact positive-distance
parity; 8S8O has exact empty-output parity. No reproducible parser/projection
bug remains in the row-level audit, and no discrepancy is unclassified.

The fail-closed gate **passes** because the scientifically reviewed coordinate
limitations are pinned by exact identity and digest rather than waived by
count. Resolution is incomplete for 41 contacts in each format, all involving
atoms absent from the partial 8R1X coordinate model. Four additional malformed
8R1X contacts per format and eight 9CCH GLN/ZN contacts per format are omitted
from executable YAML by component topology, not by coordinate presence. They
remain deterministic provenance-bearing rejections and are never remapped.
The reviewed 82-contact set has SHA-256 digest
`ce8bc2710b5068d3d1cec33794f77bb4e806019a2116a9f4127102b25cfb63e2`;
any identity, bound, provenance, addition, or removal fails CI. Projection
implication itself has 0 failures in 379,449 tested
satisfied-antecedent model cases.
