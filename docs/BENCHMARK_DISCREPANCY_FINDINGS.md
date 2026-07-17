# Paired-format discrepancy findings

The July 2026 corpus run audited all 4,195 heavy-contact discrepancies rather
than treating unequal contact counts as parser failures. The complete evidence
is in `benchmark/output/FORMAT_DISCREPANCY_AUDIT.tsv`; each row includes both
formats' source restraint and row IDs, atom expressions, canonical expansions,
physical proton sets, pseudoatom outcome, `N`, averaging policy, source upper
bound, projected heavy pair, projected terms, and final bound.

Expected differences are accepted only by tested predicates: 987 wildcard-set
versus explicit-OR rows, 1,609 x/y-assignment versus compatible physical-set
rows, 1,444 rejected geometric-pseudoatom rows, and 96 rows requiring both the
x/y and pseudoatom predicates. The verified-alias predicate is regression-tested
but is not needed by this corpus. Arbitrary representation and rejection
differences are not allowlisted.

| Entry | Audited rows | Classification | Scientific finding |
|---|---:|---|---|
| 9PQH | 45 | 43 deposition inconsistencies; 2 expected differences | The NEF restraint numbering is shifted relative to its NEF sequence. Forty-one contacts encounter an explicit residue-name conflict and two more resolve the same restraint ID to different target positions. The two common-pair bound differences are explicit expanded OR rows in NEF versus wildcard atom sets in STAR. |
| 9SGX | 813 | expected format differences | NEF uses resolvable `%` and x/y atom/assignment forms. STAR uses geometric Q/M pseudoatoms for the counterpart restraints; these remain rejected by default because no program-specific pseudoatom correction is deposited. |
| 9VUY | 734 | expected format differences | The 631 NEF-only, 7 STAR-only, and 96 different-bound rows trace to x/y or explicit atoms in NEF versus geometric pseudoatoms/canonical expansions in STAR. Replacing Q/M pseudoatoms with atom sets would be an unverified approximation. |
| 43JX | 355 | expected format differences | Both formats emit the same 1,716 heavy pairs. NEF wildcard atom sets use `N^(1/6)` under `sum-r6`; STAR deposits the canonical protons as explicit OR members, so `N=1` for each member. The different heavy bounds are mathematically expected. |
| 6M6O | 149 | expected format differences | NEF x/y rows encode assignment alternatives while STAR `#` expressions encode physical atom sets. The resulting 1 NEF-only, 59 STAR-only, and 89 different-bound contacts preserve those distinct logical statements. |

The remaining audited cases contain 160 expected differences in 8R1X, 1,920
expected differences plus 16 `ZN`/`ZN*` deposition inconsistencies in 9CCH,
and 3 expected differences in 9VQ1. 21CC, 9D99, and 9KG4 have exact positive-distance
parity; 8S8O has exact empty-output parity. No reproducible parser/projection
bug remains in the row-level audit, and no discrepancy is unclassified.

The fail-closed gate **passes** because the scientifically reviewed coordinate
limitations are pinned by exact identity and digest rather than waived by
count. Resolution is incomplete for 53 contacts in each format: 45 contacts address
residues/atoms absent from the partial 8R1X coordinate model, and 8 9CCH source
restraints place `ZN`/`ZN*` on GLN 48 while the PDB zinc is residue 101. These
are corpus/deposition limitations, not safe targets for an automatic remap.
The reviewed 106-contact set has SHA-256 digest
`60ed1b785e328fb2d972de6f5d691765012eb622afa489b7a47f41909c8eb141`;
any identity, bound, provenance, addition, or removal fails CI. Projection
implication itself has 0 failures in 390,930 tested
satisfied-antecedent model cases.
