# nmr2boltz real-data validation: PDB 6M6O

Date: 2026-07-16

## Outcome

The converter passed the real-data mathematical audit. Across 25,350 model/contact cases for which every independent source restraint was satisfied, there were zero cases where the corresponding projected heavy-atom bound failed. PDB and mmCIF coordinates agreed exactly for all 30,820 atom instances in the 10-model ensemble.

The final NOE-only run read 3,060 author-level restraint groups, retained 3,016 safely projectable groups, emitted 2,707 unique heavy-atom contacts, quarantined seven multi-pair OR groups, and recorded 37 same-parent rejections. The all-distance run read 3,271 groups and emitted 2,877 contacts.

## Deposited inputs

6M6O is a 185-residue, one-chain solution-NMR ensemble with 10 submitted conformers. The coordinates were downloaded in PDB and mmCIF formats from RCSB. The complete deposited NMR constraint file and legacy XPLOR/CNS restraint archive were also downloaded from RCSB. The linked BMRB entry is 28061.

The BMRB entry file contains chemical shifts but no general-distance restraint loop. Therefore, 6M6O_nmr-data.str is the authoritative machine-readable constraint input for this test; 6M6O.mr is the original-program cross-check. The deposited constraint file contains 3,271 distance groups: 3,060 marked NOE, 211 hydrogen-bond/general-distance groups, and 314 torsion constraints outside nmr2boltz's distance-only scope.

All input SHA-256 hashes are in workspace/input/SHA256SUMS. The BMRB-provided MD5 for bmr28061_3.str also matched.

Sources:

- RCSB structure: https://www.rcsb.org/structure/6M6O
- RCSB experiment: https://www.rcsb.org/experimental/6M6O
- RCSB coordinates and restraints: https://files.rcsb.org/download/6M6O_nmr-data.str
- BMRB entry directory: https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr28061/

## Coordinate and distance results

| Check | Result |
|---|---:|
| Coordinate models | 10 PDB and 10 mmCIF |
| Common atom instances | 30,820 |
| PDB-only / mmCIF-only atoms | 0 / 0 |
| Maximum PDB-mmCIF coordinate delta | 0.000000 A |
| Emitted heavy-atom constraints | 2,707 |
| Heavy constraint-model cases | 27,070 |
| Resolved heavy cases | 27,070 |
| Satisfied heavy cases | 27,051 (99.9298%) |
| Constraints satisfied in all 10 models | 2,705 of 2,707 (99.9261%) |
| Bound-versus-distance Pearson correlation | 0.75745 |
| Median observed distance minus bound | -3.1044 A |
| 95th percentile observed distance minus bound | -1.5490 A |
| Source restraint-model cases | 30,600 |
| Satisfied source cases | 28,814 (94.1634%) |
| Source groups satisfied in all models | 2,655 of 3,060 (86.7647%) |
| Indeterminate source cases | 0 |
| Projection implication failures | 0 of 25,350 applicable cases |

Constraint classes were 760 intraresidue, 1,157 sequential, 347 medium-range, and 443 long-range contacts. Thus the result is not explained only by local geometry.

## The 19 heavy-contact violations

All 19 violations came from two intraresidue ARG restraints in deposited list XPLOR-NIH/CNS_distance_restraints_2. They are also present verbatim in the original c-FLIPs_hbonds.tbl section of 6M6O.mr, so this is not a parser invention.

| Source group | Projected atoms | Bound | Deposited ensemble distance | Models violated |
|---|---|---:|---:|---:|
| list 2, restraint 1 | A:3:N - A:3:NH2 | 5.08 A | 4.315-6.247 A | 9 of 10 |
| list 2, restraint 2 | A:65:N - A:65:NH2 | 5.08 A | 6.460-6.645 A | 10 of 10 |

Their original HH21-N source distances also violate the deposited 4.0 A upper bounds in all 10 models. This is direct evidence of a source-restraint/deposited-ensemble inconsistency. It does not contradict the converter: the mathematical audit conditions only on source restraints that the coordinate model actually satisfies, and that audit had zero failures.

All other 2,705 emitted heavy contacts were satisfied in every conformer.

## Refinements prompted by the real test

The deposited NMR-STAR file uses the XPLOR/CNS digit wildcard # in author atom names. Topology resolution already normalized # to the NEF wildcard %, but canonical expansion-row deduplication did not recognize #. This inflated the NOE report to 7,701 apparent source alternatives.

The parser now treats # as a set expression during deduplication. A regression test was added. After rebuilding, the same physical data produce the correct 3,060 author-level alternatives while the number of safe contacts, ambiguous groups, and rejections remains unchanged. This is a provenance and hypothesis-weighting correctness fix, not a relaxation of scientific safety.

The ensemble comparator also handles protonation-specific wildcard membership explicitly. For example, HIS 92 contains HD2 but not HD1 in the deposited tautomer. The source wildcard is evaluated over atoms actually present in each hydrogen-complete conformer, while the absence is recorded in the CSV audit.

The Dockerfile now declares UID/GID 65532 as its default runtime user. A no-user-override smoke test confirmed uid=65532 and gid=65532.

The Boltz target sequence is now machine-validated instead of relying on the
manual sequence-map checklist. The checksum-pinned `workspace/benchmark.yaml`
case verified all 185 mapped residues and all 2,707 emitted constraints against
`workspace/input/6M6O_boltz.yaml`, with zero target-validation errors or warnings.

## Robustness tests

The final image passed all 39 Pytest tests and a containerized deterministic stress run containing:

- 100,000 randomized unnormalized sum-r6 implication cases;
- 100,000 constructive triangle-inequality projection cases;
- 25,000 outward-rounding cases;
- 10,000 randomized OR-max/AND-min merge-order cases;
- all 850 built-in hydrogen-to-parent topology mappings;
- NEF, NMR-STAR, compressed input, custom components, all averaging policies, and 32 reproducible ambiguity hypotheses.

## Docker result

Only one nmr2boltz image remains:

- image: nmr2boltz:0.1.0-validated
- digest: sha256:bb97a627d721afb4985caf3c5b86a48486b70f0cd1436b0129149ad172a22405
- reported size: 243 MB
- default user: 65532:65532

Every conversion, comparison, smoke, and stress container run mounted the repository workspace at /workspace. Runs used no network and a read-only root filesystem.

## Reproduction and audit files

- Exact PowerShell commands: workspace/output/RUN_COMMANDS.ps1
- Converter provenance: workspace/output/nmr2boltz_noe/conversion_report.json
- Per-model heavy distances: workspace/output/coordinate_comparison/heavy_atom_distances.csv
- Per-contact heavy summary: workspace/output/coordinate_comparison/heavy_atom_summary.csv
- Per-model source evaluation: workspace/output/coordinate_comparison/source_restraint_distances.csv
- Per-group source summary: workspace/output/coordinate_comparison/source_restraint_summary.csv
- Machine-readable result: workspace/output/coordinate_comparison/summary.json
- Stress result: workspace/output/STRESS_VALIDATION.json
- Target-aware benchmark result: workspace/output/benchmark/benchmark_summary.json

## Scientific interpretation and limit

This benchmark strongly supports parser correctness, atom-name generality, sequence mapping, conservative proton-parent projection, ambiguity quarantine, provenance, serialization, and container reproducibility for a challenging real 185-residue NMR entry.

It does not prove that restrained Boltz produces a more accurate fold than unrestrained Boltz. No GPU Boltz diffusion campaign was run here. Predictive efficacy still requires restrained/unrestrained replicate prediction, blind hold-out restraints, and independent structural-quality controls.
