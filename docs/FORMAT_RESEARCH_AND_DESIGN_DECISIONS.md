# Format research and design decisions

## 1. Data sources reviewed

The design was based on the primary format and implementation sources rather than ad-hoc examples:

- the official NEF specification repository, including its commented example and atom-name semantics;
- BMRB NMR-STAR dictionaries and the BMRB-maintained PyNMRSTAR parser/API behavior;
- wwPDB mappings between NEF distance-restraint categories and NMR-STAR general distance constraints;
- real deposited NMR-STAR restraint rows containing repeated constraint IDs, canonical atom expansion, author wildcard names, and OR-like membership;
- the wwPDB Chemical Component Dictionary atom/bond model;
- the Boltz/BoltzUI atom-contact parser and potential representation, including existing union-index behavior.

The parser intentionally selects loops by **category and tag**, not saveframe name. Saveframe names are user-controlled and are not a stable schema key.

## 2. NEF-to-NMR-STAR field correspondence used by the converter

| Concept | NEF | NMR-STAR |
|---|---|---|
| List/saveframe category | `nef_distance_restraint_list` | `general_distance_constraints` |
| Row loop category | `_nef_distance_restraint` | `_Gen_dist_constraint` |
| Restraint ID | `restraint_id` | `ID` |
| Combination ID | `restraint_combination_id` | `Combination_ID` |
| Chain | `chain_code_1/2` | `Auth_asym_ID_1/2`, fallback `Entity_assembly_ID_1/2` |
| Sequence | `sequence_code_1/2` | `Auth_seq_ID_1/2`, fallback `Comp_index_ID_1/2` |
| Residue | `residue_name_1/2` | `Auth_comp_ID_1/2`, fallback `Comp_ID_1/2` |
| Atom expression | `atom_name_1/2` | `Auth_atom_name_1/2`, then `Auth_atom_ID_1/2`, with `Atom_ID_1/2` as canonical hint |
| Upper bound | `upper_limit` | `Distance_upper_bound_val` |
| Lower bound | `lower_limit` | `Distance_lower_bound_val` |
| Target | `target_value` | `Target_val` |
| Target uncertainty | `target_value_uncertainty` | `Target_val_uncertainty` |
| Weight | `weight` | `Weight` |

The exact NMR-STAR dictionary can evolve. Tag lookup is case-normalized at the leaf-name level, and optional fields are tolerated.

## 3. Real-world patterns the code is designed to survive

### 3.1 Canonical expansion of one author atom set

A source expression such as `HG1%` may be deposited in NMR-STAR as three rows whose canonical atom IDs are `HG11`, `HG12`, and `HG13`, while the author atom-name field remains `HG1%`. Treating the three rows as three independent restraints would be wrong. The code reconstructs the author-level expression and deduplicates the canonical expansion rows.

### 3.2 Aromatic and symmetric ambiguity

Aromatic expressions such as `HD%` can map to different heavy parents. The code does not collapse them merely because they are chemically similar or symmetry-related. A Boltz structure has explicit atom names, so `CD1` and `CD2` remain distinct alternatives unless the exact component/assignment model proves otherwise.

### 3.3 Insertion-like sequence codes

NEF sequence codes can be values such as `24B`. Parsing them as integers loses information and can silently target the wrong residue. The source code remains a string; the one-based Boltz position is derived from the molecular-system sequence loop or a user map.

### 3.4 Missing bounds

Some lists contain target values but no explicit upper limit. The converter does not invent a bound by default. An expert can request a documented heuristic policy.

### 3.5 Modified components

Author atom names for nonstandard residues cannot be handled safely by a 20-amino-acid lookup table. Embedded chemical-component loops or a local CCD file are used to find direct bonds. Unknown hydrogens are rejected; unknown heavy atoms can be passed through with a warning because no proton-parent inference is needed.

### 3.6 Complex logical combinations

A non-null combination identifier may encode logic across several restraint IDs. Without a complete expression parser, flattening it as OR or AND is unsafe. Such groups are preserved in the report and rejected from ordinary YAML.

## 4. Approaches considered

## 4.1 Subtracting or adding one fixed constant

**Proposal:** convert every H-H bound to parent-parent by adding a universal value such as 2.0 Å.

**Advantages:** trivial and fast.

**Problems:** C-H, N-H, O-H, S-H, and Se-H bonds differ; H-X needs only one term; heavy-heavy needs none; pseudoatoms and atom sets remain unresolved. A single constant can be either unnecessarily weak or non-conservative.

**Decision:** rejected in favor of topology- and element-aware direct-bond envelopes.

## 4.2 Using ideal local geometry to calculate a tighter parent bound

**Proposal:** assume tetrahedral/planar geometry, known bond angles, and a rotamer to derive a tighter parent-parent distance.

**Advantages:** stronger guidance.

**Problems:** the result depends on side-chain conformation, stereochemistry, protonation, and which proton in a set produced the NOE. It can exclude a valid structure.

**Decision:** not used in the conservative converter. It could be offered later as an explicitly model-dependent mode.

## 4.3 Reconstructing hydrogens before/during Boltz inference

**Proposal:** create virtual hydrogen coordinates from the current heavy-atom state and evaluate the original restraint directly.

**Advantages:** preserves more original geometry and can represent proton-specific information.

**Problems:** it requires differentiable hydrogen placement, protonation/tautomer handling, rotamer-dependent choices, and atom-set/ensemble logic inside the model. It is substantially more invasive than the current Boltz potential API.

**Decision:** excellent long-term direction, but not the first robust integration. Use hydrogen reconstruction for post-prediction validation now.

## 4.4 Emitting all ambiguous parent pairs

**Proposal:** expand an ambiguous restraint and add every pair as an `atom_contact`.

**Advantages:** no changes to Boltz internals.

**Problems:** mathematically changes OR to AND and can force mutually incompatible contacts.

**Decision:** explicitly prohibited by the default converter.

## 4.5 Choosing the first alternative

**Proposal:** select the first deposited member.

**Advantages:** deterministic and compatible.

**Problems:** row order does not imply probability or correctness.

**Decision:** not a default. Optional hypothesis files sample/cover alternatives and label the assumption.

## 4.6 Union-aware potential

**Proposal:** share a union index across alternatives so the potential acts as a soft minimum.

**Advantages:** best match to source ambiguity; Boltz already has relevant potential machinery.

**Problems:** the parser/schema and token-conditioning path need changes; union sharpness requires calibration.

**Decision:** recommended BoltzUI extension. A proposed schema is emitted separately.

## 4.7 Full Bayesian restraint model

**Proposal:** jointly model assignment probability, calibration uncertainty, dynamics, and structural coordinates.

**Advantages:** scientifically expressive.

**Problems:** requires reliable priors, a training/evaluation program, and a much larger implementation.

**Decision:** future research direction after a transparent deterministic baseline.

## 5. Why the selected algorithm is robust

The selected method is robust in the narrow but important sense that every emitted ordinary contact follows from explicit assumptions recorded in the audit:

1. source identifiers map to a verified modeled residue;
2. each omitted proton has a unique directly bonded modeled parent;
3. the source upper bound is explicit or its derivation is labeled;
4. atom-set multiplicity uses a user-visible averaging policy;
5. every source alternative maps to the same heavy pair;
6. Boltz range adaptation never tightens the bound;
7. all discarded or ambiguous information remains available for review.

It is also deliberately modular. Sequence mapping, topology, source parsing, projection, and output can be tested independently.

## 6. BMRB/wwPDB retrieval considerations

For automated pipelines, prefer:

- downloading the exact deposited NEF/NMR-STAR file and recording its checksum;
- recording the BMRB/PDB entry identifier and retrieval date;
- using a fixed API version for production rather than an unversioned current endpoint;
- caching CCD components locally for reproducibility;
- rate-limiting BMRB API retrieval and avoiding repeated full-entry downloads;
- retaining the original file alongside `conversion_report.json`.

The converter itself does not silently download entries or CCD components. This avoids non-reproducible network behavior, credential/rate-limit surprises, and later dictionary drift. A separate acquisition layer can be added around it.

## 7. Recommended dataset-level protocol

1. Identify one or more deposited restraint lists and their original structure-calculation software.
2. Ask whether wildcard/pseudoatom corrections were applied before deposition.
3. Convert with explicit bounds only and `sum-r6` as the conservative exploratory baseline.
4. Review how many restraints become over-20 Å or ambiguous.
5. Repeat with the scientifically confirmed averaging convention.
6. Validate the source-to-Boltz sequence map.
7. Select a high-confidence long-range subset if the full set is redundant.
8. Run unrestrained, safe-restraint, and ambiguity-aware/hypothesis ensembles.
9. Score original restraints after hydrogen reconstruction.
10. Report sensitivity to policy choices rather than one best-looking model.

## 8. Paired-format benchmark and remaining coverage

The current 12-structure benchmark spans:

- proteins, RNA, DNA, and a protein-DNA complex;
- monomers, homomultimers, heterocomplexes, and interchain restraints;
- wildcard atom sets, x/y assignments, geometric pseudoatoms, and complex logic;
- non-one-based numbering, chain aliases, modified residues, caps, and ions;
- partial coordinate observations and a sequence-only/no-distance negative case;
- paired NEF and NMR-STAR representations for differential comparison.

The benchmark showed that identical restraint-group counts do not imply
equivalent projected contacts or bounds. Pseudoatom representation and atom-set
multiplicity remain the dominant parity limitations. These differences must be
reported at the physical atom-pair and bound level, not inferred from aggregate
counts.

Future additions should include intentionally corrupted inputs, independent
author-generated files, protein-RNA complexes, sparse/noisy lists with known
incorrect assignments, and systems selected for a blinded Boltz prediction
comparison.

Prediction-stage metrics should include fold recovery, held-out original-restraint satisfaction, stereochemical quality, diversity, robustness to wrong restraints, and the difference between OR-aware and incorrect-AND treatment.
