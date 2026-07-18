# Graph Report - nmr2boltz_toolkit  (2026-07-18)

## Corpus Check
- 126 files · ~17,120,885 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 569 nodes · 1502 edges · 29 communities (25 shown, 4 thin omitted)
- Extraction: 93% EXTRACTED · 7% INFERRED · 0% AMBIGUOUS · INFERRED: 101 edges (avg confidence: 0.51)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `67f094d7`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 13|Community 13]]
- [[_COMMUNITY_Community 14|Community 14]]
- [[_COMMUNITY_Community 15|Community 15]]
- [[_COMMUNITY_Community 16|Community 16]]
- [[_COMMUNITY_Community 17|Community 17]]
- [[_COMMUNITY_Community 18|Community 18]]
- [[_COMMUNITY_Community 19|Community 19]]
- [[_COMMUNITY_Community 20|Community 20]]
- [[_COMMUNITY_Community 22|Community 22]]
- [[_COMMUNITY_Community 23|Community 23]]
- [[_COMMUNITY_Community 24|Community 24]]
- [[_COMMUNITY_Community 25|Community 25]]
- [[_COMMUNITY_Community 26|Community 26]]
- [[_COMMUNITY_Community 27|Community 27]]
- [[_COMMUNITY_Community 29|Community 29]]
- [[_COMMUNITY_Community 30|Community 30]]

## God Nodes (most connected - your core abstractions)
1. `TopologyLibrary` - 58 edges
2. `ProjectionSettings` - 48 edges
3. `Scientific and computational method for projecting proton NMR distance restraints onto Boltz heavy-atom contacts` - 40 edges
4. `project_document()` - 37 edges
5. `parse_star_document()` - 35 edges
6. `BoltzAtom` - 32 edges
7. `ConversionReport` - 31 edges
8. `write_outputs()` - 30 edges
9. `SequenceRecord` - 24 edges
10. `ComponentTopology` - 23 edges

## Surprising Connections (you probably didn't know these)
- `test_yaml_rounding_is_outward_for_upper_bounds()` --calls--> `_rounded()`  [EXTRACTED]
  tests/test_robustness.py → src/nmr2boltz/output.py
- `_Saveframe` --uses--> `StarDataError`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/star.py
- `_StarLoop` --uses--> `StarDataError`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/star.py
- `_Saveframe` --uses--> `TopologyResolutionError`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/topology.py
- `_StarLoop` --uses--> `TopologyResolutionError`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/topology.py

## Import Cycles
- None detected.

## Communities (29 total, 4 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.09
Nodes (52): list, Rejection, _averaging_factor(), _canonical_pair(), _merge_independent_constraints(), _merge_or_alternatives(), project_document(), _project_group() (+44 more)

### Community 1 - "Community 1"
Cohesion: 0.14
Nodes (29): Atom, AtomKey, Point, test_structure_alignment_accepts_equivalent_nucleic_acid_names(), test_structure_alignment_rekeys_author_numbers_to_boltz_indices(), _coordinate_summary(), _align_components(), align_structure_to_sequence_map() (+21 more)

### Community 2 - "Community 2"
Cohesion: 0.05
Nodes (36): 10.1 Safe default, 10.2 Union-aware potential, 10.3 Assignment hypotheses, 10. Ambiguous groups and recommended execution strategies, 12. Suggested confidence tiers, 14. Reproducibility and auditability, 15. Questions for NMR expert review before production deployment, 16. Primary specifications and references (+28 more)

### Community 3 - "Community 3"
Cohesion: 0.10
Nodes (45): Decimal, float, Conservative NMR distance-restraint projection for Boltz atom contacts., AmbiguousGroup, BoltzAtom, ConversionReport, EmittedConstraint, ProjectedAlternative (+37 more)

### Community 4 - "Community 4"
Cohesion: 0.32
Nodes (13): Random, _add(), _distance(), main(), _projected(), Deterministic mathematical and converter stress validation for nmr2boltz.  This, _unit_vector(), validate_boolean_merge_order() (+5 more)

### Community 5 - "Community 5"
Cohesion: 0.11
Nodes (25): AtomChoice, AtomSetChoice, atom_topology_violations(), AtomTopologyValidationError, AtomTopologyViolation, build_builtin_topologies(), _compile_atom_pattern(), component_topology_snapshot() (+17 more)

### Community 6 - "Community 6"
Cohesion: 0.15
Nodes (37): as_float(), clean(), Endpoint, Return a stripped STAR value or ``None`` for STAR null tokens., RawAlternative, RestraintGroup, SequenceRecord, _category() (+29 more)

### Community 7 - "Community 7"
Cohesion: 0.08
Nodes (46): ArgumentParser, Namespace, BenchmarkCaseResult, BenchmarkManifestError, BenchmarkRunResult, _bool_option(), _case_metrics(), _expected_metrics() (+38 more)

### Community 8 - "Community 8"
Cohesion: 0.09
Nodes (22): 1. Data sources reviewed, 2. NEF-to-NMR-STAR field correspondence used by the converter, 3.1 Canonical expansion of one author atom set, 3.2 Aromatic and symmetric ambiguity, 3.3 Insertion-like sequence codes, 3.4 Missing bounds, 3.5 Modified components, 3.6 Complex logical combinations (+14 more)

### Community 9 - "Community 9"
Cohesion: 0.11
Nodes (17): Batch hypothesis fallback, Critical token-conditioning rule, Design proposal: ambiguity-aware `atom_contact_union` constraints in BoltzUI/Boltz-2, End-to-end tests, Energy tests, Feature tests, Multiple thresholds and duplicate parents, Optional weights (+9 more)

### Community 11 - "Community 11"
Cohesion: 0.10
Nodes (17): Paired-format discrepancy findings, Ambiguity rules, BoltzUI integration, Docker, Documentation, Important limitations, Installation, nmr2boltz (+9 more)

### Community 12 - "Community 12"
Cohesion: 0.18
Nodes (10): Coordinate and distance results, Deposited inputs, Docker result, nmr2boltz real-data validation: PDB 6M6O, Outcome, Refinements prompted by the real test, Reproduction and audit files, Robustness tests (+2 more)

### Community 13 - "Community 13"
Cohesion: 0.22
Nodes (8): Atom semantics, Boltz usage, Identifier mapping, NMR expert review checklist, Post-prediction validation, Projection, Sign-off, Source and calibration

### Community 14 - "Community 14"
Cohesion: 0.39
Nodes (4): inline_markup(), main(), parse_markdown(), table_from_markdown()

### Community 15 - "Community 15"
Cohesion: 0.25
Nodes (7): generated, hypotheses, requested, seed, skipped_groups_without_boltz_compatible_alternatives, total_possible_combinations, warning

### Community 16 - "Community 16"
Cohesion: 0.25
Nodes (7): generated, hypotheses, requested, seed, skipped_groups_without_boltz_compatible_alternatives, total_possible_combinations, warning

### Community 17 - "Community 17"
Cohesion: 0.33
Nodes (5): Ambiguous restraints, BoltzUI integration guide, Exact-contact schema, Recommended run protocol, Residue numbering

### Community 18 - "Community 18"
Cohesion: 0.29
Nodes (7): 11.1 Parsing validation, 11.2 Identifier validation, 11.3 Geometric validation of candidate structures, 11.4 Structural and statistical controls, 11.5 Executed validation record, 11.5 Paired-format discrepancy audit and CI gate, 11. Validation protocol

### Community 19 - "Community 19"
Cohesion: 0.33
Nodes (6): 13.1 Probabilistic union restraints, 13.2 Ensemble-aware objective, 13.3 Iterative assignment and folding, 13.4 Information-aware restraint selection, 13.5 Broader experimental data, 13. Potential methodological extensions

### Community 20 - "Community 20"
Cohesion: 0.33
Nodes (5): Included, nmr2boltz toolkit release notes, Safety defaults, Unreleased, Validation

### Community 22 - "Community 22"
Cohesion: 0.50
Nodes (4): 1.1 Included data, 1.2 Excluded or deliberately non-flattened information, 1.3 Interpretation of a Boltz contact, 1. Scope and intended use

### Community 23 - "Community 23"
Cohesion: 0.50
Nodes (4): 9.1 Lower than 2 A, 9.2 Greater than 20 A, 9.3 Precision, 9. Boltz compatibility and bound handling

### Community 24 - "Community 24"
Cohesion: 0.40
Nodes (4): Current result, Fail-closed gate, nmr2boltz paired-format benchmark, Row-level format discrepancy audit

### Community 30 - "Community 30"
Cohesion: 0.07
Nodes (71): Any, _parse_entry_file(), Path, _alternative(), _endpoint(), _gate_run(), _report(), _semantic() (+63 more)

## Knowledge Gaps
- **147 isolated node(s):** `nmr2boltz`, `warning`, `seed`, `requested`, `generated` (+142 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **4 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `TopologyLibrary` connect `Community 0` to `Community 1`, `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 7`, `Community 30`?**
  _High betweenness centrality (0.053) - this node is a cross-community bridge._
- **Why does `BoltzAtom` connect `Community 3` to `Community 0`, `Community 4`, `Community 5`, `Community 6`, `Community 7`?**
  _High betweenness centrality (0.025) - this node is a cross-community bridge._
- **Why does `project_document()` connect `Community 0` to `Community 3`, `Community 5`, `Community 6`, `Community 7`, `Community 30`?**
  _High betweenness centrality (0.023) - this node is a cross-community bridge._
- **Are the 13 inferred relationships involving `TopologyLibrary` (e.g. with `BenchmarkCaseResult` and `BenchmarkManifestError`) actually correct?**
  _`TopologyLibrary` has 13 INFERRED edges - model-reasoned connections that need verification._
- **Are the 18 inferred relationships involving `ProjectionSettings` (e.g. with `BenchmarkCaseResult` and `BenchmarkManifestError`) actually correct?**
  _`ProjectionSettings` has 18 INFERRED edges - model-reasoned connections that need verification._
- **What connects `nmr2boltz`, `Conservative NMR distance-restraint projection for Boltz atom contacts.`, `Raised when a benchmark manifest does not satisfy the versioned schema.` to the rest of the system?**
  _174 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.08644067796610169 - nodes in this community are weakly interconnected._