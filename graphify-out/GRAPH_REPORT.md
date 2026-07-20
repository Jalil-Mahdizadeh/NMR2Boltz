# Graph Report - nmr2boltz_toolkit  (2026-07-20)

## Corpus Check
- 135 files · ~17,614,864 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 691 nodes · 1930 edges · 30 communities (27 shown, 3 thin omitted)
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 120 edges (avg confidence: 0.51)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `75cc9a3c`
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
- [[_COMMUNITY_Community 10|Community 10]]
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
- [[_COMMUNITY_Community 21|Community 21]]
- [[_COMMUNITY_Community 22|Community 22]]
- [[_COMMUNITY_Community 23|Community 23]]
- [[_COMMUNITY_Community 24|Community 24]]
- [[_COMMUNITY_Community 26|Community 26]]
- [[_COMMUNITY_Community 27|Community 27]]
- [[_COMMUNITY_Community 29|Community 29]]
- [[_COMMUNITY_Community 30|Community 30]]

## God Nodes (most connected - your core abstractions)
1. `TopologyLibrary` - 75 edges
2. `ProjectionSettings` - 59 edges
3. `project_document()` - 50 edges
4. `Scientific and computational method for projecting proton NMR distance restraints onto Boltz heavy-atom contacts` - 42 edges
5. `parse_star_document()` - 41 edges
6. `BoltzAtom` - 38 edges
7. `write_outputs()` - 38 edges
8. `ConversionReport` - 35 edges
9. `SequenceRecord` - 26 edges
10. `TopologyResolutionError` - 26 edges

## Surprising Connections (you probably didn't know these)
- `test_all_intrachain_input_writes_empty_constraint_files_and_strict_audit()` --calls--> `main()`  [EXTRACTED]
  tests/test_intrachain_filter.py → src/nmr2boltz/cli.py
- `test_cli_flag_writes_only_interchain_exact_and_union_yaml()` --calls--> `main()`  [EXTRACTED]
  tests/test_intrachain_filter.py → src/nmr2boltz/cli.py
- `test_cli_flag_filters_exact_and_union_outputs()` --calls--> `main()`  [EXTRACTED]
  tests/test_intraresidue_filter.py → src/nmr2boltz/cli.py
- `test_compressed_nef_input()` --calls--> `parse_star_document()`  [EXTRACTED]
  tests/test_logic.py → src/nmr2boltz/star.py
- `_Saveframe` --uses--> `TopologyResolutionError`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/topology.py

## Import Cycles
- None detected.

## Communities (30 total, 3 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.11
Nodes (45): command_convert(), project_document(), ProjectionSettings, parse_star_document(), ParsedStarDocument, StarDataError, TopologyLibrary, test_project_document_quarantines_complete_union_above_configured_maximum() (+37 more)

### Community 1 - "Community 1"
Cohesion: 0.09
Nodes (61): Atom, AtomKey, Namespace, _parse_entry_file(), PairKey, Path, Point, _atom() (+53 more)

### Community 2 - "Community 2"
Cohesion: 0.05
Nodes (38): 10.1 Safe default, 10.2 Union-aware potential, 10.3 Assignment hypotheses, 10.4 Intraresidue exclusion, 10.5 Inter-chain-only selection, 10. Ambiguous groups and recommended execution strategies, 12. Suggested confidence tiers, 14. Reproducibility and auditability (+30 more)

### Community 3 - "Community 3"
Cohesion: 0.07
Nodes (71): Decimal, float, list, AmbiguousGroup, BoltzAtom, EmittedConstraint, ProjectedAlternative, Rejection (+63 more)

### Community 4 - "Community 4"
Cohesion: 0.32
Nodes (13): Random, _add(), _distance(), main(), _projected(), Deterministic mathematical and converter stress validation for nmr2boltz.  This, _unit_vector(), validate_boolean_merge_order() (+5 more)

### Community 5 - "Community 5"
Cohesion: 0.08
Nodes (43): AtomChoice, AtomSetChoice, clean(), Return a stripped STAR value or ``None`` for STAR null tokens., atom_topology_violations(), AtomTopologyValidationError, AtomTopologyViolation, build_builtin_topologies() (+35 more)

### Community 6 - "Community 6"
Cohesion: 0.12
Nodes (45): as_float(), Endpoint, RawAlternative, SequenceRecord, _alternative_row_key(), _canonical_expansion_family_key(), _canonical_expansion_issue(), _canonical_expansion_semantics_key() (+37 more)

### Community 7 - "Community 7"
Cohesion: 0.25
Nodes (7): `benchmark`, Command structure, `convert`, `convert` effective-default notes, Exit statuses, Global flags, NMR2Boltz command-line reference

### Community 8 - "Community 8"
Cohesion: 0.09
Nodes (22): 1. Data sources reviewed, 2. NEF-to-NMR-STAR field correspondence used by the converter, 3.1 Canonical expansion of one author atom set, 3.2 Aromatic and symmetric ambiguity, 3.3 Insertion-like sequence codes, 3.4 Missing bounds, 3.5 Modified components, 3.6 Complex logical combinations (+14 more)

### Community 9 - "Community 9"
Cohesion: 0.11
Nodes (17): Batch hypothesis fallback, Critical token-conditioning rule, Design proposal: ambiguity-aware `atom_contact_union` constraints in BoltzUI/Boltz-2, End-to-end tests, Energy tests, Feature tests, Multiple thresholds and duplicate parents, Optional weights (+9 more)

### Community 10 - "Community 10"
Cohesion: 0.44
Nodes (8): _convert(), _filtered_rejection(), test_cli_flag_filters_exact_and_union_outputs(), test_exclude_intraresidue_leaves_interresidue_union_unaffected(), test_exclude_intraresidue_preserves_exact_contact_behavior_and_provenance(), test_exclude_intraresidue_rejects_complete_all_intraresidue_union(), test_exclude_intraresidue_rejects_complete_mixed_scope_union(), test_writer_fails_closed_if_intraresidue_excluded_report_is_mutated()

### Community 11 - "Community 11"
Cohesion: 0.08
Nodes (20): Paired NEF/NMR-STAR benchmark, Exact-contact distance check, Paired-format discrepancy findings, Ambiguity rules, BoltzUI integration, Docker, Documentation, Important limitations (+12 more)

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
Cohesion: 0.25
Nodes (7): Ambiguous restraints, BoltzUI integration guide, Exact-contact schema, Inter-chain-only constraints, Intraresidue-excluded constraints, Recommended run protocol, Residue numbering

### Community 18 - "Community 18"
Cohesion: 0.25
Nodes (8): 11.1 Parsing validation, 11.2 Identifier validation, 11.3 Geometric validation of candidate structures, 11.4 Structural and statistical controls, 11.5 Paired-format discrepancy audit and CI gate, 11.6 Exact-contact PDB-model distance matrix, 11.7 Executed validation record, 11. Validation protocol

### Community 19 - "Community 19"
Cohesion: 0.33
Nodes (6): 13.1 Probabilistic union restraints, 13.2 Ensemble-aware objective, 13.3 Iterative assignment and folding, 13.4 Information-aware restraint selection, 13.5 Broader experimental data, 13. Potential methodological extensions

### Community 20 - "Community 20"
Cohesion: 0.33
Nodes (5): Included, nmr2boltz toolkit release notes, Safety defaults, Unreleased, Validation

### Community 21 - "Community 21"
Cohesion: 0.06
Nodes (55): ArgumentParser, BenchmarkCaseResult, BenchmarkManifestError, BenchmarkRunResult, _bool_option(), _case_metrics(), _expected_metrics(), _file_spec() (+47 more)

### Community 22 - "Community 22"
Cohesion: 0.50
Nodes (4): 1.1 Included data, 1.2 Excluded or deliberately non-flattened information, 1.3 Interpretation of a Boltz contact, 1. Scope and intended use

### Community 23 - "Community 23"
Cohesion: 0.50
Nodes (4): 9.1 Lower than 2 A, 9.2 Greater than 20 A, 9.3 Precision, 9. Boltz compatibility and bound handling

### Community 24 - "Community 24"
Cohesion: 0.33
Nodes (5): Constraint output split, Current result, Fail-closed gate, nmr2boltz paired-format benchmark, Row-level format discrepancy audit

### Community 30 - "Community 30"
Cohesion: 0.07
Nodes (74): Any, _alternative(), _constraint(), _endpoint(), _gate_run(), _report(), _semantic(), test_arbitrary_rejection_is_not_an_expected_format_difference() (+66 more)

## Knowledge Gaps
- **160 isolated node(s):** `nmr2boltz`, `warning`, `seed`, `requested`, `generated` (+155 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **3 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `TopologyLibrary` connect `Community 0` to `Community 1`, `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 10`, `Community 21`, `Community 30`?**
  _High betweenness centrality (0.059) - this node is a cross-community bridge._
- **Why does `project_document()` connect `Community 0` to `Community 3`, `Community 5`, `Community 6`, `Community 10`, `Community 21`, `Community 30`?**
  _High betweenness centrality (0.030) - this node is a cross-community bridge._
- **Why does `ProjectionSettings` connect `Community 0` to `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 10`, `Community 21`, `Community 30`?**
  _High betweenness centrality (0.025) - this node is a cross-community bridge._
- **Are the 17 inferred relationships involving `TopologyLibrary` (e.g. with `BenchmarkCaseResult` and `BenchmarkManifestError`) actually correct?**
  _`TopologyLibrary` has 17 INFERRED edges - model-reasoned connections that need verification._
- **Are the 19 inferred relationships involving `ProjectionSettings` (e.g. with `BenchmarkCaseResult` and `BenchmarkManifestError`) actually correct?**
  _`ProjectionSettings` has 19 INFERRED edges - model-reasoned connections that need verification._
- **What connects `nmr2boltz`, `Conservative NMR distance-restraint projection for Boltz atom contacts.`, `Raised when a benchmark manifest does not satisfy the versioned schema.` to the rest of the system?**
  _192 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.10708898944193061 - nodes in this community are weakly interconnected._