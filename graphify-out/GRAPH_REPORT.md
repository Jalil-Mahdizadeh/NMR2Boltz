# Graph Report - nmr2boltz_toolkit  (2026-07-19)

## Corpus Check
- 128 files · ~17,606,713 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 617 nodes · 1689 edges · 32 communities (28 shown, 4 thin omitted)
- Extraction: 93% EXTRACTED · 7% INFERRED · 0% AMBIGUOUS · INFERRED: 115 edges (avg confidence: 0.51)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `29a55b14`
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
- [[_COMMUNITY_Community 25|Community 25]]
- [[_COMMUNITY_Community 26|Community 26]]
- [[_COMMUNITY_Community 27|Community 27]]
- [[_COMMUNITY_Community 29|Community 29]]
- [[_COMMUNITY_Community 30|Community 30]]
- [[_COMMUNITY_Community 31|Community 31]]

## God Nodes (most connected - your core abstractions)
1. `TopologyLibrary` - 71 edges
2. `ProjectionSettings` - 53 edges
3. `project_document()` - 44 edges
4. `Scientific and computational method for projecting proton NMR distance restraints onto Boltz heavy-atom contacts` - 41 edges
5. `parse_star_document()` - 38 edges
6. `BoltzAtom` - 33 edges
7. `write_outputs()` - 33 edges
8. `ConversionReport` - 32 edges
9. `SequenceRecord` - 26 edges
10. `TopologyResolutionError` - 26 edges

## Surprising Connections (you probably didn't know these)
- `test_yaml_rounding_is_outward_for_upper_bounds()` --calls--> `_rounded()`  [EXTRACTED]
  tests/test_robustness.py → src/nmr2boltz/output.py
- `_Saveframe` --uses--> `StarDataError`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/star.py
- `_StarLoop` --uses--> `StarDataError`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/star.py
- `test_compressed_nef_input()` --calls--> `parse_star_document()`  [EXTRACTED]
  tests/test_logic.py → src/nmr2boltz/star.py
- `_Saveframe` --uses--> `TopologyResolutionError`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/topology.py

## Import Cycles
- None detected.

## Communities (32 total, 4 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.13
Nodes (37): command_convert(), project_document(), ProjectionSettings, parse_star_document(), ParsedStarDocument, TopologyLibrary, convert(), test_embedded_custom_component_topology() (+29 more)

### Community 1 - "Community 1"
Cohesion: 0.15
Nodes (26): Atom, AtomKey, Point, test_structure_alignment_accepts_equivalent_nucleic_acid_names(), test_structure_alignment_rekeys_author_numbers_to_boltz_indices(), _align_components(), align_structure_to_sequence_map(), atom_key() (+18 more)

### Community 2 - "Community 2"
Cohesion: 0.05
Nodes (37): 10.1 Safe default, 10.2 Union-aware potential, 10.3 Assignment hypotheses, 10.4 Inter-chain-only selection, 10. Ambiguous groups and recommended execution strategies, 12. Suggested confidence tiers, 14. Reproducibility and auditability, 15. Questions for NMR expert review before production deployment (+29 more)

### Community 3 - "Community 3"
Cohesion: 0.10
Nodes (44): Decimal, float, list, AmbiguousGroup, ConversionReport, EmittedConstraint, asdict_atom(), _atom_list() (+36 more)

### Community 4 - "Community 4"
Cohesion: 0.11
Nodes (35): ProjectedAlternative, Rejection, RestraintGroup, _averaging_factor(), _canonical_pair(), _intrachain_filter_rejection(), _merge_independent_constraints(), _merge_or_alternatives() (+27 more)

### Community 5 - "Community 5"
Cohesion: 0.11
Nodes (27): AtomChoice, AtomSetChoice, clean(), Return a stripped STAR value or ``None`` for STAR null tokens., atom_topology_violations(), AtomTopologyValidationError, AtomTopologyViolation, build_builtin_topologies() (+19 more)

### Community 6 - "Community 6"
Cohesion: 0.16
Nodes (33): as_float(), Endpoint, RawAlternative, SequenceRecord, _alternative_row_key(), _canonical_expansion_family_key(), _canonical_expansion_issue(), _canonical_expansion_semantics_key() (+25 more)

### Community 7 - "Community 7"
Cohesion: 0.08
Nodes (33): ArgumentParser, Namespace, build_parser(), command_benchmark(), main(), parse_chain_map(), Conservative NMR distance-restraint projection for Boltz atom contacts., BoltzAtom (+25 more)

### Community 8 - "Community 8"
Cohesion: 0.09
Nodes (22): 1. Data sources reviewed, 2. NEF-to-NMR-STAR field correspondence used by the converter, 3.1 Canonical expansion of one author atom set, 3.2 Aromatic and symmetric ambiguity, 3.3 Insertion-like sequence codes, 3.4 Missing bounds, 3.5 Modified components, 3.6 Complex logical combinations (+14 more)

### Community 9 - "Community 9"
Cohesion: 0.11
Nodes (17): Batch hypothesis fallback, Critical token-conditioning rule, Design proposal: ambiguity-aware `atom_contact_union` constraints in BoltzUI/Boltz-2, End-to-end tests, Energy tests, Feature tests, Multiple thresholds and duplicate parents, Optional weights (+9 more)

### Community 10 - "Community 10"
Cohesion: 0.44
Nodes (11): emitted_atom_topology_violations(), _document(), _project(), test_available_modified_ligand_and_ion_topology_emits(), test_deposited_corpus_atom_defects_are_quarantined(), test_final_output_validator_fails_before_writing_invalid_yaml(), test_glutamine_zinc_atom_is_quarantined_with_complete_provenance(), test_protein_residue_rejects_projected_nucleotide_atoms() (+3 more)

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
Cohesion: 0.29
Nodes (6): Ambiguous restraints, BoltzUI integration guide, Exact-contact schema, Inter-chain-only constraints, Recommended run protocol, Residue numbering

### Community 18 - "Community 18"
Cohesion: 0.29
Nodes (7): 11.1 Parsing validation, 11.2 Identifier validation, 11.3 Geometric validation of candidate structures, 11.4 Structural and statistical controls, 11.5 Executed validation record, 11.5 Paired-format discrepancy audit and CI gate, 11. Validation protocol

### Community 19 - "Community 19"
Cohesion: 0.33
Nodes (6): 13.1 Probabilistic union restraints, 13.2 Ensemble-aware objective, 13.3 Iterative assignment and folding, 13.4 Information-aware restraint selection, 13.5 Broader experimental data, 13. Potential methodological extensions

### Community 20 - "Community 20"
Cohesion: 0.33
Nodes (5): Included, nmr2boltz toolkit release notes, Safety defaults, Unreleased, Validation

### Community 21 - "Community 21"
Cohesion: 0.18
Nodes (23): BenchmarkCaseResult, BenchmarkManifestError, BenchmarkRunResult, _bool_option(), _case_metrics(), _expected_metrics(), _file_spec(), FileSpec (+15 more)

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
Nodes (79): Any, _parse_entry_file(), Path, _alternative(), _constraint(), _endpoint(), _gate_run(), _report() (+71 more)

### Community 31 - "Community 31"
Cohesion: 0.44
Nodes (12): _cartesian_rows(), _endpoint(), _nef_group(), _report(), _star_group(), test_equivalent_nef_and_star_exact_and_union_yaml_match(), test_incomplete_or_expansion_is_rejected_fail_closed(), test_incomplete_two_sided_cartesian_expansion_is_rejected() (+4 more)

## Knowledge Gaps
- **150 isolated node(s):** `nmr2boltz`, `warning`, `seed`, `requested`, `generated` (+145 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **4 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `TopologyLibrary` connect `Community 0` to `Community 1`, `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 7`, `Community 10`, `Community 21`, `Community 30`, `Community 31`?**
  _High betweenness centrality (0.063) - this node is a cross-community bridge._
- **Why does `project_document()` connect `Community 0` to `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 7`, `Community 10`, `Community 21`, `Community 30`, `Community 31`?**
  _High betweenness centrality (0.027) - this node is a cross-community bridge._
- **Why does `ProjectionSettings` connect `Community 0` to `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 7`, `Community 10`, `Community 21`, `Community 30`, `Community 31`?**
  _High betweenness centrality (0.024) - this node is a cross-community bridge._
- **Are the 17 inferred relationships involving `TopologyLibrary` (e.g. with `BenchmarkCaseResult` and `BenchmarkManifestError`) actually correct?**
  _`TopologyLibrary` has 17 INFERRED edges - model-reasoned connections that need verification._
- **Are the 19 inferred relationships involving `ProjectionSettings` (e.g. with `BenchmarkCaseResult` and `BenchmarkManifestError`) actually correct?**
  _`ProjectionSettings` has 19 INFERRED edges - model-reasoned connections that need verification._
- **What connects `nmr2boltz`, `Conservative NMR distance-restraint projection for Boltz atom contacts.`, `Raised when a benchmark manifest does not satisfy the versioned schema.` to the rest of the system?**
  _179 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.12956810631229235 - nodes in this community are weakly interconnected._