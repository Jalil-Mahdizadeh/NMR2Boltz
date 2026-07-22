# Graph Report - nmr2boltz_toolkit  (2026-07-22)

## Corpus Check
- 127 files · ~16,649,409 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 728 nodes · 2175 edges · 26 communities (24 shown, 2 thin omitted)
- Extraction: 93% EXTRACTED · 7% INFERRED · 0% AMBIGUOUS · INFERRED: 142 edges (avg confidence: 0.51)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `c0017a73`
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
- [[_COMMUNITY_Community 13|Community 13]]
- [[_COMMUNITY_Community 14|Community 14]]
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
- [[_COMMUNITY_Community 29|Community 29]]
- [[_COMMUNITY_Community 30|Community 30]]

## God Nodes (most connected - your core abstractions)
1. `TopologyLibrary` - 80 edges
2. `ProjectionSettings` - 63 edges
3. `project_document()` - 55 edges
4. `BoltzAtom` - 51 edges
5. `write_outputs()` - 45 edges
6. `parse_star_document()` - 45 edges
7. `Scientific and computational method for projecting proton NMR distance restraints onto Boltz heavy-atom contacts` - 42 edges
8. `ConversionReport` - 39 edges
9. `SequenceRecord` - 37 edges
10. `TopologyResolutionError` - 26 edges

## Surprising Connections (you probably didn't know these)
- `_Saveframe` --uses--> `SequenceRecord`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/model.py
- `_StarLoop` --uses--> `SequenceRecord`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/model.py
- `_Saveframe` --uses--> `Endpoint`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/model.py
- `_StarLoop` --uses--> `Endpoint`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/model.py
- `_Saveframe` --uses--> `SequenceResolver`  [INFERRED]
  tests/test_robustness.py → src/nmr2boltz/star.py

## Import Cycles
- None detected.

## Communities (26 total, 2 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.07
Nodes (66): ArgumentParser, Namespace, build_parser(), command_benchmark(), command_convert(), main(), parse_chain_map(), project_document() (+58 more)

### Community 1 - "Community 1"
Cohesion: 0.09
Nodes (60): Atom, AtomKey, _parse_entry_file(), PairKey, Path, Point, _atom(), test_build_case_writes_union_of_exact_pairs_with_model_distances() (+52 more)

### Community 2 - "Community 2"
Cohesion: 0.05
Nodes (38): 10.1 Safe default, 10.2 Union-aware potential, 10.3 Assignment hypotheses, 10.4 Intraresidue exclusion, 10.5 Inter-chain-only selection, 10. Ambiguous groups and recommended execution strategies, 12. Suggested confidence tiers, 14. Reproducibility and auditability (+30 more)

### Community 3 - "Community 3"
Cohesion: 0.07
Nodes (65): Decimal, float, list, Conservative NMR distance-restraint projection for Boltz atom contacts., AmbiguousGroup, BoltzToken, ConversionReport, EmittedConstraint (+57 more)

### Community 4 - "Community 4"
Cohesion: 0.32
Nodes (13): Random, _add(), _distance(), main(), _projected(), Deterministic mathematical and converter stress validation for nmr2boltz.  This, _unit_vector(), validate_boolean_merge_order() (+5 more)

### Community 5 - "Community 5"
Cohesion: 0.09
Nodes (32): AtomChoice, AtomSetChoice, atom_topology_violations(), AtomTopologyValidationError, AtomTopologyViolation, build_builtin_topologies(), _compile_atom_pattern(), component_topology_snapshot() (+24 more)

### Community 6 - "Community 6"
Cohesion: 0.10
Nodes (53): as_float(), clean(), Endpoint, Return a stripped STAR value or ``None`` for STAR null tokens., RawAlternative, SequenceRecord, _alternative_row_key(), _canonical_expansion_family_key() (+45 more)

### Community 7 - "Community 7"
Cohesion: 0.22
Nodes (8): `benchmark`, Command structure, `convert`, `convert` effective-default notes, `convert` output bundle, Exit statuses, Global flags, NMR2Boltz command-line reference

### Community 8 - "Community 8"
Cohesion: 0.09
Nodes (22): 1. Data sources reviewed, 2. NEF-to-NMR-STAR field correspondence used by the converter, 3.1 Canonical expansion of one author atom set, 3.2 Aromatic and symmetric ambiguity, 3.3 Insertion-like sequence codes, 3.4 Missing bounds, 3.5 Modified components, 3.6 Complex logical combinations (+14 more)

### Community 9 - "Community 9"
Cohesion: 0.11
Nodes (17): Batch hypothesis fallback, Critical token-conditioning rule, Design proposal: ambiguity-aware `atom_contact_union` constraints in BoltzUI/Boltz-2, End-to-end tests, Energy tests, Feature tests, Multiple thresholds and duplicate parents, Optional weights (+9 more)

### Community 10 - "Community 10"
Cohesion: 0.30
Nodes (20): BoltzAtom, _alternative(), _exact(), _project(), _report(), _sequence(), test_above_maximum_bound_is_omitted_and_never_clipped(), test_collapsible_union_above_maximum_is_omitted_as_one_semantic_unit() (+12 more)

### Community 11 - "Community 11"
Cohesion: 0.08
Nodes (21): Paired NEF/NMR-STAR benchmark, Exact-contact distance check, Token-contact distance check, Paired-format discrepancy findings, Ambiguity rules, BoltzUI integration, Docker, Documentation (+13 more)

### Community 13 - "Community 13"
Cohesion: 0.22
Nodes (8): Atom semantics, Boltz usage, Identifier mapping, NMR expert review checklist, Post-prediction validation, Projection, Sign-off, Source and calibration

### Community 14 - "Community 14"
Cohesion: 0.39
Nodes (4): inline_markup(), main(), parse_markdown(), table_from_markdown()

### Community 17 - "Community 17"
Cohesion: 0.22
Nodes (8): Ambiguous restraints, BoltzUI integration guide, Exact-contact schema, Inter-chain-only constraints, Intraresidue-excluded constraints, Recommended run protocol, Residue numbering, Standalone token-contact schema

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
Cohesion: 0.09
Nodes (45): BenchmarkCaseResult, BenchmarkManifestError, BenchmarkRunResult, _bool_option(), _case_metrics(), _expected_metrics(), _file_spec(), FileSpec (+37 more)

### Community 22 - "Community 22"
Cohesion: 0.50
Nodes (4): 1.1 Included data, 1.2 Excluded or deliberately non-flattened information, 1.3 Interpretation of a Boltz contact, 1. Scope and intended use

### Community 23 - "Community 23"
Cohesion: 0.40
Nodes (5): 9.1 Lower than 2 A, 9.2 Greater than 20 A, 9.3 Precision, 9.4 Token-domain projection, 9. Boltz compatibility and bound handling

### Community 24 - "Community 24"
Cohesion: 0.33
Nodes (5): Constraint output split, Current result, Fail-closed gate, nmr2boltz paired-format benchmark, Row-level format discrepancy audit

### Community 25 - "Community 25"
Cohesion: 0.11
Nodes (40): ProjectedAlternative, Rejection, RestraintGroup, _adapt_ambiguous_group_bounds(), _averaging_factor(), _canonical_pair(), _intrachain_filter_rejection(), _intraresidue_filter_rejection() (+32 more)

### Community 30 - "Community 30"
Cohesion: 0.07
Nodes (78): Any, SimpleNamespace, _alternative(), _constraint(), _endpoint(), _gate_run(), _report(), _semantic() (+70 more)

## Knowledge Gaps
- **142 isolated node(s):** `nmr2boltz`, `TokenContribution`, `TokenProjectionOmission`, `Installation`, `Docker` (+137 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **2 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `TopologyLibrary` connect `Community 0` to `Community 1`, `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 10`, `Community 21`, `Community 25`, `Community 30`?**
  _High betweenness centrality (0.064) - this node is a cross-community bridge._
- **Why does `BoltzAtom` connect `Community 10` to `Community 0`, `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 21`, `Community 25`, `Community 30`?**
  _High betweenness centrality (0.041) - this node is a cross-community bridge._
- **Why does `project_document()` connect `Community 0` to `Community 3`, `Community 5`, `Community 6`, `Community 21`, `Community 25`, `Community 30`?**
  _High betweenness centrality (0.032) - this node is a cross-community bridge._
- **Are the 17 inferred relationships involving `TopologyLibrary` (e.g. with `BenchmarkCaseResult` and `BenchmarkManifestError`) actually correct?**
  _`TopologyLibrary` has 17 INFERRED edges - model-reasoned connections that need verification._
- **Are the 19 inferred relationships involving `ProjectionSettings` (e.g. with `BenchmarkCaseResult` and `BenchmarkManifestError`) actually correct?**
  _`ProjectionSettings` has 19 INFERRED edges - model-reasoned connections that need verification._
- **What connects `nmr2boltz`, `Conservative NMR distance-restraint projection for Boltz atom contacts.`, `Raised when a benchmark manifest does not satisfy the versioned schema.` to the rest of the system?**
  _182 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.06783511846802986 - nodes in this community are weakly interconnected._