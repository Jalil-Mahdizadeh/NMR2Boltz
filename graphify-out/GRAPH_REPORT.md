# Graph Report - .  (2026-07-15)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 118 nodes · 410 edges · 13 communities
- Extraction: 90% EXTRACTED · 10% INFERRED · 0% AMBIGUOUS · INFERRED: 43 edges (avg confidence: 0.51)
- Token cost: 0 input · 0 output

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

## God Nodes (most connected - your core abstractions)
1. `ProjectionSettings` - 19 edges
2. `ComponentTopology` - 19 edges
3. `TopologyLibrary` - 17 edges
4. `project_document()` - 16 edges
5. `clean()` - 15 edges
6. `BoltzAtom` - 13 edges
7. `StarDataError` - 13 edges
8. `SequenceResolver` - 13 edges
9. `parse_star_document()` - 13 edges
10. `RestraintGroup` - 12 edges

## Surprising Connections (you probably didn't know these)
- `FlowList` --uses--> `BoltzAtom`  [INFERRED]
  src/nmr2boltz/output.py → src/nmr2boltz/model.py
- `ProjectionSettings` --uses--> `BoltzAtom`  [INFERRED]
  src/nmr2boltz/project.py → src/nmr2boltz/model.py
- `ParsedStarDocument` --uses--> `SequenceRecord`  [INFERRED]
  src/nmr2boltz/star.py → src/nmr2boltz/model.py
- `ParsedStarDocument` --uses--> `Endpoint`  [INFERRED]
  src/nmr2boltz/star.py → src/nmr2boltz/model.py
- `SequenceResolver` --uses--> `Endpoint`  [INFERRED]
  src/nmr2boltz/star.py → src/nmr2boltz/model.py

## Import Cycles
- None detected.

## Communities (13 total, 0 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.28
Nodes (13): Rejection, RestraintGroup, _averaging_factor(), _canonical_pair(), _merge_independent_constraints(), _merge_or_alternatives(), project_document(), _project_group() (+5 more)

### Community 1 - "Community 1"
Cohesion: 0.26
Nodes (7): Conservative NMR distance-restraint projection for Boltz atom contacts., AmbiguousGroup, BoltzAtom, EmittedConstraint, ProjectedAlternative, BoltzYamlDumper, _combine_hypothesis_constraints()

### Community 2 - "Community 2"
Cohesion: 0.33
Nodes (9): AtomChoice, AtomSetChoice, _compile_atom_pattern(), _leading_digit_alias(), _looks_like_hydrogen(), _normalize_legacy_expression(), _pseudoatom_to_expression(), TopologyResolutionError (+1 more)

### Community 3 - "Community 3"
Cohesion: 0.38
Nodes (3): _optional_float(), TopologyLibrary, Path

### Community 4 - "Community 4"
Cohesion: 0.47
Nodes (10): Any, clean(), Return a stripped STAR value or ``None`` for STAR null tokens., _category(), extract_embedded_topologies(), extract_restraint_groups(), loop_rows(), pick() (+2 more)

### Community 5 - "Community 5"
Cohesion: 0.47
Nodes (8): SequenceRecord, _extract_nef_sequence(), _extract_nmrstar_sequence(), extract_sequence_map(), overlay_residue_map(), SequenceResolver, StarDataError, ValueError

### Community 6 - "Community 6"
Cohesion: 0.47
Nodes (7): ConversionReport, asdict_atom(), _dump_yaml(), _summary_text(), write_hypotheses(), write_outputs(), _write_tsv()

### Community 7 - "Community 7"
Cohesion: 0.39
Nodes (6): ArgumentParser, Namespace, build_parser(), command_convert(), main(), parse_chain_map()

### Community 8 - "Community 8"
Cohesion: 0.36
Nodes (6): as_float(), Endpoint, RawAlternative, _parse_nef_row(), _parse_nmrstar_row(), _select_upper_bound()

### Community 9 - "Community 9"
Cohesion: 0.38
Nodes (3): build_builtin_topologies(), ComponentTopology, infer_element()

### Community 10 - "Community 10"
Cohesion: 0.40
Nodes (5): list, _flow_list_representer(), FlowList, Node, SafeDumper

### Community 11 - "Community 11"
Cohesion: 0.70
Nodes (4): infer_sequence_map(), _natural_sequence_key(), _parse_entry_file(), parse_star_document()

### Community 12 - "Community 12"
Cohesion: 0.67
Nodes (4): _atom_list(), _constraint_payload(), _rounded(), _union_payload()

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `ComponentTopology` connect `Community 9` to `Community 0`, `Community 2`, `Community 3`, `Community 4`, `Community 5`, `Community 11`?**
  _High betweenness centrality (0.084) - this node is a cross-community bridge._
- **Why does `TopologyLibrary` connect `Community 3` to `Community 0`, `Community 2`, `Community 7`?**
  _High betweenness centrality (0.068) - this node is a cross-community bridge._
- **Why does `clean()` connect `Community 4` to `Community 1`, `Community 2`, `Community 3`, `Community 5`, `Community 8`, `Community 9`, `Community 11`?**
  _High betweenness centrality (0.058) - this node is a cross-community bridge._
- **Are the 11 inferred relationships involving `ProjectionSettings` (e.g. with `AmbiguousGroup` and `BoltzAtom`) actually correct?**
  _`ProjectionSettings` has 11 INFERRED edges - model-reasoned connections that need verification._
- **Are the 5 inferred relationships involving `ComponentTopology` (e.g. with `ParsedStarDocument` and `SequenceResolver`) actually correct?**
  _`ComponentTopology` has 5 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `TopologyLibrary` (e.g. with `ProjectionSettings` and `AtomChoice`) actually correct?**
  _`TopologyLibrary` has 3 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Conservative NMR distance-restraint projection for Boltz atom contacts.`, `Return a stripped STAR value or ``None`` for STAR null tokens.`, `Merge duplicate heavy pairs inside one OR group.      For the same pair, ``d <=` to the rest of the system?**
  _4 weakly-connected nodes found - possible documentation gaps or missing edges._