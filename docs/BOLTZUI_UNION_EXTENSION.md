# Design proposal: ambiguity-aware `atom_contact_union` constraints in BoltzUI/Boltz-2

## Problem

An NMR restraint often means:

```text
(A:10:CG1 -- A:42:CB <= 7.0 A)
OR
(A:10:CG2 -- A:42:CB <= 7.0 A)
```

The existing `atom_contact` schema represents one exact pair. Writing both alternatives as two entries changes OR into AND. This is especially harmful for aromatic symmetry, non-stereospecific assignments, and ambiguous NOE peak assignments.

The Boltz potential layer already has the useful abstraction of a `union_index`: multiple terms with the same union index are combined as one soft alternative group. The clean extension is therefore to expose that grouping in the BoltzUI input schema and preserve it through preprocessing.

## Proposed YAML schema

```yaml
constraints:
  - atom_contact_union:
      alternatives:
        - atom1: [A, 10, CG1]
          atom2: [A, 42, CB]
          max_distance: 7.0
        - atom1: [A, 10, CG2]
          atom2: [A, 42, CB]
          max_distance: 7.0
      force: true
      union_strength: 5.0       # optional; maps to union_lambda
      label: noe_list_1_187     # optional provenance label
```

Each alternative may have its own threshold. The parser should reject:

- fewer than two alternatives;
- an alternative whose two selectors resolve to the same atom;
- duplicate alternatives after canonical ordering;
- values outside the accepted atom-contact range;
- missing `force: true` if that remains a policy requirement;
- unsupported Boltz versions.

Normal `atom_contact` entries remain backward compatible.

## Suggested internal model

```python
@dataclass(frozen=True)
class AtomContactAlternative:
    atom1_index: int
    atom2_index: int
    max_distance: float

@dataclass(frozen=True)
class AtomContactUnion:
    alternatives: tuple[AtomContactAlternative, ...]
    union_strength: float
    label: str | None = None
```

Canonicalize each atom pair so `(i,j)` and `(j,i)` are identical. Keep a stable group ID for logs and provenance.

## Potential construction

For each union group:

1. allocate one new `union_index`;
2. append every alternative atom pair to the potential pair array;
3. append its individual threshold to the threshold array;
4. append the same `union_index` for every alternative;
5. use the requested/default `union_lambda` in the union aggregation.

Conceptual pseudocode:

```python
union_id = next_union_id()
for alternative in group.alternatives:
    pair_indices.append([alternative.atom1_index, alternative.atom2_index])
    thresholds.append(alternative.max_distance)
    union_indices.append(union_id)
```

The union energy should behave like a soft minimum across alternative violation energies. With sufficiently large positive `union_lambda`, one well-satisfied alternative dominates. The value should be exposed cautiously because an overly sharp soft minimum can cause unstable gradients, while an overly diffuse value can let several poor alternatives average together.

## Critical token-conditioning rule

Do **not** set the binary token-contact feature to 1 for every alternative in an OR group. That feature communicates that all listed token pairs are contacts, recreating AND semantics in the conditioning path even if the coordinate potential is union-aware.

Safe initial choices are:

1. apply the union only in the coordinate potential and leave token-contact conditioning unset for ambiguous groups; or
2. add a new group-aware probabilistic conditioning representation rather than reusing the binary contact matrix.

Option 1 is the minimal scientifically correct implementation.

## Multiple thresholds and duplicate parents

Before constructing the union:

- if multiple alternatives resolve to the same exact atom pair, replace them with one alternative using the **maximum** threshold, because they are inside one OR group;
- if independent union/contact groups later target the same pair, their effects remain conjunctive and should not be merged using the same rule.

## Optional weights

A future schema could include assignment probabilities:

```yaml
alternatives:
  - atom1: [A, 10, CG1]
    atom2: [A, 42, CB]
    max_distance: 7.0
    prior: 0.7
  - atom1: [A, 10, CG2]
    atom2: [A, 42, CB]
    max_distance: 7.0
    prior: 0.3
```

This should not be implemented by treating arbitrary NEF/NMR-STAR `weight` fields as probabilities automatically. The semantic meaning of weight is format/list dependent. Priors should be normalized, validated, and explicitly labeled as assignment probabilities.

## UI proposal

For each union group, display:

- group label/source restraint ID;
- alternatives as editable rows;
- maximum distance per row;
- an OR badge between alternatives;
- a warning that the token contact feature is not applied to ambiguous groups;
- an expandable provenance panel from the conversion report.

Provide three import modes:

1. **exact only**: import `atom_constraints_exact.yaml`;
2. **safe + unions**: import union groups using the new schema;
3. **hypothesis batch**: create one Boltz job per sampled assignment hypothesis.

## Validation tests

### Parser tests

- two valid alternatives resolve to one union group;
- invalid chain/residue/atom selector reports the exact alternative and group;
- same-atom alternatives are rejected;
- duplicate reversed pairs are deduplicated;
- individual thresholds are retained;
- exact `atom_contact` behavior is unchanged.

### Feature tests

- all alternatives share one union index;
- two separate groups receive distinct union indices;
- ambiguous alternatives do not all set token-contact conditioning;
- potential arrays remain aligned after batching/padding;
- thresholds survive device transfer and dtype conversion.

### Energy tests

Construct two alternatives where only one is within threshold and verify that:

- the union energy is low;
- enforcing both as separate contacts gives a higher energy;
- swapping alternative order does not change energy;
- duplicate alternatives do not change the result after deduplication;
- the gradient points toward satisfying the nearest plausible alternative.

### End-to-end tests

- parse one `atom_contact` plus one union group from YAML;
- run a minimal Boltz-2 inference job with potentials enabled;
- verify no token-conditioning AND leak;
- preserve provenance in output logs;
- compare a synthetic two-fold system under unrestrained, incorrect-AND, and union-aware guidance.

## Batch hypothesis fallback

Until this extension exists, `nmr2boltz --hypotheses N` is the safest compatible workaround. Each job contains one chosen pair per OR group. The result set should be pooled and scored against the original restraints. This costs more inference runs but does not require modifying Boltz internals.

## Recommended implementation order

1. Add schema/parser/dataclass support with no UI.
2. Route alternatives to the shared `union_index` in the potential layer.
3. Explicitly disable binary token-contact conditioning for union groups.
4. Add unit and synthetic energy tests.
5. Add CLI/import support for `atom_constraints_union.yaml`.
6. Add UI editing and provenance.
7. Evaluate `union_lambda` and potential strength on blinded NMR-guided folding cases.
8. Consider weighted priors only after the unweighted behavior is stable.

## Scientific acceptance criterion

The extension is successful only if one satisfied alternative can guide the structure without simultaneously attracting every alternative pair. Lower loss alone is not sufficient; the implementation must preserve the source OR semantics in every conditioning and potential pathway.
