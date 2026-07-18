# Scientific and computational method for projecting proton NMR distance restraints onto Boltz heavy-atom contacts

## Executive summary

The practical problem is a mismatch of representations. Nuclear magnetic resonance distance information is commonly expressed between protons, or between one proton and another nucleus, while the present Boltz folding workflow does not explicitly model hydrogens. A useful conversion must therefore answer three separate questions:

1. **Which physical nuclei does the deposited restraint represent?**
2. **Which modeled heavy atom is directly bonded to each omitted proton?**
3. **What heavy-heavy upper bound follows without introducing a false geometric assumption?**

The conservative answer for one explicit proton pair is the triangle-inequality bound

\[
D(P_1,P_2) \le U(H_1,H_2) + L(P_1-H_1) + L(P_2-H_2) + m,
\]

where `P1/P2` are directly bonded heavy parents, `U` is the source proton-distance upper bound, `L` is a conservative covalent-bond-length upper envelope, and `m` is optional user-selected slack. For a proton-heavy restraint, only one X-H term is added. For a heavy-heavy restraint, no X-H term is added.

This formula is rigorous for an explicit pair, but many NMR restraint rows are
not explicit pairs. They may denote an r^-6 atom set, a non-stereospecific
assignment, a peak-assignment alternative, or a geometric pseudoatom. The
conversion therefore preserves **OR logic**. It emits a normal Boltz
`atom_contact` only when every alternative reduces to the same heavy-atom pair;
otherwise a fully projectable group is written as one metadata-free
`atom_contact_union`. Unsupported or partially projectable groups remain
rejected, and explicit assignment hypotheses are available as a separate
model-selection workflow.

The method is intentionally one-sided: it is designed not to make an experimental upper bound tighter accidentally. It does not claim that the resulting heavy-atom contact is information-equivalent to the original NMR restraint.

Validation record: this document describes toolkit version 0.1.0 as audited on 2026-07-17. The executed software and mathematical validation is reported in Section 11.5. That validation establishes converter behavior, paired-format differences, and target-schema compatibility; it does not establish the empirical accuracy of Boltz predictions on experimental benchmark structures.

---

## 1. Scope and intended use

### 1.1 Included data

The implementation reads:

- NEF distance-restraint saveframes containing `_nef_distance_restraint` loops;
- NMR-STAR general-distance-constraint saveframes containing `_Gen_dist_constraint` loops;
- molecular-system/sequence loops needed to map source identifiers to one-based Boltz sequence positions;
- NMR-STAR chemical-component atom and bond loops when present;
- optional local wwPDB Chemical Component Dictionary mmCIF files for modified residues and ligands.

The principal target is upper-bound information derived from NOE-like measurements, but the parser can retain any general distance restraint. An origin/type filter can be applied at the command line.

### 1.2 Excluded or deliberately non-flattened information

The current implementation does not attempt to translate:

- lower bounds into Boltz atom contacts;
- full restraint potential shapes, force constants, or linear-tail behavior;
- NOE intensity calibration itself;
- spin diffusion corrections;
- conformational-ensemble or time averaging as a learned ensemble potential;
- nested restraint-combination logic with non-null combination identifiers or non-OR NMR-STAR member logic;
- geometric pseudoatoms without an explicit opt-in approximation;
- probabilistic chemical-shift/peak assignment networks.

These are preserved in provenance or reported as unresolved rather than guessed.

### 1.3 Interpretation of a Boltz contact

A Boltz distance potential is guidance during generative structure prediction. It is not a proof that the final coordinates obey a hard bound, and it is not a substitute for standard NMR restraint validation. The recommended workflow is to generate an ensemble of candidate structures, reconstruct hydrogens, and score the original restraints after prediction.

---

## 2. Source data models

## 2.1 NEF

A typical NEF distance-restraint list contains saveframe-level metadata and a loop with fields such as:

- `restraint_id` and optional `restraint_combination_id`;
- chain, sequence, residue, and atom identifiers for endpoint 1 and endpoint 2;
- weight, target value, uncertainty, lower limit, upper limit, and linear limits.

Important semantic properties are:

1. **`sequence_code` is a string.** It can contain insertion-like or otherwise non-integer identifiers. It must not be converted directly to a Boltz residue index.
2. **Rows with the same `restraint_id` are alternatives.** They normally represent an OR relation, not independent contacts that should all be enforced.
3. **`%` and `*` are atom-name wildcards.** `%` denotes a sequence of digits and `*` a general non-whitespace string.
4. **`x` and `y` denote non-stereospecific assignments.** Their consistency may be global across a dataset; treating each occurrence independently is only a local representation.
5. **IUPAC Q/M pseudoatoms are geometric pseudoatoms.** They are not synonyms for wildcard atom sets and need their own correction convention.
6. **Identifiers are formally case-sensitive.** The software permits a unique case-insensitive fallback only when it is unambiguous and records a warning.

The `_nef_sequence` loop is used in file order to assign one-based positions within each chain. Thus a source sequence `24`, `24A`, `25` can map to Boltz positions `24`, `25`, `26` if that is the actual sequence order.

## 2.2 NMR-STAR

The NMR-STAR counterpart is usually:

- saveframe category `general_distance_constraints`;
- loop category `_Gen_dist_constraint`;
- source author identifiers (`Auth_*`) and translated canonical identifiers (`Entity_assembly_ID`, `Comp_index_ID`, `Comp_ID`, `Atom_ID`);
- `Distance_upper_bound_val` and related fields;
- member IDs and logic codes for alternatives.

A translated NMR-STAR entry may expand one author-level atom set into several
canonical OR rows. For example, author-level `HG1%` or `HG1` can appear as
canonical `HG11`, `HG12`, and `HG13` rows with the same restraint ID and upper
bound. Before projection, the converter groups rows only when author endpoint
identities and all bounds, target, uncertainty, weight, origin, combination,
and explicit OR semantics agree. Component topology must prove that the
canonical atoms form the complete proton set on each heavy parent, and the
observed rows must form the complete one- or two-sided Cartesian product. The
reconstructed atom sets determine `N` once before heavy-parent projection.
Branches on different heavy parents remain disjunctive alternatives. Missing,
incomplete, inconsistent, or topology-unverified expansions are rejected in
full rather than inferred from lexical prefixes.

The NMR-STAR dictionary defines `Combination_ID` and `Member_logic_code` in addition to the constraint ID. The converter flattens only explicit `OR` member logic. Any non-null combination identifier, `AND` member logic, or unknown non-OR code marks the group as complex and prevents emission. This conservative rule avoids converting a structured Boolean expression into an accidental conjunction of Boltz contacts.

Residue mapping is built from `_Chem_comp_assembly` when available. `_Entity_poly_seq` plus `_Entity_assembly` is used as a fallback. Canonical component indices are accepted directly only when they form a unique contiguous sequence `1..N`; otherwise loop order is used and a warning is emitted.

## 2.3 Null values and optional fields

STAR null tokens `.` and `?` are treated as missing. An explicit upper bound is required by default. Optional policies can derive a heuristic value from:

- upper linear limit;
- target plus uncertainty;
- target alone.

Such a derived value is labeled in the report because it is not equivalent to a deposited hard upper limit.
An uncertainty used by the target-plus-uncertainty policy must be non-negative. Non-finite numeric values are treated as unavailable or rejected at the settings boundary.

---

## 3. Internal logical representation

Every source list is normalized into:

- a **restraint group**, identified by list/saveframe and restraint ID;
- one or more **source alternatives** within that group;
- two endpoint expressions per alternative;
- the selected upper bound and its provenance;
- source rows, origin/type, weights, member IDs, and warnings.

The group is interpreted as a disjunction:

\[
A_1 \lor A_2 \lor \cdots \lor A_K.
\]

This distinction is essential. Converting the group into K independent Boltz contacts would instead produce:

\[
A_1 \land A_2 \land \cdots \land A_K,
\]

which can be much stronger and physically false.

The software supports a safe simplification only when all alternatives ultimately target the same heavy-atom pair. If two alternatives for the same pair have thresholds `u1` and `u2`, then

\[
(D \le u_1) \lor (D \le u_2) \equiv D \le \max(u_1,u_2).
\]

By contrast, independent restraint groups are conjunctive. If two independent restraints target the same heavy pair,

\[
(D \le u_1) \land (D \le u_2) \equiv D \le \min(u_1,u_2).
\]

This max-inside-OR/min-across-AND rule is implemented explicitly.

Completeness is also required. If any source OR alternative, or any explicit pair within an atom-set branch, cannot be projected conservatively, the converter emits none of the remaining alternatives from that group. Keeping only the successful alternatives would narrow the disjunction and strengthen the source restraint.

---

## 4. Residue and chain mapping

## 4.1 Why source numbering cannot be copied blindly

The Boltz atom selector uses the chain ID and the one-based residue position in the exact input sequence. NMR files may use:

- author residue numbering with gaps;
- insertion codes;
- negative or zero numbering;
- entity/component indices;
- different chain aliases;
- separate author and canonical namespaces.

Therefore residue mapping is an independent, auditable transformation.

## 4.2 Mapping precedence

The converter resolves a source endpoint in this order:

1. exact author chain, sequence code, and residue name;
2. exact canonical assembly/component index and component name;
3. unique chain-plus-sequence match with a residue-name mismatch warning;
4. globally unique sequence match when the chain is absent;
5. explicit user residue-map override.

A user map can override extracted records. It should be mandatory in production whenever the sequence used for Boltz differs from the deposited molecular system.

## 4.3 Inference policy

If no sequence loop exists, inference from restraint identifiers is disabled by default. `--allow-inferred-sequence-map` sorts identifiers by numeric prefix and suffix within each chain, but every inferred record is marked. This option is convenient for exploration and unsafe for unattended production unless verified against the exact Boltz input.

---

## 5. Chemical topology and proton-parent mapping

## 5.1 Direct-bond rule

A proton is mapped only to the heavy atom directly bonded to it. Examples include:

- backbone `H/HN -> N`;
- `HA -> CA`;
- `VAL HG11/HG12/HG13 -> CG1`;
- `LYS HZ1/HZ2/HZ3 -> NZ`;
- `SER HG -> OG`;
- nucleic-acid `H1' -> C1'`, `H5'/H5'' -> C5'`, and base-specific mappings.

A hydrogen with no unique parent is not projected.

## 5.2 Topology sources

Topology is resolved in the following hierarchy:

1. chemical-component atom/bond loops embedded in the NMR-STAR file;
2. user-supplied wwPDB CCD mmCIF data;
3. a bundled atom inventory and proton-parent map for the 20 standard amino
   acids, common protonation variants, MSE/SEC, and common RNA/DNA residues.

Embedded and external bond tables are data-driven. A hydrogen is accepted only when it has a unique non-hydrogen neighbor.

## 5.3 Fail-closed atom membership

After sequence mapping and proton-to-heavy projection, but before OR/pair
deduplication, both projected atoms must occur in the exact topology of their
mapped residue/component. A missing atom or unavailable component topology
quarantines the complete contact as `atom_not_present_in_mapped_residue`; no
atom or residue is guessed or remapped. The rejection preserves source rows,
endpoint expressions, mapped chain/residue/component/atom, restraint group, and
original bounds.

The conversion report freezes the component atom dictionaries used for that
decision. Immediately before writing, an independent validator checks every
emitted endpoint against this target-topology snapshot and raises on any
violation before creating files. This is a chemical-topology invariant:
coordinate absence is neither necessary nor sufficient evidence that an atom is
invalid.

## 5.4 X-H upper envelopes

The default conservative envelopes are:

| Parent element | X-H upper envelope (A) |
|---|---:|
| C | 1.12 |
| N | 1.08 |
| O | 1.02 |
| S | 1.36 |
| Se | 1.49 |
| P | 1.45 |
| unknown | 1.55 |

When a CCD ideal distance and estimated standard deviation are available, the implementation uses the larger of the element default and `ideal + max(0.03, 3*esd)`.

These numbers are software defaults, not universal physical constants. They can be overridden with a JSON file:

```json
{
  "element_upper": {"C": 1.13, "N": 1.09},
  "components": {"CYS": {"HG": 1.37}}
}
```

An additional global non-negative `--projection-margin` may be used. Increasing a bond envelope or margin weakens the heavy-atom constraint and therefore preserves conservativeness.

## 5.5 Isotopes and exchangeable protons

Names beginning with D/T are interpreted as hydrogen isotopes when topology supports them. Exchangeable hydrogens attached to N/O/S/Se are topology-dependent. Their observation and assignment can depend on protonation, tautomer, pH, solvent, temperature, and exchange kinetics. A chemically valid parent mapping does not guarantee that the deposited atom was present in every conformer.

---

## 6. Mathematical projection for explicit pairs

Let the observed/modelled nuclei be `h1` and `h2`, with directly bonded modeled parents `p1` and `p2`. By the triangle inequality along the path `p1 -> h1 -> h2 -> p2`,

\[
\|p_1-p_2\| \le \|p_1-h_1\| + \|h_1-h_2\| + \|h_2-p_2\|.
\]

If the source states

\[
\|h_1-h_2\| \le U,
\]

and the covalent bond lengths obey

\[
\|p_1-h_1\| \le L_1, \qquad \|p_2-h_2\| \le L_2,
\]

then

\[
\boxed{\|p_1-p_2\| \le U + L_1 + L_2.}
\]

With explicit margin `m >= 0`, the emitted bound is

\[
U_{PP}=U+L_1+L_2+m.
\]

Special cases:

- H-H: add both bond envelopes;
- H-X, where X is already a modeled heavy atom: add one envelope;
- X-Y: add none;
- both hydrogens on the same heavy parent: the projected pair is the same atom and is not emitted.

The formula uses no rotamer, bond-angle, tetrahedral-geometry, or peptide-planarity assumption. A tighter conversion is possible only by adding structural assumptions and should be treated as a different model.

---

## 7. Atom sets and r^-6 averaging

## 7.1 Unnormalized sum convention

Suppose one source alternative expands to `N` explicit nucleus pairs with distances `r1...rN`, and the effective distance is

\[
r_{\mathrm{eff}} = \left(\sum_{i=1}^{N} r_i^{-6}\right)^{-1/6}.
\]

If `r_eff <= U`, then

\[
\sum_i r_i^{-6} \ge U^{-6}.
\]

Let `r_min = min_i r_i`. Because every `r_i >= r_min`,

\[
\sum_i r_i^{-6} \le N r_{\min}^{-6}.
\]

Combining inequalities gives

\[
N r_{\min}^{-6} \ge U^{-6}
\]

and therefore

\[
\boxed{r_{\min} \le N^{1/6}U.}
\]

Thus at least one explicit pair satisfies an upper bound `alpha*U`, where

\[
\alpha=N^{1/6}.
\]

The heavy-parent result is a disjunction over explicit pairs:

\[
\bigvee_i \left[D(P_{1i},P_{2i}) \le N^{1/6}U + L_{1i}+L_{2i}+m\right].
\]

## 7.2 Normalized mean convention

For

\[
r_{\mathrm{eff}} = \left(\frac{1}{N}\sum_i r_i^{-6}\right)^{-1/6},
\]

the same argument gives

\[
r_{\min}\le U.
\]

The `mean-r6` policy therefore uses factor 1.

## 7.3 Hard existential assignment

If a restraint means simply that at least one assignment is correct and that assignment obeys `r <= U`, the `hard-or` policy also uses factor 1.

## 7.4 Calibration caveat

Deposited upper bounds may already contain pseudoatom or multiplicity corrections imposed by the structure-calculation program. Applying `N^(1/6)` again can make the projected contact unnecessarily weak. Conversely, omitting it for a true unnormalized atom-set sum can be non-conservative. This is a dataset-level scientific choice and should be confirmed from deposition metadata, the original structure-calculation protocol, or the authors.

---

## 8. Wildcards, stereochemistry, and pseudoatoms

## 8.1 Wildcard atom sets

Examples:

- `VAL HG1%` expands to `HG11/HG12/HG13`; all share parent `CG1` and can collapse after projection.
- `ALA HB%` expands to three protons sharing `CB`.
- `TYR HD%` expands to `HD1/HD2`, whose parents are `CD1/CD2`; this remains a two-parent OR.

The factor `N` is the number of explicit cross-products between endpoint atom sets. A 3-proton methyl to another 3-proton methyl gives `N=9` under the sum-r6 interpretation.

## 8.2 x/y non-stereospecific identifiers

`x` and `y` are treated as assignment branches, not as an r^-6 set by themselves. For example, `VAL HGx%` generates a branch for the `CG1` methyl and a branch for the `CG2` methyl. The converter does not solve global x/y consistency across all restraints. That requires a dataset-level assignment model.

## 8.3 Geometric pseudoatoms

Names such as `QB` or `MG1` may represent geometric pseudoatom positions with convention-specific correction values. A geometric pseudoatom is not a real nucleus and is not mathematically identical to the wildcard set of constituent protons. The default is therefore rejection.

`--pseudoatom-policy atomset` replaces a simple Q/M name with an approximate H wildcard and reports the approximation. This mode should not be used without expert confirmation of the originating program's pseudoatom definition and any correction already applied to the bound.

---

## 9. Boltz compatibility and bound handling

The target BoltzUI atom-contact implementation, validated at revision `c3e5c7f6ae80d9261c357a0951a0929c50a1115d`, uses exact selectors:

```yaml
atom1: [CHAIN_ID, RESIDUE_INDEX, ATOM_NAME]
atom2: [CHAIN_ID, RESIDUE_INDEX, ATOM_NAME]
max_distance: 6.5
force: true
```

The accepted interval in that custom fork is 2-20 A. Each atom pair retains its exact atom-level distance potential. When several exact atom restraints share one token pair, token-level contact conditioning uses the minimum threshold while the exact atom potentials remain separate. This schema is specific to the audited BoltzUI fork and must not be assumed interchangeable with another Boltz or BoltzUI revision.

### 9.1 Lower than 2 A

Raising a bound from, for example, 1.8 A to 2.0 A weakens it. This is logically safe and is recorded as a Boltz adjustment.

### 9.2 Greater than 20 A

Reducing 21 A to 20 A would strengthen the constraint. The converter never does this. The value is retained in the audit report and rejected from normal Boltz YAML.

### 9.3 Precision

Calculations use finite Python floating-point values. Executable YAML and tabular upper bounds are rounded upward to six decimal places. Upward rounding can weaken an upper bound by less than one millionth of an angstrom but can never tighten it. Full unrounded formula inputs remain in JSON provenance.

---

## 10. Ambiguous groups and recommended execution strategies

## 10.1 Safe default

Only groups that collapse to one heavy pair are emitted in
`atom_constraints_exact.yaml`. This is the recommended first-stage guidance
set. Multi-pair groups are emitted separately in
`atom_constraints_union.yaml`; exact contacts never appear there merely as a
way to share a file.

## 10.2 Union-aware potential

A correct representation of a multi-pair ambiguity is a group of alternatives
sharing one union identifier in the distance potential. NMR2Boltz writes that
minimal schema to `atom_constraints_union.yaml`; the corresponding consumer
implementation plan is described in `BOLTZUI_UNION_EXTENSION.md`.

A critical restriction is that token-level conditioning must not mark every alternative as an independent contact. Doing so communicates an AND relation to the model even if the coordinate potential is union-aware.

## 10.3 Assignment hypotheses

Until union-aware parsing is available, `--hypotheses N` writes multiple YAML files. Each file chooses one alternative from each ambiguous group and combines those choices with all safe contacts. This is a model-selection strategy, not a conservative logical conversion.

Recommended use:

1. generate several hypotheses and multiple diffusion samples per hypothesis;
2. include an unrestrained control;
3. post-score every result against the original disjunctive restraints;
4. retain structures that satisfy the largest consistent restraint subset without severe stereochemical problems.

The number of assignment combinations grows exponentially, so hypotheses are sampled deterministically from a user seed.

---

## 11. Validation protocol

A production validation should include at least four levels.

### 11.1 Parsing validation

- compare group and row counts with the source file;
- confirm the detected format and restraint origin;
- inspect missing/null upper bounds;
- verify that canonical expansion rows were reconstructed only when their
  complete author/source semantics, topology-proven proton sets, and Cartesian
  product match;
- review all non-null combination identifiers.
- review every NMR-STAR member logic code and confirm that only explicit OR groups were flattened.

### 11.2 Identifier validation

- compare `sequence_map.tsv` with the exact Boltz YAML sequence;
- spot-check insertion codes, chain aliases, and modified residues;
- verify that every emitted atom exists in the Boltz residue template;
- reject any mapping derived only from restraint order unless manually checked.

### 11.3 Geometric validation of candidate structures

- calculate every emitted heavy-heavy distance;
- reconstruct hydrogens using chemically appropriate protonation and geometry;
- calculate the original H-H/H-X distances;
- evaluate the original atom-set averaging convention;
- report per-restraint and per-residue violations, not only a global score.

### 11.4 Structural and statistical controls

- compare restrained and unrestrained Boltz ensembles;
- inspect whether one restraint dominates a fold transition;
- assess backbone/side-chain stereochemistry and clashes;
- test sensitivity to projection margin and averaging policy;
- cross-validate by withholding a subset of restraints;
- where available, validate against RDCs, PREs, J couplings, SAXS, cryo-EM density, or independent biochemical data.

### 11.5 Paired-format discrepancy audit and CI gate

The corpus validator treats format parity as row-level scientific evidence, not
as a count comparison. Every NEF-only, NMR-STAR-only, or common heavy pair with
a different final bound records the source restraint/group IDs and row IDs,
author atom expressions, canonical atom expansions, resolved physical proton
sets, pseudoatom policy outcome, explicit proton-pair count `N`, averaging
policy, deposited upper bound, projected heavy pair, projected terms, and final
bound. Corresponding source groups are joined by normalized list name plus
restraint ID, including groups rejected on one side.

Each row is classified as:

- `expected_format_difference` only when a tested predicate proves wildcard
  atom set versus explicit OR members, x/y assignment versus a compatible
  physical set, rejected geometric Q/M pseudoatom handling, or a verified
  canonical naming alias with identical physical atoms and destinations;
- `deposition_inconsistency` when a corresponding restraint is absent, its bound
  differs, or its sequence/residue identifiers resolve inconsistently;
- `parser_projection_bug` only after a reproducible implementation defect is
  verified; or
- `unresolved` whenever the available evidence is insufficient.

The normal corpus command fails closed if any implication fails, any audit row
is unresolved or marks a parser/projection bug, the reviewed audit digest or
metric snapshot changes, or the exact reviewed missing-coordinate set changes.
Known corpus limitations are pinned by contact identity, bound, provenance, and
digest rather than waived by count. A run still writes its evidence before
returning a nonzero exit. This prevents denominator omission and baseline drift
from being mistaken for robustness.

### 11.5 Executed validation record

The following checks were executed on 2026-07-18 against the current source tree:

- all 91 Pytest regression, format, topology, logic, target-validation,
  ensemble-alignment, constraint-serialization, and robustness tests passed;
- Python byte compilation passed for source, tests, and the stress harness;
- 100,000 randomized sum-r6 implication cases and 100,000 constructive triangle-inequality cases passed in the final Docker image;
- 25,000 outward-rounding cases and 10,000 randomized OR-max/AND-min order-invariance cases passed;
- all 845 explicit hydrogens in the built-in protein and nucleic-acid topologies resolved to one finite, positive parent-bond envelope;
- NEF, NMR-STAR, compressed NEF, embedded custom-component topology, all three averaging policies, and 32 deterministic ambiguity hypotheses were exercised;
- strict-mode and non-finite-input failure paths returned their documented exit codes, 3 and 2;
- sequence-only NEF/NMR-STAR files produce empty, auditable distance conversions;
- source residue identities are checked against the resolved sequence record before topology lookup;
- every executable atom is proven against its mapped component both before
  deduplication and again before output serialization;
- every conversion writes a polymer-only FASTA sequence file;
- PDB author numbering is aligned to Boltz one-based sequence positions before coordinate evaluation.

External semantics were checked at fixed revisions: BoltzUI `c3e5c7f6ae80d9261c357a0951a0929c50a1115d`, NEF `9ab6bc023a406c87df407837597efecaf289fe55`, and the NMR-STAR dictionary `35c6e32a4c948de8d9bf8b367dfab1217a216c0a`.

A paired-format benchmark was run for 12 deposited NMR structures. All 24 NEF
and NMR-STAR conversions completed, including two valid empty distance
conversions for 8S8O. Conservative defaults emitted 12,998 NEF and 11,829
NMR-STAR contacts. Resolved contact/model satisfaction against the deposited
ensembles was 99.88% and 99.86%, respectively. The projected implication had
zero failures in 379,449 cases with satisfied source antecedents. The 4,177-row
format audit contains 4,134 allowlisted expected differences, 43 deposition
inconsistencies, zero unresolved rows, and zero remaining parser/projection
bugs. Exact
pair-and-bound parity was observed for three positive-distance cases; the
remaining discrepancies are retained as explicit audit evidence rather than
being silently approximated.

No GPU Boltz structure-prediction campaign was run. Accordingly, these results
support robustness and generality of parsing, logical projection, topology
resolution, serialization, and coordinate auditing, but do not establish
improvement over unrestrained Boltz. The PDF was intentionally not regenerated.

---

## 12. Suggested confidence tiers

A useful production filter can classify projected information as follows:

| Tier | Description | Recommended action |
|---|---|---|
| A | Explicit atoms, explicit upper bound, verified sequence map, standard/CCD topology, one heavy pair | Emit directly |
| B | Wildcard/atom set but all explicit atoms share one heavy pair; averaging policy confirmed | Emit with audit |
| C | Several heavy-pair assignment alternatives | Union potential or hypothesis ensemble |
| D | Derived upper bound, inferred residue map, pseudoatom approximation, unknown topology, or complex logic | Expert review; usually do not use for initial folding |

The current software records the components needed to implement such filtering, although it does not assign a single opaque confidence score.

---

## 13. Potential methodological extensions

### 13.1 Probabilistic union restraints

Use peak-assignment probabilities or restraint weights to define a weighted soft minimum across alternatives rather than a uniform union. This requires careful calibration because source `weight` fields do not always encode assignment probabilities.

### 13.2 Ensemble-aware objective

NMR observables often report ensemble/time averages. A future multi-sample Boltz objective could apply r^-6 averaging across generated conformers rather than forcing each individual conformer to satisfy every contact.

### 13.3 Iterative assignment and folding

An EM-like workflow could alternate between:

1. predicting structures under currently plausible assignments;
2. reweighting ambiguous assignments from structural compatibility;
3. regenerating structures with updated union weights.

This can be powerful but risks confirmation bias and should use held-out restraints.

### 13.4 Information-aware restraint selection

Redundant local restraints can overwhelm a smaller number of informative long-range contacts. Selection can consider sequence separation, graph coverage, ambiguity, upper-bound tightness, and clustering. Any down-selection should remain visible and reproducible.

### 13.5 Broader experimental data

The same architecture can host other guidance types, but their physics differs:

- PRE: broad r^-6-sensitive upper information with paramagnetic-center uncertainty;
- FRET: distribution- and dye-linker-dependent distances;
- crosslinks: chemistry-specific accessible-volume bounds;
- RDC: orientational, not distance, information;
- hydrogen bonds: donor/acceptor and angular geometry.

They should not be forced into the proton-parent formula without a modality-specific model.

---

## 14. Reproducibility and auditability

Every run writes:

- exact command settings;
- sequence mapping and source;
- source row IDs and author endpoint expressions;
- complete normalized source restraint groups, source observations, upper/lower/target/uncertainty values, weight, origin, member logic, and selected fallback policy;
- atom-set size and averaging factor;
- X-H offset and final projected value;
- ambiguity and rejection reasons;
- Boltz min/max adjustments.

The JSON report is intended to be the primary reproducibility object. The compact YAML is only the execution artifact.

---

## 15. Questions for NMR expert review before production deployment

1. Were the upper bounds calibrated as explicit pair bounds, normalized r^-6 means, unnormalized r^-6 sums, or program-specific pseudoatom-corrected values?
2. Do wildcard names denote physical atom sets, unresolved assignments, or exported pseudoatoms in this dataset?
3. Are x/y labels globally correlated across the restraint list?
4. Are repeated restraint-ID rows true OR alternatives, and are any nested combination IDs used?
5. Which sequence namespace matches the Boltz input exactly?
6. Are modified residues and protonation/tautomer states represented correctly by the selected CCD topology?
7. Should exchangeable-proton restraints be retained under the experimental conditions?
8. What violation tolerance and ensemble-averaging convention was used in the original NMR structure calculation?
9. Which restraints should be withheld for cross-validation?
10. Is the scientific goal fold rescue, local refinement, ranking, or production of a restraint-satisfying NMR ensemble? The optimal potential strength and selection differ among these goals.

---

## 16. Primary specifications and references

1. NMR Exchange Format specification and commented examples, NMRExchangeFormat/NEF, audited revision `9ab6bc023a406c87df407837597efecaf289fe55`: https://github.com/NMRExchangeFormat/NEF
2. Gutmanas A. et al. *NMR Exchange Format: a unified and open standard for representation of NMR restraint data.* Nature Structural & Molecular Biology 22, 433-434 (2015). DOI: 10.1038/nsmb.3041
3. Biological Magnetic Resonance Data Bank: https://bmrb.io/
4. PyNMRSTAR, BMRB-maintained STAR parser: https://github.com/uwbmrb/PyNMRSTAR
5. wwPDB file-format and NMR data-exchange documentation: https://www.wwpdb.org/documentation/file-format-content
6. wwPDB Chemical Component Dictionary: https://www.wwpdb.org/data/ccd
7. Boltz source repository and potential implementation: https://github.com/jwohlwend/boltz
8. BoltzUI target atom-contact implementation, audited revision `c3e5c7f6ae80d9261c357a0951a0929c50a1115d`: https://github.com/Jalil-Mahdizadeh/BoltzUI
9. BMRB NMR-STAR dictionary, audited revision `35c6e32a4c948de8d9bf8b367dfab1217a216c0a`: https://github.com/bmrb-io/nmr-star-dictionary

---

## 17. Concise algorithm

```text
parse STAR syntax with PyNMRSTAR
identify NEF or NMR-STAR distance loop by category
build an explicit source-to-Boltz residue map
normalize rows into restraint-ID OR groups
prefer author atom expressions; retain canonical atom hints
select only an explicit upper bound unless a fallback policy was requested
resolve each endpoint expression into one or more atom-set/assignment branches
map every proton to its directly bonded heavy parent from topology
for each atom-set branch pair:
    N = number of explicit nucleus-pair combinations
    alpha = N^(1/6) for sum-r6, else 1
    for every explicit pair:
        U_parent = alpha * U_source + L1 + L2 + margin
        group by resulting heavy-parent pair
merge duplicate heavy pairs within one OR group using max(U_parent)
if one heavy pair remains:
    add it to the safe set
else:
    retain the OR group; do not emit every alternative
merge the same heavy pair across independent groups using min(U_parent)
raise values below Boltz minimum; never clip values above Boltz maximum
write safe YAML plus complete audit, ambiguity, rejection, and mapping files
```
