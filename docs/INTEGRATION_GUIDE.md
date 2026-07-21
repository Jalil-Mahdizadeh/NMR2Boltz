# BoltzUI integration guide

## Exact-contact schema

The converter's `atom_constraints_exact.yaml` contains only independent,
unambiguous heavy-atom contacts:

```yaml
constraints:
- atom_contact:
    atom1: [A, 12, CB]
    atom2: [A, 44, CG1]
    max_distance: 6.720000
    force: true
```

Run Boltz 2 with `force: true`. Atom-contact guidance remains active whether
or not `--use_potentials` is selected; that option controls additional
FK/physical steering and is optional. Keep the converter's JSON audit beside
every prediction so that chain and residue mapping, atom expansion, bound
calculation, merge decisions, and omissions remain reviewable.

## Ambiguous restraints

`atom_constraints_union.yaml` contains only `atom_contact_union` groups. Each
alternative retains its own upper bound, and the file deliberately contains no
provenance metadata. A consumer must implement the alternatives as one OR
potential; see `BOLTZUI_UNION_EXTENSION.md`. The complete source provenance
remains in `conversion_report.json` and `ambiguous_groups.tsv`.

The same executable interval applies to exact and union output. A
sub-minimum union alternative is raised to the minimum, which weakens it. If
any alternative exceeds the maximum, NMR2Boltz quarantines the complete union
group rather than dropping one branch and strengthening the source OR.

## Standalone token-contact schema

Every successful conversion also writes `token_constraints.yaml` as a
coarse-grained, standalone alternative:

```yaml
constraints:
- contact:
    token1: [A, 12]
    token2: [B, 44]
    max_distance: 6.720000
    force: false
```

The projector operates on the same resolved exact constraints and OR groups as
the atom outputs; it does not parse the source again or read either generated
YAML file. An exact heavy-atom pair becomes its unordered residue-token pair.
A union is collapsed only when every alternative maps to the same non-self
token pair, using the maximum alternative bound. A union spanning several
token pairs, or containing any self-token alternative, is omitted in full so
OR is never converted into AND. Independent exact and collapsed-union
candidates for one token pair are then merged with the minimum bound.

Native token contacts require finite bounds in 4-20 Å. Sub-4 Å values are
raised to 4 Å and audited; over-20 Å semantic units are omitted rather than
clipped. This means native token contacts cannot exactly reproduce sub-4 Å
token thresholds that exact atom constraints can induce internally in the
patched BoltzUI. `force: false` deliberately avoids an additional forced
token-contact potential.

Exact atom constraints already activate token conditioning in the patched
BoltzUI, so loading `atom_constraints_exact.yaml` and `token_constraints.yaml`
together is partly redundant. Collapsed union contacts are the important
difference: they add token conditioning that `atom_constraints_union.yaml`
intentionally does not add. Treat token-only, atom-only, and hybrid runs as
distinct experimental arms and report which files were loaded.

`token_constraints.tsv` and the `token_constraints`,
`token_projection_omissions`, and `token_projection_statistics` sections of
`conversion_report.json` contain review provenance. Executable token YAML stays
minimal and metadata-free.

`hypotheses/`, when requested, provides one or more ordinary current-schema YAML files. Each file chooses one branch per ambiguous group. These are alternative calculations, not restraints to combine in one run. Rank or filter the resulting structures against the original proton-level restraints after adding hydrogens.

## Residue numbering

Boltz atom selectors use one-based positions in the sequence supplied to Boltz. NEF `sequence_code` is a string and may contain insertion-like suffixes. Supply an explicit mapping table whenever the NMR numbering is not a simple one-to-one sequence order. The converter records every mapping decision in the audit output.

`sequences.fasta` is a clean polymer-only view of the resolved sequence map.
Caps, ions, and non-polymer ligands are intentionally omitted from FASTA and
must be represented separately when they are part of the Boltz target.

## Intraresidue-excluded constraints

Use `--exclude-intraresidue` to retain only contacts whose projected heavy
atoms do not belong to the same mapped Boltz chain and residue index. The
policy applies to exact and union output. All-intraresidue unions and unions
mixing intraresidue with inter-residue alternatives are omitted in full and
recorded as `intraresidue_filtered`; trimming only the local branches would
strengthen the original OR. Confirm both the chain and residue-index mapping in
`sequence_map.tsv`.

## Inter-chain-only constraints

Use `--exclude-intrachain` when a protein, DNA, or RNA complex should be guided
only by contacts between different mapped Boltz chain IDs. The policy applies
to both exact and union outputs. Mixed intra/inter OR groups are omitted in
full, not trimmed, and appear as `intrachain_filtered` audit records. Confirm
the intended chain mapping in `sequence_map.tsv` before using the resulting
constraints.

## Recommended run protocol

1. Convert with strict defaults and inspect all rejected, ambiguous, and
   token-projection-omission records.
2. Confirm sequence mapping and non-standard component topology. Every emitted
   atom is checked against the mapped component dictionary; do not substitute a
   coordinate-presence check for this topology invariant.
   For an inter-chain-only run, also confirm that every exact endpoint pair and
   every union alternative has different mapped chain IDs.
   For an intraresidue-excluded run, confirm that no exact pair or union
   alternative has identical mapped chain and residue-index endpoints.
3. Ask an NMR expert to confirm the atom-set averaging convention used to calibrate each restraint list.
4. Select a documented token-only, atom-only, or hybrid experimental arm, then
   run multiple Boltz diffusion samples per safe or hypothesis input.
5. Add hydrogens to generated structures using an appropriate protonation model.
6. Re-evaluate the original NMR restraint logic at the proton or atom-set level; do not validate only the projected heavy-atom inequalities.
7. Compare unconstrained and guided ensembles to detect over-steering.
