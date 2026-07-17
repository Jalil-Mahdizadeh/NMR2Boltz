# BoltzUI integration guide

## Current patched schema

The converter's `boltz_constraints.yaml` contains only independent, unambiguous heavy-atom contacts that the present BoltzUI patch can represent:

```yaml
constraints:
  - atom_contact:
      atom1: [A, 12, CB]
      atom2: [A, 44, CG1]
      max_distance: 6.72
      force: true
```

Run Boltz 2 with inference potentials enabled. Keep the converter's JSON audit beside every prediction so that chain and residue mapping, atom expansion, bound calculation, merge decisions, and omissions remain reviewable.

## Ambiguous restraints

`proposed_atom_contact_unions.yaml` preserves alternatives under a shared OR group. It is not accepted by the current patch; it is the recommended input to the extension described in `BOLTZUI_UNION_EXTENSION.md`.

`hypotheses/`, when requested, provides one or more ordinary current-schema YAML files. Each file chooses one branch per ambiguous group. These are alternative calculations, not restraints to combine in one run. Rank or filter the resulting structures against the original proton-level restraints after adding hydrogens.

## Residue numbering

Boltz atom selectors use one-based positions in the sequence supplied to Boltz. NEF `sequence_code` is a string and may contain insertion-like suffixes. Supply an explicit mapping table whenever the NMR numbering is not a simple one-to-one sequence order. The converter records every mapping decision in the audit output.

`sequences.fasta` is a clean polymer-only view of the resolved sequence map.
Caps, ions, and non-polymer ligands are intentionally omitted from FASTA and
must be represented separately when they are part of the Boltz target.

## Recommended run protocol

1. Convert with strict defaults and inspect all rejected and ambiguous groups.
2. Confirm sequence mapping and non-standard component topology.
3. Ask an NMR expert to confirm the atom-set averaging convention used to calibrate each restraint list.
4. Run multiple Boltz diffusion samples per safe or hypothesis input.
5. Add hydrogens to generated structures using an appropriate protonation model.
6. Re-evaluate the original NMR restraint logic at the proton or atom-set level; do not validate only the projected heavy-atom inequalities.
7. Compare unconstrained and guided ensembles to detect over-steering.
