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

Run Boltz 2 with inference potentials enabled. Keep the converter's JSON audit beside every prediction so that chain and residue mapping, atom expansion, bound calculation, merge decisions, and omissions remain reviewable.

## Ambiguous restraints

`atom_constraints_union.yaml` contains only `atom_contact_union` groups. Each
alternative retains its own upper bound, and the file deliberately contains no
provenance metadata. A consumer must implement the alternatives as one OR
potential; see `BOLTZUI_UNION_EXTENSION.md`. The complete source provenance
remains in `conversion_report.json` and `ambiguous_groups.tsv`.

`hypotheses/`, when requested, provides one or more ordinary current-schema YAML files. Each file chooses one branch per ambiguous group. These are alternative calculations, not restraints to combine in one run. Rank or filter the resulting structures against the original proton-level restraints after adding hydrogens.

## Residue numbering

Boltz atom selectors use one-based positions in the sequence supplied to Boltz. NEF `sequence_code` is a string and may contain insertion-like suffixes. Supply an explicit mapping table whenever the NMR numbering is not a simple one-to-one sequence order. The converter records every mapping decision in the audit output.

`sequences.fasta` is a clean polymer-only view of the resolved sequence map.
Caps, ions, and non-polymer ligands are intentionally omitted from FASTA and
must be represented separately when they are part of the Boltz target.

## Recommended run protocol

1. Convert with strict defaults and inspect all rejected and ambiguous groups.
2. Confirm sequence mapping and non-standard component topology. Every emitted
   atom is checked against the mapped component dictionary; do not substitute a
   coordinate-presence check for this topology invariant.
3. Ask an NMR expert to confirm the atom-set averaging convention used to calibrate each restraint list.
4. Run multiple Boltz diffusion samples per safe or hypothesis input.
5. Add hydrogens to generated structures using an appropriate protonation model.
6. Re-evaluate the original NMR restraint logic at the proton or atom-set level; do not validate only the projected heavy-atom inequalities.
7. Compare unconstrained and guided ensembles to detect over-steering.
