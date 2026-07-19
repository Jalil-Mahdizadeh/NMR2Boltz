# NMR expert review checklist

Use this checklist before treating converted contacts as production restraints.

## Source and calibration

- [ ] Confirm which restraint lists are NOE-derived and which are hydrogen bonds, PRE, crosslink, or other distance-like data.
- [ ] Confirm whether upper bounds are explicit, target-derived, or calibrated from intensity classes.
- [ ] Determine whether atom-set distances use unnormalized r^-6 sums, normalized r^-6 means, hard assignment alternatives, center averaging, or software-specific handling.
- [ ] Determine whether pseudoatom/multiplicity corrections are already included in deposited upper bounds.
- [ ] Review any non-null `restraint_combination_id`/`Combination_ID` and nested logic.

## Identifier mapping

- [ ] Compare every chain in `sequence_map.tsv` with the exact Boltz YAML.
- [ ] Compare `sequences.fasta` with the intended polymer sequence; confirm that omitted caps, ions, and ligands are represented separately in the Boltz target.
- [ ] Resolve every `sequence_residue_mismatch`; do not repair it implicitly through atom topology.
- [ ] Verify insertion codes, residue-number gaps, negative numbering, and chain aliases.
- [ ] Verify all modified residues and ligands against the selected CCD chemistry.
- [ ] Ensure the Boltz residue index is the one-based position in the input sequence, not the author residue number.
- [ ] Reject or manually verify every mapping labeled inferred.
- [ ] For `--exclude-intrachain`, verify that mapped Boltz chain IDs represent
  the intended protein/DNA/RNA chains; source author-chain labels are not the
  filtering boundary.

## Atom semantics

- [ ] Confirm that `%`/`*` expressions are intended atom sets rather than exporter-specific pseudoatom notation.
- [ ] Confirm whether x/y labels require global stereospecific consistency.
- [ ] Review every Q/M pseudoatom; do not use `--pseudoatom-policy atomset` without a program-specific correction rationale.
- [ ] Review exchangeable N/O/S/Se protons under the experimental pH, solvent, temperature, and tautomer/protonation state.
- [ ] Spot-check direct proton-parent mappings for each residue class present.
- [ ] Review every `atom_not_present_in_mapped_residue` quarantine against the
  source row and mapped component; never repair it by guessing an atom or residue.

## Projection

- [ ] Approve the X-H bond-length upper envelopes and any additional projection margin.
- [ ] Approve the selected `sum-r6`, `mean-r6`, or `hard-or` policy.
- [ ] Inspect every projected bound near or above the 20 Å Boltz limit.
- [ ] Confirm that no over-20 Å value was clipped.
- [ ] Confirm that multiple alternatives in one restraint ID were not emitted as simultaneous contacts.
- [ ] Confirm that mixed intra/inter OR groups were omitted in full rather than
  reduced to only their inter-chain alternatives.
- [ ] When paired NEF/NMR-STAR files exist, inspect exact atom-pair and bound parity rather than only total counts.

## Boltz usage

- [ ] Use Boltz-2 with potentials enabled and `force: true` as required by the atom-contact implementation.
- [ ] Start with the safe unambiguous set; retain an unrestrained control.
- [ ] Generate multiple diffusion samples and, for ambiguous data, multiple assignment hypotheses or union-aware runs.
- [ ] Avoid flooding the model with redundant local restraints; document any selection criteria.
- [ ] Record the exact Boltz/BoltzUI commit, model, seeds, sampling settings, and potential settings.

## Post-prediction validation

- [ ] Reconstruct hydrogens with an NMR-aware method and appropriate protonation/tautomer states.
- [ ] Evaluate the original source restraints—not only the projected heavy-atom contacts.
- [ ] Apply the original atom-set/pseudoatom averaging convention.
- [ ] Report individual violations, RMS violation, largest violation, and fraction satisfied.
- [ ] Check covalent geometry, chirality, clashes, Ramachandran/rotamer quality, and chain topology.
- [ ] Compare restrained and unrestrained ensembles.
- [ ] Withhold a restraint subset for cross-validation where data volume permits.
- [ ] Validate against independent observables when available.

## Sign-off

- NMR expert:
- Structural-modeling expert:
- Dataset/entry identifiers:
- Converter version/commit:
- BoltzUI/Boltz commit:
- Approved policies:
- Known limitations accepted:
