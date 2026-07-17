# Validation status

Status: PASS WITH DOCUMENTED FORMAT-PARITY LIMITATIONS
Date: 2026-07-17

## Regression and stress validation

- 44 Pytest tests passed.
- 100,000 randomized sum-r6 implication cases passed.
- 100,000 constructive triangle-inequality cases passed.
- 25,000 outward-rounding cases passed.
- 10,000 randomized OR-max/AND-min merge-order cases passed.
- All 850 built-in hydrogen-parent mappings resolved.
- Deterministic NEF, NMR-STAR, compressed-input, custom-component, averaging-policy, and 32-hypothesis paths passed.
- Target-YAML validation checks chain existence, one-based index bounds, canonical residue identity, declared modifications, mapping collisions, and emitted-contact positions.
- Sequence-only NEF/NMR-STAR entries produce auditable empty distance conversions.
- Sequence/residue conflicts are rejected explicitly before atom-topology resolution.
- PDB coordinates are aligned to one-based Boltz sequence indices before distance evaluation.
- NMR-STAR sequence-alias provenance is sorted deterministically rather than depending on Python set iteration.
- Sequence-map ordering remains deterministic across Python hash seeds.

## Paired-format deposited-data validation

Twelve deposited NMR ensembles were converted from both NEF and NMR-STAR.
All 24 conversions completed; the two 8S8O inputs are valid empty distance
conversions containing sequence, shifts, and torsion data but no distance loop.
The default safe projection emitted 13,004 NEF and 11,841 NMR-STAR contacts.
Resolved contact/model satisfaction was 99.88% and 99.86%, respectively.
Across 390,810 model/contact cases with satisfied source antecedents, the
projected heavy-atom implication had zero failures.

Only three positive-distance cases have exact NEF/STAR pair-and-bound parity
(21CC, 9D99, and 9KG4); 8S8O also has exact empty-output parity. Large differences
caused by pseudoatoms, atom-set multiplicity, and source inconsistencies remain
visible and are not silently approximated. The complete current report is under
`benchmark/output`. This benchmark validates conversion behavior, not the
predictive accuracy of a GPU Boltz folding campaign.

## Container status

The previously validated container predates the current sequence-validation,
FASTA-output, and corpus-runner changes. Rebuild and revalidate the image before
using a container digest as the release artifact.
