# Validation status

Status: PASS
Date: 2026-07-16

## Regression and stress validation

- 31 Pytest tests passed.
- 100,000 randomized sum-r6 implication cases passed.
- 100,000 constructive triangle-inequality cases passed.
- 25,000 outward-rounding cases passed.
- 10,000 randomized OR-max/AND-min merge-order cases passed.
- All 850 built-in hydrogen-parent mappings resolved.
- Deterministic NEF, NMR-STAR, compressed-input, custom-component, averaging-policy, and 32-hypothesis paths passed.

## Real deposited-data validation

PDB 6M6O and linked BMRB 28061 were tested across all 10 deposited conformers.
The PDB and mmCIF files agreed exactly for 30,820 atom instances. The NOE-only
conversion emitted 2,707 heavy-atom contacts. In 25,350 cases where all source
antecedents were satisfied, the projected contact implication had zero failures.

The detailed report, commands, raw data, CSV evidence, and JSON summary are under
workspace/output. This benchmark validates conversion behavior, not the predictive
accuracy of a GPU Boltz folding campaign.

## Container

- Image: nmr2boltz:0.1.0-validated
- Digest: sha256:ee0bb2c1c79cf569fe274e363add28f4d39ec0f3d6201d22c68a0576c05fa74b
- Default user: 65532:65532
- Runtime checks used no network and a read-only root filesystem.
- Only this nmr2boltz image remains.
