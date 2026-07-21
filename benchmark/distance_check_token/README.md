# Token-contact distance check

Each `<PDB-ID>.csv` contains the union of the ordinary token contacts
emitted from that entry's NEF and NMR-STAR conversions. NEF and STAR
bounds remain separate; a blank bound means that format did not emit the
residue pair. Contacts may originate from exact atom restraints or from
safely collapsed single-token-pair union restraints.

For every deposited PDB model, the reported geometric value is the minimum
Euclidean distance between any non-hydrogen atom in the first residue and
any non-hydrogen atom in the second residue. A blank model cell means that
one or both residues have no observed heavy atom in that model.

PDB author numbering is aligned independently through each format's generated
sequence map to Boltz one-based residue indices. For common pairs, generation
fails if NEF- and STAR-driven distances differ by more than 1e-06 A
or if only one format resolves the coordinate. Token YAML is cross-checked
against conversion-report provenance, including canonical unique pairs,
4-20 A bounds, conservative six-decimal outward rounding, and `force: false`.
Geometric PDB distances are reported to six decimals.

The satisfaction columns below are descriptive comparisons using
`minimum heavy-atom distance <= max_distance + 0.01 A`;
they are not model-prediction accuracy estimates.

| PDB | models | rows | common | NEF only | STAR only | different bounds | missing distances | NEF satisfied | STAR satisfied |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 21CC | 10 | 18 | 18 | 0 | 0 | 0 | 0 | 100.00% | 100.00% |
| 43JX | 20 | 442 | 442 | 0 | 0 | 59 | 0 | 99.97% | 99.97% |
| 6M6O | 10 | 828 | 828 | 0 | 0 | 21 | 0 | 100.00% | 100.00% |
| 8R1X | 20 | 286 | 286 | 0 | 0 | 46 | 220 | 99.64% | 99.64% |
| 8S8O | 20 | 0 | 0 | 0 | 0 | 0 | 0 | n/a | n/a |
| 9CCH | 20 | 678 | 678 | 0 | 0 | 271 | 0 | 99.98% | 100.00% |
| 9D99 | 10 | 146 | 146 | 0 | 0 | 0 | 0 | 100.00% | 100.00% |
| 9KG4 | 15 | 53 | 53 | 0 | 0 | 0 | 0 | 100.00% | 100.00% |
| 9PQH | 10 | 33 | 2 | 1 | 30 | 2 | 0 | 33.33% | 100.00% |
| 9SGX | 20 | 399 | 285 | 114 | 0 | 72 | 0 | 100.00% | 100.00% |
| 9VQ1 | 20 | 54 | 54 | 0 | 0 | 0 | 0 | 100.00% | 100.00% |
| 9VUY | 20 | 440 | 287 | 153 | 0 | 82 | 0 | 99.77% | 100.00% |

Reproduce from the repository root after running the corpus benchmark:

```bash
mkdir -p benchmark/distance_check_token
docker run --rm --network none --read-only \
  -v "$PWD:/work:ro" \
  -v "$PWD/benchmark/distance_check_token:/output" \
  --entrypoint python -w /work nmr2boltz:latest \
  validation/distance_check_token.py benchmark/input \
  --conversion-output benchmark/output --output-directory /output
```

`distance_check_token_summary.json` records model IDs, row counts,
missing-coordinate counts, descriptive satisfaction totals, and a
SHA-256 digest for every CSV.
