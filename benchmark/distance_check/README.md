# Exact-contact distance check

Each `<PDB-ID>.csv` contains the union of the topology-validated exact
heavy-atom pairs emitted from that entry's NEF and NMR-STAR inputs.
A blank format-bound cell means that format did not emit the pair. A blank
model-distance cell means that at least one endpoint is not observed in that
deposited PDB model. Ambiguous OR/union groups are intentionally excluded:
flattening their alternatives into independent restraints would change their
scientific meaning.

PDB author numbering is aligned independently through each format's generated
sequence map to Boltz one-based residue indices. For common pairs, generation
fails if NEF- and STAR-driven aligned distances differ by more than 1e-06 A
or if only one format resolves the coordinate. Executable YAML bounds are
cross-checked against conversion provenance and retain conservative six-decimal
outward rounding. Geometric PDB distances are reported to six decimals.

The satisfaction columns below are descriptive comparisons using
`distance <= max_distance + 0.01 A`; they are not model-prediction
accuracy estimates.

| PDB | models | rows | common | NEF only | STAR only | different bounds | missing distances | NEF satisfied | STAR satisfied |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 21CC | 10 | 176 | 176 | 0 | 0 | 0 | 0 | 99.72% | 99.72% |
| 43JX | 20 | 1716 | 1716 | 0 | 0 | 353 | 0 | 99.96% | 99.96% |
| 6M6O | 10 | 2878 | 2818 | 1 | 59 | 89 | 0 | 99.72% | 99.60% |
| 8R1X | 20 | 1114 | 1114 | 0 | 0 | 160 | 820 | 99.74% | 99.74% |
| 8S8O | 20 | 0 | 0 | 0 | 0 | 0 | 0 | n/a | n/a |
| 9CCH | 20 | 3516 | 3350 | 0 | 166 | 1754 | 0 | 99.98% | 100.00% |
| 9D99 | 10 | 456 | 456 | 0 | 0 | 0 | 0 | 100.00% | 100.00% |
| 9KG4 | 15 | 406 | 406 | 0 | 0 | 0 | 0 | 99.74% | 99.74% |
| 9PQH | 10 | 45 | 2 | 1 | 42 | 2 | 0 | 3.33% | 100.00% |
| 9SGX | 20 | 1578 | 765 | 813 | 0 | 0 | 0 | 100.00% | 100.00% |
| 9VQ1 | 20 | 56 | 53 | 0 | 3 | 0 | 0 | 95.00% | 94.20% |
| 9VUY | 20 | 1334 | 696 | 631 | 7 | 96 | 0 | 99.96% | 100.00% |

Reproduce from the repository root after running the corpus benchmark:

```bash
python validation/distance_check.py benchmark/input \
  --conversion-output benchmark/output \
  --output-directory benchmark/distance_check
```

`distance_check_summary.json` records model IDs, row counts, missing-coordinate
counts, descriptive satisfaction totals, and a SHA-256 digest for every CSV.
