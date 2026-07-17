# Paired NEF/NMR-STAR benchmark

`input/` contains one directory per PDB ID. Each case must contain exactly one
`.pdb`, one `.nef`, and one `.str` file. Run the complete coordinate audit from
the repository root:

```bash
python validation/benchmark_corpus.py benchmark/input \
  --output-directory benchmark/output \
  --reviewed-baseline benchmark/reviewed_baseline.json
```

Results are written to:

```text
benchmark/output/
  <PDB-ID>/
    nef/
      conversion_report.json
      coordinate_summary.json
      sequences.fasta
      ...
    star/
      conversion_report.json
      coordinate_summary.json
      sequences.fasta
      ...
    format_parity.json
    format_discrepancy_audit.tsv
    format_discrepancy_summary.json
  FORMAT_DISCREPANCY_AUDIT.tsv
  format_discrepancy_summary.json
  benchmark_summary.json
  BENCHMARK_REPORT.md
  RUN_COMMAND.txt
```

The runner uses conservative defaults (`sum-r6`, explicit upper bounds only,
and pseudoatom rejection) and evaluates every deposited conformer with a 0.01 Å
coordinate tolerance. The PDB sequence is aligned to Boltz one-based positions,
so author numbering, gaps, and partial coordinate observations do not silently
mis-index the audit.

Exact format parity means both inputs emitted the same heavy-atom pairs with the
same bounds within 1e-6 Å. High PDB satisfaction does not replace parity review
and is not an independent Boltz prediction-accuracy result.

The command is fail-closed: execution errors, implication failures, unresolved
discrepancies, missing heavy-atom coordinate resolution, and changes from
`reviewed_baseline.json` all produce a nonzero exit. The audit files are still
written on a scientific gate failure so the cause remains inspectable. Baseline
replacement requires the explicit `--write-reviewed-baseline` flag and should
only follow review of every changed row and metric.
