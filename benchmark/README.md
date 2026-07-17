# Paired NEF/NMR-STAR benchmark

`input/` contains one directory per PDB ID. Each case must contain exactly one
`.pdb`, one `.nef`, and one `.str` file. Run the complete coordinate audit from
the repository root:

```bash
python validation/benchmark_corpus.py benchmark/input \
  --output-directory benchmark/output
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
