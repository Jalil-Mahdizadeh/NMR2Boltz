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

The exact-contact/PDB-model distance matrix is generated after the corpus gate:

```bash
python validation/distance_check.py benchmark/input \
  --conversion-output benchmark/output \
  --output-directory benchmark/distance_check
```

This creates one CSV per PDB ID under `distance_check/`. Each row is the union
of the NEF- and NMR-STAR-emitted exact heavy-atom pairs, followed by the
corresponding six-decimal bounds and one Euclidean distance column per PDB
model. A blank bound means that format did not emit the pair; a blank distance
means that the deposited model does not observe an endpoint. Ambiguous union
groups are excluded because treating their alternatives as independent
constraints would change OR semantics. See
[`distance_check/README.md`](distance_check/README.md) for current counts,
scientific limitations, and artifact digests.

The runner uses conservative defaults (`sum-r6`, explicit upper bounds only,
and pseudoatom rejection) and evaluates every deposited conformer with a 0.01 Å
coordinate tolerance. The PDB sequence is aligned to Boltz one-based positions,
so author numbering, gaps, and partial coordinate observations do not silently
mis-index the audit.

Exact format parity means both inputs emitted the same heavy-atom pairs with the
same bounds within 1e-6 Å. High PDB satisfaction does not replace parity review
and is not an independent Boltz prediction-accuracy result.

The command is fail-closed: execution errors, implication failures, unresolved
discrepancies, parser/projection bugs, discrepancy-digest changes, scientific
metric changes, and changes to the exact reviewed missing-coordinate set all
produce a nonzero exit. The reviewed set contains full contact identities,
bounds, source provenance, and a SHA-256 digest rather than a count-only waiver.
Audit files are still written on failure so the cause remains inspectable.
Baseline replacement requires the explicit `--write-reviewed-baseline` flag
and should only follow review of every changed row and metric.
