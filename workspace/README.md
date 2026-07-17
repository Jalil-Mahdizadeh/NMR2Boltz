# Development validation workspace

This directory retains checksum, target-validation, and container-development
fixtures used by the versioned manifest runner.

- `input/`: authoritative deposited coordinate and NMR data plus checksums.
- `output/`: converter results, coordinate comparisons, commands, logs, and the final report.
- `benchmark.yaml`: checksum-pinned, target-aware manifest fixture.

Every `nmr2boltz` container run bind-mounts this directory at `/workspace`.

The current multi-structure scientific result is generated separately under
`benchmark/output`; see `benchmark/README.md` for the paired-format command and
output layout.
