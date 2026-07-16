# 6M6O real-data validation workspace

All downloaded inputs and generated outputs for the 6M6O validation live below this directory.

- `input/`: authoritative deposited coordinate and NMR data plus checksums.
- `output/`: converter results, coordinate comparisons, commands, logs, and the final report.

Every `nmr2boltz` container run bind-mounts this directory at `/workspace`.

Start with `output/REAL_TEST_6M6O.md`; exact commands are in
`output/RUN_COMMANDS.ps1` and machine-readable metrics are in
`output/coordinate_comparison/summary.json`.
