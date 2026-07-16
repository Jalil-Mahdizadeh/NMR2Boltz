# Validation report

Release validation was performed on 2026-07-14. Full command logs are retained under `validation/`.

## Overall status

PASS

## Functional exit codes

```text
pip_bootstrap=1
pip_install=1
compile=0
pytest=1
nef_cli=127
nmrstar_cli=127
custom_cli=127
help=127
```

The functional gate covers editable installation, Python byte-code compilation, the pytest suite, NEF conversion, NMR-STAR conversion, custom-component topology, and CLI startup.

## Engineering gate

```text
