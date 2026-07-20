# NMR2Boltz command-line reference

This page lists every command-line flag exposed by the current
`nmr2boltz` parser. Boolean flags are disabled unless supplied. Repeatable
flags may be given more than once.

## Command structure

```text
nmr2boltz [global flags] <command> [command flags] <required arguments>
```

The required command is one of:

| Command | Purpose |
|---|---|
| `convert` | Convert one NEF or NMR-STAR restraint file. |
| `benchmark` | Run cases from a versioned benchmark-manifest YAML file. |

## Global flags

Global flags must appear before the command.

| Flag | Brief description | Default |
|---|---|---|
| `-h`, `--help` | Print top-level help and exit. | Not invoked. |
| `--version` | Print the installed NMR2Boltz version and exit. | Not invoked. |

## `convert`

```text
nmr2boltz convert INPUT [flags]
```

`INPUT` is required and must be a NEF or NMR-STAR file. Compressed inputs
supported by the parser may also be used.

| Flag | Brief description | Default |
|---|---|---|
| `-h`, `--help` | Print `convert` help and exit. | Not invoked. |
| `-o PATH`, `--output-dir PATH` | Directory for the conversion bundle. | `<input_stem>_nmr2boltz` beside the input. |
| `--format {auto,nef,nmr-star}` | Select the input format or detect it from STAR loop categories. | `auto` |
| `--averaging {sum-r6,mean-r6,hard-or}` | Select atom-set averaging semantics. | `sum-r6` |
| `--projection-margin ANGSTROM` | Add non-negative slack after triangle-inequality projection. | `0.0` Å |
| `--missing-upper-policy {reject,upper-linear,target-plus-uncertainty,target}` | Select how rows without an explicit upper bound are handled. | `reject` |
| `--pseudoatom-policy {reject,atomset}` | Reject geometric Q/M pseudoatoms or explicitly approximate them as atom sets. | `reject` |
| `--residue-map PATH` | Load a CSV, TSV, or JSON source-to-Boltz residue mapping. | Not set; use sequence-loop mapping. |
| `--target-yaml PATH` | Validate mapped chains, residues, modifications, and every exact/union-alternative position against an exact Boltz target before writing. | Not set; target-YAML validation is skipped. |
| `--chain-map SOURCE=BOLTZ` | Map one source chain code to a Boltz chain ID. Repeat for multiple mappings. | No mappings; preserve resolved source chain IDs. |
| `--allow-inferred-sequence-map` | Permit residue-order inference from restraints when no sequence loop exists. | `false` |
| `--ccd PATH` | Add an external CCD mmCIF file or directory for modified residues, ligands, or ions. Repeat as needed. | No external paths; use built-in and embedded topology only. |
| `--bond-length-config PATH` | Load JSON overrides for conservative X-H bond-length upper envelopes. | Not set; use built-in envelopes. |
| `--origin TYPE` | Include only a case-insensitive restraint origin/type such as `NOE`. Repeat to allow several types. | No filter; process all general distance-restraint origins. |
| `--boltz-min-distance ANGSTROM` | Set the minimum executable BoltzUI atom-contact distance. Smaller exact and union-alternative bounds are raised, which weakens them. | `2.0` Å |
| `--boltz-max-distance ANGSTROM` | Set the maximum executable BoltzUI atom-contact distance. Larger exact bounds are rejected; a union with any larger alternative is quarantined in full. Bounds are never clipped. | `20.0` Å |
| `--min-sequence-separation N` | Filter same-chain projected contacts whose residue-index separation is less than `N`. | `0` (disabled) |
| `--exclude-intraresidue` | Omit contacts within one mapped residue; all-intraresidue and mixed intra/inter-residue OR groups are omitted in full. | `false` |
| `--exclude-intrachain` | Emit only contacts between different mapped Boltz chain IDs; mixed intra/inter OR groups are omitted in full. | `false` |
| `--hypotheses N` | Generate up to `N` deterministic assignment-hypothesis YAML files from ambiguous OR groups. | `0` (none) |
| `--seed SEED` | Set the random seed used to choose assignment hypotheses. | `0` |
| `--strict` | After writing audits, return status 3 if any ambiguity or rejection remains. | `false` |

### `convert` effective-default notes

- `--format auto` detects NEF versus NMR-STAR from parsed categories, not from
  the filename extension alone.
- An empty `--origin` selection means all general distance-restraint origins;
  it does not mean NOE-only.
- `--min-sequence-separation 0`, `--exclude-intraresidue=false`, and
  `--exclude-intrachain=false` preserve all otherwise safe local and
  same-chain contacts.
- `--seed` has no effect unless `--hypotheses` is greater than zero.
- `--strict` does not suppress output. It writes the audit bundle first and
  then returns status 3 when unresolved ambiguity or rejection records exist.

## `benchmark`

```text
nmr2boltz benchmark MANIFEST [flags]
```

`MANIFEST` is required and must be a schema-version-1 benchmark YAML file.

| Flag | Brief description | Default |
|---|---|---|
| `-h`, `--help` | Print `benchmark` help and exit. | Not invoked. |
| `-o PATH`, `--output-dir PATH` | Directory for case outputs and `benchmark_summary.json`. | `<manifest_stem>_results` beside the manifest. |
| `--case ID` | Run only the named case. Repeat to select several cases. | No selection; run every manifest case. |

## Exit statuses

| Status | Meaning |
|---:|---|
| `0` | Requested operation completed successfully. |
| `2` | CLI input, parsing, topology, mapping, or target validation failed. |
| `3` | `convert --strict` completed its audit output but ambiguity or rejection remains. |
| `4` | One or more `benchmark` cases failed. |

The parser implementation in `src/nmr2boltz/cli.py` is the source of truth.
Run `nmr2boltz --help`, `nmr2boltz convert --help`, or
`nmr2boltz benchmark --help` to inspect the installed version directly.
