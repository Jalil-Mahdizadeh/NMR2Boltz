from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import yaml

from nmr2boltz.output import _outward_decimal
from validation.distance_check import (
    _merge_distances,
    build_case,
    load_exact_bounds,
    run,
)


def _atom(chain: str, residue: int, name: str) -> dict[str, object]:
    return {"chain": chain, "residue_index": residue, "atom_name": name}


def _write_format(
    root: Path,
    format_name: str,
    contacts: list[tuple[dict[str, object], dict[str, object], float]],
) -> None:
    directory = root / format_name
    directory.mkdir(parents=True)
    sequence_map = [
        {
            "source_chain": "A",
            "source_sequence_code": "10",
            "residue_name": "ALA",
            "boltz_chain": "A",
            "boltz_residue_index": 1,
            "aliases": [["A", "10", "ALA"]],
        },
        {
            "source_chain": "A",
            "source_sequence_code": "11",
            "residue_name": "GLY",
            "boltz_chain": "A",
            "boltz_residue_index": 2,
            "aliases": [["A", "11", "GLY"]],
        },
    ]
    emitted = [
        {"atom1": atom1, "atom2": atom2, "max_distance": bound}
        for atom1, atom2, bound in contacts
    ]
    (directory / "conversion_report.json").write_text(
        json.dumps({"sequence_map": sequence_map, "emitted_constraints": emitted}),
        encoding="utf-8",
    )
    constraints = [
        {
            "atom_contact": {
                "atom1": [
                    contact[0]["chain"],
                    contact[0]["residue_index"],
                    contact[0]["atom_name"],
                ],
                "atom2": [
                    contact[1]["chain"],
                    contact[1]["residue_index"],
                    contact[1]["atom_name"],
                ],
                "max_distance": float(_outward_decimal(contact[2])),
                "force": True,
            }
        }
        for contact in contacts
    ]
    (directory / "atom_constraints_exact.yaml").write_text(
        yaml.safe_dump({"constraints": constraints}, sort_keys=False),
        encoding="utf-8",
    )


def _write_pdb(path: Path) -> None:
    path.write_text(
        """MODEL        1
ATOM      1  CA  ALA A  10       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  CB  ALA A  10       0.000   3.000   0.000  1.00 20.00           C
ATOM      3  CA  GLY A  11       4.000   0.000   0.000  1.00 20.00           C
ATOM      7  N   GLY A  11       3.000   0.000   0.000  1.00 20.00           N
ENDMDL
MODEL        2
ATOM      4  CA  ALA A  10       0.000   0.000   0.000  1.00 20.00           C
ATOM      5  CB  ALA A  10       0.000   6.000   0.000  1.00 20.00           C
ATOM      6  CA  GLY A  11       8.000   0.000   0.000  1.00 20.00           C
ATOM      8  N   GLY A  11       7.000   0.000   0.000  1.00 20.00           N
ENDMDL
END
""",
        encoding="utf-8",
    )


def test_build_case_writes_union_of_exact_pairs_with_model_distances(tmp_path: Path) -> None:
    case = tmp_path / "input" / "TEST"
    case.mkdir(parents=True)
    _write_pdb(case / "TEST.pdb")
    conversion = tmp_path / "conversion" / "TEST"
    common = (_atom("A", 1, "CA"), _atom("A", 2, "CA"), 5.0000001)
    nef_only = (_atom("A", 1, "CB"), _atom("A", 2, "CA"), 6.0)
    star_only = (_atom("A", 1, "CA"), _atom("A", 1, "CB"), 4.0)
    _write_format(conversion, "nef", [common, nef_only])
    _write_format(conversion, "star", [common, star_only])

    result = build_case(
        case,
        tmp_path / "conversion",
        tmp_path / "distance_check",
        mapping_tolerance=1e-6,
        satisfaction_tolerance=0.01,
    )

    with (tmp_path / "distance_check" / "TEST.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["heavy_atom_pair"] for row in rows] == sorted(
        row["heavy_atom_pair"] for row in rows
    )
    assert rows == [
        {
            "heavy_atom_pair": "A:1:CA--A:1:CB",
            "nef_max_distance_A": "",
            "star_max_distance_A": "4.000000",
            "model_1_distance_A": "3.000000",
            "model_2_distance_A": "6.000000",
        },
        {
            "heavy_atom_pair": "A:1:CA--A:2:CA",
            "nef_max_distance_A": "5.000001",
            "star_max_distance_A": "5.000001",
            "model_1_distance_A": "4.000000",
            "model_2_distance_A": "8.000000",
        },
        {
            "heavy_atom_pair": "A:1:CB--A:2:CA",
            "nef_max_distance_A": "6.000000",
            "star_max_distance_A": "",
            "model_1_distance_A": "5.000000",
            "model_2_distance_A": "10.000000",
        },
    ]
    assert result["model_count"] == 2
    assert result["common_pairs"] == 1
    assert result["nef_only_pairs"] == 1
    assert result["star_only_pairs"] == 1
    assert result["missing_distance_cells"] == 0


def test_empty_exact_files_still_write_model_columns(tmp_path: Path) -> None:
    case = tmp_path / "input" / "EMPTY"
    case.mkdir(parents=True)
    _write_pdb(case / "EMPTY.pdb")
    conversion = tmp_path / "conversion" / "EMPTY"
    _write_format(conversion, "nef", [])
    _write_format(conversion, "star", [])

    result = build_case(
        case,
        tmp_path / "conversion",
        tmp_path / "distance_check",
        mapping_tolerance=1e-6,
        satisfaction_tolerance=0.01,
    )

    text = (tmp_path / "distance_check" / "EMPTY.csv").read_text(encoding="utf-8")
    assert text == (
        "heavy_atom_pair,nef_max_distance_A,star_max_distance_A,"
        "model_1_distance_A,model_2_distance_A\n"
    )
    assert result["row_count"] == 0


def test_missing_coordinate_is_blank_not_invalid_atom_evidence(tmp_path: Path) -> None:
    case = tmp_path / "input" / "PARTIAL"
    case.mkdir(parents=True)
    _write_pdb(case / "PARTIAL.pdb")
    pdb_text = (case / "PARTIAL.pdb").read_text(encoding="utf-8")
    (case / "PARTIAL.pdb").write_text(
        pdb_text.replace(
            "ATOM      6  CA  GLY A  11       8.000   0.000   0.000  1.00 20.00           C\n",
            "",
        ),
        encoding="utf-8",
    )
    conversion = tmp_path / "conversion" / "PARTIAL"
    contact = (_atom("A", 1, "CA"), _atom("A", 2, "CA"), 5.0)
    _write_format(conversion, "nef", [contact])
    _write_format(conversion, "star", [contact])

    result = build_case(
        case,
        tmp_path / "conversion",
        tmp_path / "distance_check",
        mapping_tolerance=1e-6,
        satisfaction_tolerance=0.01,
    )

    with (tmp_path / "distance_check" / "PARTIAL.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        row = next(csv.DictReader(handle))
    assert row["model_1_distance_A"] == "4.000000"
    assert row["model_2_distance_A"] == ""
    assert result["missing_distance_cells"] == 1
    assert result["pairs_missing_any_model"] == 1


def test_exact_yaml_and_report_must_match(tmp_path: Path) -> None:
    root = tmp_path / "nef"
    contact = (_atom("A", 1, "CA"), _atom("A", 2, "CA"), 5.0)
    _write_format(tmp_path, "nef", [contact])
    payload = yaml.safe_load((root / "atom_constraints_exact.yaml").read_text())
    payload["constraints"][0]["atom_contact"]["max_distance"] = 4.999999
    (root / "atom_constraints_exact.yaml").write_text(yaml.safe_dump(payload))

    try:
        load_exact_bounds(root)
    except ValueError as exc:
        assert "YAML/report bound mismatch" in str(exc)
    else:
        raise AssertionError("stale or tightened executable YAML was accepted")


def test_common_pair_mapping_mismatch_fails_closed() -> None:
    pair = (("A", 1, "CA"), ("A", 2, "CB"))
    try:
        _merge_distances(pair, [3.0, None], [3.0, 4.0], 1e-6)
    except ValueError as exc:
        assert "coordinate-resolution mismatch" in str(exc)
    else:
        raise AssertionError("format-dependent coordinate resolution was accepted")


def test_run_requires_exact_expected_corpus_size(tmp_path: Path) -> None:
    input_directory = tmp_path / "input"
    input_directory.mkdir()
    args = SimpleNamespace(
        input_directory=input_directory,
        conversion_output=tmp_path / "conversion",
        output_directory=tmp_path / "output",
        expected_case_count=12,
        mapping_tolerance=1e-6,
        satisfaction_tolerance=0.01,
    )
    try:
        run(args)
    except ValueError as exc:
        assert "Expected 12 benchmark entries" in str(exc)
    else:
        raise AssertionError("incomplete benchmark corpus was accepted")
