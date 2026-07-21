from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import yaml

from nmr2boltz.output import _outward_decimal
from validation.distance_check_token import build_case, load_token_bounds, run


Token = tuple[str, int]
Contact = tuple[Token, Token, float]


def _write_format(root: Path, format_name: str, contacts: list[Contact]) -> None:
    directory = root / format_name
    directory.mkdir(parents=True)
    sequence_map = [
        {
            "source_chain": "A",
            "source_sequence_code": str(9 + index),
            "residue_name": residue_name,
            "boltz_chain": "A",
            "boltz_residue_index": index,
            "aliases": [["A", str(9 + index), residue_name]],
        }
        for index, residue_name in enumerate(("ALA", "GLY", "SER"), start=1)
    ]
    token_constraints = [
        {
            "token1": {"chain": token1[0], "residue_index": token1[1]},
            "token2": {"chain": token2[0], "residue_index": token2[1]},
            "max_distance": bound,
        }
        for token1, token2, bound in contacts
    ]
    (directory / "conversion_report.json").write_text(
        json.dumps(
            {"sequence_map": sequence_map, "token_constraints": token_constraints}
        ),
        encoding="utf-8",
    )
    constraints = [
        {
            "contact": {
                "token1": list(token1),
                "token2": list(token2),
                "max_distance": float(_outward_decimal(bound)),
                "force": False,
            }
        }
        for token1, token2, bound in contacts
    ]
    (directory / "token_constraints.yaml").write_text(
        yaml.safe_dump({"constraints": constraints}, sort_keys=False),
        encoding="utf-8",
    )


def _write_pdb(path: Path) -> None:
    path.write_text(
        """MODEL        1
ATOM      1  CA  ALA A  10       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  H   ALA A  10       3.900   0.000   0.000  1.00 20.00           H
ATOM      3  CA  GLY A  11       4.000   0.000   0.000  1.00 20.00           C
ATOM      4  CA  SER A  12      10.000   0.000   0.000  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      5  CA  ALA A  10       0.000   0.000   0.000  1.00 20.00           C
ATOM      6  H   ALA A  10       7.900   0.000   0.000  1.00 20.00           H
ATOM      7  CA  GLY A  11       8.000   0.000   0.000  1.00 20.00           C
ATOM      8  CA  SER A  12      12.000   0.000   0.000  1.00 20.00           C
ENDMDL
END
""",
        encoding="utf-8",
    )


def test_build_case_writes_token_pairs_with_minimum_heavy_atom_distances(
    tmp_path: Path,
) -> None:
    case = tmp_path / "input" / "TEST"
    case.mkdir(parents=True)
    _write_pdb(case / "TEST.pdb")
    conversion = tmp_path / "conversion" / "TEST"
    common_nef = (("A", 1), ("A", 2), 5.0000001)
    common_star = (("A", 1), ("A", 2), 4.5)
    nef_only = (("A", 1), ("A", 3), 11.0)
    star_only = (("A", 2), ("A", 3), 7.0)
    _write_format(conversion, "nef", [common_nef, nef_only])
    _write_format(conversion, "star", [common_star, star_only])

    result = build_case(
        case,
        tmp_path / "conversion",
        tmp_path / "distance_check_token",
        mapping_tolerance=1e-6,
        satisfaction_tolerance=0.01,
    )

    csv_path = tmp_path / "distance_check_token" / "TEST.csv"
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows == [
        {
            "residue_pair": "A:1--A:2",
            "nef_max_distance_A": "5.000001",
            "star_max_distance_A": "4.500000",
            "model_1_min_heavy_atom_distance_A": "4.000000",
            "model_2_min_heavy_atom_distance_A": "8.000000",
        },
        {
            "residue_pair": "A:1--A:3",
            "nef_max_distance_A": "11.000000",
            "star_max_distance_A": "",
            "model_1_min_heavy_atom_distance_A": "10.000000",
            "model_2_min_heavy_atom_distance_A": "12.000000",
        },
        {
            "residue_pair": "A:2--A:3",
            "nef_max_distance_A": "",
            "star_max_distance_A": "7.000000",
            "model_1_min_heavy_atom_distance_A": "6.000000",
            "model_2_min_heavy_atom_distance_A": "4.000000",
        },
    ]
    assert result["model_count"] == 2
    assert result["common_pairs"] == 1
    assert result["nef_only_pairs"] == 1
    assert result["star_only_pairs"] == 1
    assert result["common_pairs_with_different_bounds"] == 1
    assert result["missing_distance_cells"] == 0


def test_empty_token_files_still_write_model_columns(tmp_path: Path) -> None:
    case = tmp_path / "input" / "EMPTY"
    case.mkdir(parents=True)
    _write_pdb(case / "EMPTY.pdb")
    conversion = tmp_path / "conversion" / "EMPTY"
    _write_format(conversion, "nef", [])
    _write_format(conversion, "star", [])

    result = build_case(
        case,
        tmp_path / "conversion",
        tmp_path / "distance_check_token",
        mapping_tolerance=1e-6,
        satisfaction_tolerance=0.01,
    )

    text = (tmp_path / "distance_check_token" / "EMPTY.csv").read_text(
        encoding="utf-8"
    )
    assert text == (
        "residue_pair,nef_max_distance_A,star_max_distance_A,"
        "model_1_min_heavy_atom_distance_A,"
        "model_2_min_heavy_atom_distance_A\n"
    )
    assert result["row_count"] == 0


def test_token_yaml_and_report_must_match(tmp_path: Path) -> None:
    contact = (("A", 1), ("A", 2), 5.0)
    _write_format(tmp_path, "nef", [contact])
    yaml_path = tmp_path / "nef" / "token_constraints.yaml"
    payload = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    payload["constraints"][0]["contact"]["max_distance"] = 4.999999
    yaml_path.write_text(yaml.safe_dump(payload), encoding="utf-8")

    try:
        load_token_bounds(tmp_path / "nef")
    except ValueError as exc:
        assert "YAML/report bound mismatch" in str(exc)
    else:
        raise AssertionError("stale or tightened executable token YAML was accepted")


def test_token_yaml_must_use_force_false(tmp_path: Path) -> None:
    contact = (("A", 1), ("A", 2), 5.0)
    _write_format(tmp_path, "nef", [contact])
    yaml_path = tmp_path / "nef" / "token_constraints.yaml"
    payload = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    payload["constraints"][0]["contact"]["force"] = True
    yaml_path.write_text(yaml.safe_dump(payload), encoding="utf-8")

    try:
        load_token_bounds(tmp_path / "nef")
    except ValueError as exc:
        assert "Token contact is forced" in str(exc)
    else:
        raise AssertionError("forced token contact was accepted")


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
