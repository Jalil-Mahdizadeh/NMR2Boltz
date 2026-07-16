from pathlib import Path

import pytest
import yaml

from nmr2boltz.cli import main
from nmr2boltz.model import ConversionReport, SequenceRecord
from nmr2boltz.project import ProjectionSettings, project_document
from nmr2boltz.star import parse_star_document
from nmr2boltz.target import validate_report_against_target
from nmr2boltz.topology import TopologyLibrary


FIXTURES = Path(__file__).parent / "fixtures"


def _report():
    path = FIXTURES / "example.nef"
    return project_document(
        parse_star_document(path),
        input_file=str(path),
        topology_library=TopologyLibrary(),
        settings=ProjectionSettings(),
    )


def _write_target(path: Path, sequence: str) -> None:
    path.write_text(
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: [A, B]\n"
        f"      sequence: {sequence}\n"
        "      msa: empty\n",
        encoding="utf-8",
    )


def _minimal_report(residue_names: list[str]) -> ConversionReport:
    return ConversionReport(
        input_file="synthetic",
        detected_format="nef",
        settings={},
        sequence_map=[
            SequenceRecord(
                source_chain="A",
                source_sequence_code=str(index),
                residue_name=residue,
                boltz_chain="A",
                boltz_residue_index=index,
                source="synthetic",
            )
            for index, residue in enumerate(residue_names, start=1)
        ],
        emitted_constraints=[],
        ambiguous_groups=[],
        rejections=[],
        warnings=[],
        statistics={},
    )


def test_matching_boltz_target_validates_sequence_and_constraints(tmp_path):
    target = tmp_path / "target.yaml"
    _write_target(target, "VAYLG")

    result = validate_report_against_target(_report(), target)

    assert not result.errors
    assert result.checked_sequence_records == 5
    assert result.checked_constraints > 0
    assert result.mapped_positions == 5
    assert result.target_chains == ["A", "B"]


def test_target_residue_mismatch_is_an_error(tmp_path):
    target = tmp_path / "target.yaml"
    _write_target(target, "VAFLG")

    result = validate_report_against_target(_report(), target)

    mismatch = next(issue for issue in result.errors if issue.code == "target_residue_mismatch")
    assert mismatch.chain == "A"
    assert mismatch.residue_index == 3
    assert "TYR" in mismatch.message


def test_convert_does_not_write_constraints_when_target_validation_fails(tmp_path):
    target = tmp_path / "target.yaml"
    output = tmp_path / "output"
    _write_target(target, "VAFLG")

    status = main(
        [
            "convert",
            str(FIXTURES / "example.nef"),
            "--target-yaml",
            str(target),
            "--output-dir",
            str(output),
        ]
    )

    assert status == 2
    assert not (output / "boltz_constraints.yaml").exists()


@pytest.mark.parametrize(
    ("entity_type", "sequence", "residues"),
    [
        ("rna", "ACGU", ["A", "C", "G", "U"]),
        ("dna", "ACGT", ["DA", "DC", "DG", "DT"]),
    ],
)
def test_nucleic_acid_target_identity_validation(
    tmp_path, entity_type, sequence, residues
):
    target = tmp_path / "target.yaml"
    target.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "sequences": [{entity_type: {"id": "A", "sequence": sequence}}],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    result = validate_report_against_target(_minimal_report(residues), target)

    assert not result.errors


def test_declared_polymer_modification_must_match_source_ccd(tmp_path):
    target = tmp_path / "target.yaml"
    target.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "sequences": [
                    {
                        "protein": {
                            "id": "A",
                            "sequence": "M",
                            "modifications": [{"position": 1, "ccd": "MSE"}],
                        }
                    }
                ],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    matching = validate_report_against_target(_minimal_report(["MSE"]), target)
    mismatching = validate_report_against_target(_minimal_report(["MET"]), target)

    assert not matching.errors
    assert [issue.code for issue in mismatching.errors] == ["target_modification_mismatch"]
