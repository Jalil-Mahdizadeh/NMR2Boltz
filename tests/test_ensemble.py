from validation.compare_ensemble import align_structure_to_sequence_map


def test_structure_alignment_rekeys_author_numbers_to_boltz_indices():
    sequence_map = [
        {"boltz_chain": "A", "boltz_residue_index": 1, "residue_name": "ALA"},
        {"boltz_chain": "A", "boltz_residue_index": 2, "residue_name": "GLY"},
        {"boltz_chain": "A", "boltz_residue_index": 3, "residue_name": "SER"},
    ]
    models = [
        {
            ("A", 10, "CA"): (1.0, 2.0, 3.0),
            ("A", 12, "CA"): (2.0, 3.0, 4.0),
            ("A", 13, "CA"): (3.0, 4.0, 5.0),
        }
    ]
    residues = [{("A", 10): "ALA", ("A", 12): "GLY", ("A", 13): "SER"}]

    aligned, summary = align_structure_to_sequence_map(sequence_map, models, residues)

    assert set(aligned[0]) == {("A", 1, "CA"), ("A", 2, "CA"), ("A", 3, "CA")}
    assert summary[0]["exact_component_matches"] == 3
    assert summary[0]["unmapped_target_residues"] == 0


def test_structure_alignment_accepts_equivalent_nucleic_acid_names():
    sequence_map = [
        {"boltz_chain": "X", "boltz_residue_index": 1, "residue_name": "DA"},
        {"boltz_chain": "X", "boltz_residue_index": 2, "residue_name": "DC"},
    ]
    models = [{("X", 5, "P"): (0.0, 0.0, 0.0), ("X", 6, "P"): (1.0, 0.0, 0.0)}]
    residues = [{("X", 5): "A", ("X", 6): "C"}]

    _aligned, summary = align_structure_to_sequence_map(sequence_map, models, residues)

    assert summary[0]["exact_component_matches"] == 0
    assert summary[0]["compatible_component_matches"] == 2
