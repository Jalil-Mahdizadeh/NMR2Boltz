from __future__ import annotations

import argparse
import re
from pathlib import Path

from nmr2boltz.cli import build_parser


ROOT = Path(__file__).resolve().parents[1]
CLI_REFERENCE = ROOT / "docs" / "CLI_REFERENCE.md"


def _parser_flags() -> set[str]:
    parser = build_parser()
    parsers = [parser]
    parsers.extend(
        child
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
        for child in action.choices.values()
    )
    return {
        flag
        for current in parsers
        for action in current._actions
        for flag in action.option_strings
    }


def _documented_flag_rows(text: str) -> tuple[set[str], list[list[str]]]:
    rows: list[list[str]] = []
    flags: set[str] = set()
    for line in text.splitlines():
        if not line.startswith("| `"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) != 3:
            continue
        row_flags = set(re.findall(r"`(-{1,2}[a-z][a-z0-9-]*)", cells[0]))
        if not row_flags:
            continue
        flags.update(row_flags)
        rows.append(cells)
    return flags, rows


def test_cli_reference_covers_every_parser_flag_and_states_a_default():
    text = CLI_REFERENCE.read_text(encoding="utf-8")
    documented, rows = _documented_flag_rows(text)

    assert documented == _parser_flags()
    assert rows
    assert all(row[1] and row[2] for row in rows)


def test_cli_reference_covers_commands_positionals_and_exit_statuses():
    text = CLI_REFERENCE.read_text(encoding="utf-8")

    for token in (
        "`convert`",
        "`benchmark`",
        "`INPUT`",
        "`MANIFEST`",
        "| `0` |",
        "| `2` |",
        "| `3` |",
        "| `4` |",
    ):
        assert token in text
