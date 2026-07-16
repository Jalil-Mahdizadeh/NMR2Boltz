#!/usr/bin/env python3
from __future__ import annotations

import html
import re
import sys
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    HRFlowable,
    PageTemplate,
    Paragraph,
    Preformatted,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs" / "SCIENTIFIC_METHOD.md"
OUTPUT = ROOT / "docs" / "SCIENTIFIC_METHOD.pdf"


def register_fonts() -> tuple[str, str, str]:
    candidates = [
        (
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"),
        ),
        (
            Path("C:/Windows/Fonts/arial.ttf"),
            Path("C:/Windows/Fonts/arialbd.ttf"),
            Path("C:/Windows/Fonts/consola.ttf"),
        ),
    ]
    for regular, bold, mono in candidates:
        if regular.exists() and bold.exists() and mono.exists():
            pdfmetrics.registerFont(TTFont("DocSans", str(regular)))
            pdfmetrics.registerFont(TTFont("DocSansBold", str(bold)))
            pdfmetrics.registerFont(TTFont("DocMono", str(mono)))
            return "DocSans", "DocSansBold", "DocMono"
    return "Helvetica", "Helvetica-Bold", "Courier"


REGULAR, BOLD, MONO = register_fonts()


def inline_markup(text: str) -> str:
    code_spans: list[str] = []

    def protect_code(match: re.Match[str]) -> str:
        token = f"@@NMR2BOLTZCODE{len(code_spans)}@@"
        code_spans.append(match.group(1))
        return token

    protected = re.sub(r"`([^`]+)`", protect_code, text)
    escaped = html.escape(protected)
    escaped = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", escaped)
    escaped = re.sub(r"\*([^*]+)\*", r"<i>\1</i>", escaped)
    for index, code in enumerate(code_spans):
        escaped = escaped.replace(
            f"@@NMR2BOLTZCODE{index}@@",
            f'<font name="{MONO}">{html.escape(code)}</font>',
        )
    return escaped


def make_styles():
    base = getSampleStyleSheet()
    body = ParagraphStyle(
        "Body",
        parent=base["BodyText"],
        fontName=REGULAR,
        fontSize=9.1,
        leading=12.4,
        spaceAfter=5,
        alignment=TA_LEFT,
        textColor=colors.HexColor("#1E293B"),
    )
    styles = {
        "title": ParagraphStyle(
            "Title",
            parent=base["Title"],
            fontName=BOLD,
            fontSize=20,
            leading=25,
            alignment=TA_CENTER,
            textColor=colors.HexColor("#0F172A"),
            spaceAfter=16,
        ),
        "h1": ParagraphStyle(
            "H1",
            parent=body,
            fontName=BOLD,
            fontSize=15,
            leading=19,
            textColor=colors.HexColor("#0F3D66"),
            spaceBefore=12,
            spaceAfter=7,
            keepWithNext=True,
        ),
        "h2": ParagraphStyle(
            "H2",
            parent=body,
            fontName=BOLD,
            fontSize=12.2,
            leading=15.5,
            textColor=colors.HexColor("#174A75"),
            spaceBefore=9,
            spaceAfter=5,
            keepWithNext=True,
        ),
        "h3": ParagraphStyle(
            "H3",
            parent=body,
            fontName=BOLD,
            fontSize=10.5,
            leading=13.5,
            textColor=colors.HexColor("#245B82"),
            spaceBefore=7,
            spaceAfter=4,
            keepWithNext=True,
        ),
        "body": body,
        "bullet": ParagraphStyle(
            "Bullet",
            parent=body,
            leftIndent=12,
            firstLineIndent=-7,
            bulletIndent=4,
            spaceAfter=2.5,
        ),
        "number": ParagraphStyle(
            "Number",
            parent=body,
            leftIndent=14,
            firstLineIndent=-10,
            spaceAfter=2.5,
        ),
        "code": ParagraphStyle(
            "Code",
            parent=body,
            fontName=MONO,
            fontSize=7.7,
            leading=10.2,
            leftIndent=5,
            rightIndent=5,
            borderColor=colors.HexColor("#CBD5E1"),
            borderWidth=0.5,
            borderPadding=6,
            backColor=colors.HexColor("#F8FAFC"),
            spaceBefore=4,
            spaceAfter=7,
        ),
        "equation": ParagraphStyle(
            "Equation",
            parent=body,
            fontName=MONO,
            fontSize=8.2,
            leading=11,
            leftIndent=12,
            rightIndent=12,
            alignment=TA_CENTER,
            textColor=colors.HexColor("#0F172A"),
            spaceBefore=4,
            spaceAfter=7,
        ),
        "table": ParagraphStyle(
            "TableCell",
            parent=body,
            fontSize=7.7,
            leading=9.7,
            spaceAfter=0,
        ),
        "table_header": ParagraphStyle(
            "TableHeader",
            parent=body,
            fontName=BOLD,
            fontSize=7.7,
            leading=9.7,
            textColor=colors.white,
            spaceAfter=0,
        ),
    }
    return styles


STYLES = make_styles()


def table_from_markdown(lines: list[str], available_width: float):
    rows: list[list[str]] = []
    for line in lines:
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        rows.append(cells)
    if len(rows) >= 2 and all(re.fullmatch(r":?-{3,}:?", cell.replace(" ", "")) for cell in rows[1]):
        rows.pop(1)
    column_count = max(len(row) for row in rows)
    normalized = [row + [""] * (column_count - len(row)) for row in rows]
    data = []
    for row_index, row in enumerate(normalized):
        style = STYLES["table_header"] if row_index == 0 else STYLES["table"]
        data.append([Paragraph(inline_markup(cell), style) for cell in row])
    widths = [available_width / column_count] * column_count
    table = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#315B78")),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#94A3B8")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F8FAFC")]),
            ]
        )
    )
    return table


def parse_markdown(text: str, available_width: float):
    lines = text.splitlines()
    story = []
    paragraph: list[str] = []
    code_lines: list[str] = []
    equation_lines: list[str] = []
    in_code = False
    in_equation = False

    def flush_paragraph():
        nonlocal paragraph
        if paragraph:
            joined = " ".join(part.strip() for part in paragraph)
            story.append(Paragraph(inline_markup(joined), STYLES["body"]))
            paragraph = []

    index = 0
    first_heading = True
    while index < len(lines):
        line = lines[index]
        stripped = line.strip()
        if in_code:
            if stripped.startswith("```"):
                story.append(Preformatted("\n".join(code_lines), STYLES["code"], maxLineLength=120))
                code_lines = []
                in_code = False
            else:
                code_lines.append(line)
            index += 1
            continue
        if in_equation:
            if stripped == r"\]":
                story.append(Preformatted("\n".join(equation_lines), STYLES["equation"], maxLineLength=120))
                equation_lines = []
                in_equation = False
            else:
                equation_lines.append(line)
            index += 1
            continue
        if stripped.startswith("```"):
            flush_paragraph()
            in_code = True
            index += 1
            continue
        if stripped == r"\[":
            flush_paragraph()
            in_equation = True
            index += 1
            continue
        if stripped.startswith("|") and stripped.endswith("|"):
            flush_paragraph()
            table_lines = []
            while index < len(lines) and lines[index].strip().startswith("|"):
                table_lines.append(lines[index])
                index += 1
            story.append(table_from_markdown(table_lines, available_width))
            story.append(Spacer(1, 6))
            continue
        heading = re.match(r"^(#{1,4})\s+(.*)$", stripped)
        if heading:
            flush_paragraph()
            level = len(heading.group(1))
            title = heading.group(2)
            if first_heading:
                story.append(Paragraph(inline_markup(title), STYLES["title"]))
                story.append(
                    Paragraph(
                        "Prepared for scientific review of NMR-guided Boltz structure prediction",
                        ParagraphStyle(
                            "Subtitle",
                            parent=STYLES["body"],
                            alignment=TA_CENTER,
                            textColor=colors.HexColor("#475569"),
                            spaceAfter=14,
                        ),
                    )
                )
                first_heading = False
            else:
                style = STYLES["h1"] if level == 1 else STYLES["h2"] if level == 2 else STYLES["h3"]
                story.append(Paragraph(inline_markup(title), style))
            index += 1
            continue
        if stripped in {"---", "***"}:
            flush_paragraph()
            story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#CBD5E1"), spaceBefore=5, spaceAfter=7))
            index += 1
            continue
        bullet = re.match(r"^-\s+(.*)$", stripped)
        if bullet:
            flush_paragraph()
            story.append(Paragraph(inline_markup(bullet.group(1)), STYLES["bullet"], bulletText="-"))
            index += 1
            continue
        numbered = re.match(r"^(\d+)\.\s+(.*)$", stripped)
        if numbered:
            flush_paragraph()
            story.append(
                Paragraph(
                    f"<b>{numbered.group(1)}.</b> {inline_markup(numbered.group(2))}",
                    STYLES["number"],
                )
            )
            index += 1
            continue
        if not stripped:
            flush_paragraph()
            index += 1
            continue
        paragraph.append(line)
        index += 1

    flush_paragraph()
    if code_lines:
        story.append(Preformatted("\n".join(code_lines), STYLES["code"]))
    if equation_lines:
        story.append(Preformatted("\n".join(equation_lines), STYLES["equation"]))
    return story


def header_footer(canvas, doc):
    canvas.saveState()
    width, height = A4
    canvas.setStrokeColor(colors.HexColor("#CBD5E1"))
    canvas.setLineWidth(0.4)
    canvas.line(18 * mm, height - 13 * mm, width - 18 * mm, height - 13 * mm)
    canvas.setFont("Helvetica", 7.3)
    canvas.setFillColor(colors.HexColor("#64748B"))
    canvas.drawRightString(width - 18 * mm, 10 * mm, f"Page {doc.page}")
    canvas.restoreState()


def main() -> int:
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else SOURCE
    output = Path(sys.argv[2]) if len(sys.argv) > 2 else OUTPUT
    output.parent.mkdir(parents=True, exist_ok=True)
    page_width, page_height = A4
    left = right = 18 * mm
    top = 18 * mm
    bottom = 16 * mm
    frame = Frame(left, bottom, page_width - left - right, page_height - top - bottom, id="normal")
    document = BaseDocTemplate(
        str(output),
        pagesize=A4,
        leftMargin=left,
        rightMargin=right,
        topMargin=top,
        bottomMargin=bottom,
        title="Scientific method for NMR-to-Boltz heavy-atom restraint projection",
        author="nmr2boltz project",
        subject="NMR-guided Boltz structure prediction",
    )
    document.addPageTemplates([PageTemplate(id="main", frames=[frame], onPageEnd=header_footer)])
    story = parse_markdown(source.read_text(encoding="utf-8"), frame._width)
    document.build(story)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
