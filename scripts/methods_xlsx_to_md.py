#!/usr/bin/env python3
"""Render the organizers' filled methods workbook as Markdown.

The organizers supply the methods description as an **Excel template**, and their
instructions say to *"use the provided methods description template to organize
your analysis, then export your report as a PDF or Markdown file for
submission."* So the workbook is the scaffold and Markdown is an accepted export.

🔴 **This script exists so the two cannot drift.** Keeping a hand-written Markdown
copy alongside a filled spreadsheet guarantees they disagree eventually, and the
disagreement is invisible to anyone who can only open one of them. The Markdown is
generated *from* the workbook, every time, and the workbook is the single source.

Args:
    --xlsx: the filled methods workbook.
    --sheet: which sheet to render (the workbook carries one per track).
    --out: the Markdown file to write.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import openpyxl


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", type=Path, required=True)
    ap.add_argument("--sheet", default="Methods")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    workbook = openpyxl.load_workbook(args.xlsx)
    sheet = workbook[args.sheet]

    title = sheet["A1"].value or args.sheet

    lines: list[str] = [f"# {title}", ""]
    lines.append("---")
    lines.append("")

    answered = 0
    for row in range(7, sheet.max_row + 1):
        question = sheet[f"A{row}"].value
        if not question:
            continue
        answer = sheet[f"B{row}"].value
        required = (sheet[f"A{row}"].font.color or None) and getattr(
            sheet[f"A{row}"].font.color, "rgb", None
        ) == "FFFF0000"
        heading = str(question).strip()
        lines.append(f"## {'🔴 ' if required else ''}{heading}")
        lines.append("")
        if answer:
            answered += 1
            lines.append(str(answer).strip())
        else:
            lines.append("*(not answered)*")
        lines.append("")

    args.out.write_text("\n".join(lines).rstrip() + "\n")
    print(f"{answered} answers -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
