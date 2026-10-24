#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Build Oral Oncology submission package:
  submission_package/
    Fig1.tif ... Fig4.tif          (300 dpi LZW TIFF, Elsevier)
    Supplementary_FigS1.tif
    Highlights.docx                (3-5 bullets, <=85 chars each)
    Table1.docx / Table2.docx / Table3.docx
    TitlePage.docx
    README_submission.md           (upload checklist)
No manuscript text is generated here - only mechanical packaging of existing content.
"""
import os
import re
import sys
from pathlib import Path

from PIL import Image
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT = Path(r"C:/Users/tidus/PycharmProjects/Stomatology-PVL")
FIGS = ROOT / "results" / "manuscript_figs"
MS = ROOT / "manuscript" / "submission_OralOncology_v1.md"
OUT = ROOT / "manuscript" / "submission_package"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- figures
FIG_MAP = [
    ("Fig1_severity_umap.png", "Fig1.tif"),
    ("Fig2_lincs_reversal.png", "Fig2.tif"),
    ("Fig3_bulk_hallmark.png", "Fig3.tif"),
    ("Fig4_second_cohort.png", "Fig4.tif"),
    ("PVL10_cluster_sigscore_bar.png", "Supplementary_FigS1.tif"),
]

def convert_figs():
    report = []
    for src_name, dst_name in FIG_MAP:
        src = FIGS / src_name
        if not src.exists():
            report.append(f"MISSING {src_name}")
            continue
        im = Image.open(src)
        info = im.info.get("dpi")
        if im.mode not in ("RGB", "L"):
            im = im.convert("RGB")
        dst = OUT / dst_name
        im.save(dst, format="TIFF", compression="tiff_lzw", dpi=(300, 300))
        chk = Image.open(dst)
        report.append(
            f"{dst_name}: {chk.size[0]}x{chk.size[1]} px, dpi={chk.info.get('dpi')}, "
            f"mode={chk.mode}, {dst.stat().st_size/1024:.0f} KB (source dpi={info})"
        )
    return report

# ---------------------------------------------------------------- manuscript parsing
def read_ms():
    return MS.read_text(encoding="utf-8")

PIPE = "\x00PIPE\x00"

def parse_tables(md):
    """Return list of (num, caption, header, rows, note) for the three tables.

    Escaped pipes (\\|logFC\\|) must be protected BEFORE splitting on '|',
    otherwise a cell containing an escaped pipe is silently split in two.
    A bare paragraph directly after a table is treated as that table's note.
    """
    tables = []
    pattern = r"\*\*Table (\d)\.\*\*\s*(.+?)\n\n((?:\|.*\n)+)(?:\n([^\n|*#][^\n]*)\n)?"
    for m in re.finditer(pattern, md):
        num = m.group(1)
        caption = m.group(2).strip()
        note = (m.group(4) or "").strip()
        # a horizontal rule or a stray marker is not a table note
        if note and set(note) <= {"-", "—", "*", "_", " "}:
            note = ""
        block = m.group(3).strip().replace("\\|", PIPE).split("\n")
        cells = [[c.strip().replace(PIPE, "|") for c in ln.strip().strip("|").split("|")]
                 for ln in block]
        widths = {len(r) for r in cells}
        if len(cells) < 3:
            continue
        if len(widths) != 1:
            raise SystemExit(
                f"Table {num}: inconsistent column counts {sorted(widths)} - refusing to build"
            )
        header, _sep, rows = cells[0], cells[1], cells[2:]
        tables.append((num, caption, header, rows, note))
    return tables

def parse_highlights(md):
    """Highlights live in the early draft (kept there); the submission draft
    no longer carries them because Elsevier wants a standalone file."""
    for src in (md, (ROOT / "manuscript" / "manuscript_v1.md").read_text(encoding="utf-8")):
        m = re.search(r"## Highlights\n\n((?:- .+\n)+)", src)
        if m:
            return [ln[2:].strip() for ln in m.group(1).strip().split("\n") if ln.startswith("- ")]
    return []

def clean_md(s):
    s = s.replace("\\|", "|")
    s = re.sub(r"\*\*(.+?)\*\*", r"\1", s)
    s = re.sub(r"\*(.+?)\*", r"\1", s)
    s = s.replace("<sup>", "").replace("</sup>", "")
    return s.strip()

SUP_RE = re.compile(r"<sup>(.*?)</sup>")

def add_rich_text(paragraph, text, bold=False, size=None):
    """Add text to a paragraph, rendering <sup>..</sup> as true superscript runs."""
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    pos = 0
    for m in SUP_RE.finditer(text):
        before = text[pos:m.start()]
        if before:
            r = paragraph.add_run(before)
            r.bold = bold
            if size:
                r.font.size = size
        r = paragraph.add_run(m.group(1))
        r.font.superscript = True
        r.bold = bold
        if size:
            r.font.size = size
        pos = m.end()
    tail = text[pos:]
    if tail:
        r = paragraph.add_run(tail)
        r.bold = bold
        if size:
            r.font.size = size

# ---------------------------------------------------------------- docx helpers
def set_normal(doc):
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(12)
    pf = st.paragraph_format
    pf.space_after = Pt(0)
    pf.line_spacing = 2.0

def add_table_doc(num, caption, header, rows, note, path):
    doc = Document()
    set_normal(doc)
    p = doc.add_paragraph()
    r = p.add_run(f"Table {num}. {clean_md(caption)}")
    r.bold = True
    doc.add_paragraph()
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    for i, h in enumerate(header):
        cell = t.rows[0].cells[i]
        cell.text = ""
        run = cell.paragraphs[0].add_run(clean_md(h))
        run.bold = True
        run.font.size = Pt(10)
    for row in rows:
        cells = t.add_row().cells
        for i, val in enumerate(row[: len(header)]):
            cells[i].text = ""
            run = cells[i].paragraphs[0].add_run(clean_md(val))
            run.font.size = Pt(10)
    if note:
        doc.add_paragraph()
        np_ = doc.add_paragraph()
        run = np_.add_run("Note. " + clean_md(note))
        run.font.size = Pt(10)
    doc.save(path)
    return f"Table{num}.docx: {len(rows)} data rows x {len(header)} cols" + (
        ", with note" if note else ", NO NOTE"
    )

def add_highlights(bullets, path):
    doc = Document()
    set_normal(doc)
    p = doc.add_paragraph()
    r = p.add_run("Highlights")
    r.bold = True
    r.font.size = Pt(14)
    doc.add_paragraph()
    for b in bullets:
        para = doc.add_paragraph(clean_md(b), style="List Bullet")
        para.paragraph_format.line_spacing = 1.5
        para.paragraph_format.space_after = Pt(6)
    doc.save(path)
    return [(b, len(clean_md(b))) for b in bullets]

def add_title_page(md, path):
    """Extract the TITLE PAGE block and write it as a separate Word file."""
    m = re.search(r"## TITLE PAGE.*?\n(.*?)\n---\n", md, re.S)
    body = m.group(1).strip() if m else ""
    body = re.sub(r"<!--.*?-->", "", body, flags=re.S)  # drop author-only HTML notes
    doc = Document()
    set_normal(doc)
    for block in [b for b in body.split("\n\n") if b.strip()]:
        lines = [ln for ln in block.split("\n") if ln.strip()]
        for ln in lines:
            txt = ln.strip()
            if not txt:
                continue
            p = doc.add_paragraph()
            is_title = txt.startswith("**Title:**")
            txt = txt.replace("**Title:**", "Title:") if is_title else txt
            add_rich_text(p, txt, bold=is_title, size=Pt(14) if is_title else None)
            p.paragraph_format.line_spacing = 1.5
            p.paragraph_format.space_after = Pt(6)
    doc.save(path)
    return path.name

# ---------------------------------------------------------------- main
def main():
    md = read_ms()
    lines = []
    lines.append("### Figures -> TIFF (300 dpi, LZW)")
    lines += ["- " + r for r in convert_figs()]

    tables = parse_tables(md)
    lines.append("\n### Tables -> DOCX")
    if not tables:
        lines.append("- NO TABLES PARSED")
    for num, caption, header, rows, note in tables:
        msg = add_table_doc(num, caption, header, rows, note, OUT / f"Table{num}.docx")
        lines.append("- " + msg)

    bullets = parse_highlights(md)
    lines.append("\n### Highlights -> DOCX")
    if not bullets:
        lines.append("- NO HIGHLIGHTS PARSED")
    else:
        info = add_highlights(bullets, OUT / "Highlights.docx")
        for b, n in info:
            flag = "OK" if n <= 85 else "OVER LIMIT"
            lines.append(f"- {n:>2} chars [{flag}] {b}")

    lines.append("\n### Title page -> DOCX")
    lines.append("- " + add_title_page(md, OUT / "TitlePage.docx"))

    print("\n".join(lines))

if __name__ == "__main__":
    main()
