# -*- coding: utf-8 -*-
"""
Build Editorial Manager-ready Manuscript.docx from submission_OralOncology_v1.md.

FORMAT-ONLY conversion: no wording, numbers, or punctuation are changed.
Deletions (author-facing material, per instruction):
  - file-header submission metadata block (Manuscript type / Target journal / Draft / Date)
  - the whole '## TITLE PAGE (submit as a separate file)' block (submitted separately)
  - the trailing '*Submission notes (for the author, delete before upload):*' block
  - all HTML comments (none present in this source, but handled)
"""
import re
import sys

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

SRC = r'C:/Users/tidus/PycharmProjects/Stomatology-PVL/manuscript/submission_OralOncology_v1.md'
OUT = r'C:/Users/tidus/PycharmProjects/Stomatology-PVL/manuscript/submission_package/Manuscript.docx'

FONT = 'Times New Roman'

# exact Times New Roman metrics for column sizing (advance widths from the
# installed font files; base-14 metrics mishandle en dashes, ×, ρ, superscripts)
_METRICS = {}


def _load_metrics(path):
    from fontTools.ttLib import TTFont
    f = TTFont(path, lazy=True)
    return f.getBestCmap(), f['head'].unitsPerEm, f['hmtx']


try:
    _METRICS['r'] = _load_metrics(r'C:\Windows\Fonts\times.ttf')
    _METRICS['b'] = _load_metrics(r'C:\Windows\Fonts\timesbd.ttf')

    def _text_len(s, bold, size=10.0):
        cmap, upem, hmtx = _METRICS['b' if bold else 'r']
        tot = 0
        for ch in s:
            gn = cmap.get(ord(ch))
            tot += hmtx[gn][0] if gn else int(0.5 * upem)
        return tot * size / upem
except Exception:                                     # pragma: no cover
    def _text_len(s, bold, size=10.0):
        return len(s) * (5.6 if bold else 5.3)

# --------------------------------------------------------------------------
# 1. Read source and apply deletions
# --------------------------------------------------------------------------
raw = open(SRC, encoding='utf-8').read()
raw = re.sub(r'<!--.*?-->', '', raw, flags=re.S)          # HTML comments

lines = raw.split('\n')
title_line = lines[0].strip()
assert title_line.startswith('# '), 'unexpected first line'
DOC_TITLE = title_line[2:].strip()

# body starts at '## Abstract' (drops metadata block + TITLE PAGE block)
body = raw.split('## Abstract', 1)[1]
body = '## Abstract' + body

# drop trailing submission-notes block
cut = body.find('*Submission notes (for the author, delete before upload):*')
assert cut != -1, 'submission notes marker not found'
body = body[:cut]
body = body.rstrip()
# strip a trailing horizontal rule left behind
body = re.sub(r'\n-{3,}\s*$', '', body).rstrip()

for forbidden in ('TITLE PAGE', 'Manuscript type:', 'Target journal:', 'Submission notes', '<!--'):
    assert forbidden not in body, forbidden

# --------------------------------------------------------------------------
# 2. Block parser
# --------------------------------------------------------------------------
MD_HR = re.compile(r'^-{3,}$')

def parse_blocks(text):
    """This source has no wrapped lines: every non-blank line is its own block."""
    blocks = []
    ls = text.split('\n')
    i = 0
    while i < len(ls):
        s = ls[i].strip()
        if not s or MD_HR.match(s):
            i += 1
            continue
        if s.startswith('### '):
            blocks.append(('h2', s[4:].strip()))
            i += 1
            continue
        if s.startswith('## '):
            blocks.append(('h1', s[3:].strip()))
            i += 1
            continue
        if s.startswith('# '):
            blocks.append(('title', s[2:].strip()))
            i += 1
            continue
        if s.startswith('|'):
            rows = []
            while i < len(ls) and ls[i].strip().startswith('|'):
                row = ls[i].strip()
                if not re.fullmatch(r'\|[\s\-:|]+\|', row):
                    cells = re.split(r'(?<!\\)\|', row)
                    cells = cells[1:-1] if row.startswith('|') and row.endswith('|') else cells
                    cells = [c.strip().replace('\\|', '|') for c in cells]
                    rows.append(cells)
                i += 1
            blocks.append(('table', rows))
            continue
        if s.startswith('- '):
            blocks.append(('bullet', s[2:].strip()))
            i += 1
            continue
        blocks.append(('p', s))
        i += 1
    return blocks

blocks = parse_blocks(body)

# --------------------------------------------------------------------------
# 3. Inline tokenizer -> runs
# --------------------------------------------------------------------------
TOKEN = re.compile(r'(\*\*.+?\*\*|\*[^*\n]+?\*|<sub>.+?</sub>|<sup>.+?</sup>)')
TAG = re.compile(r'<(sub|sup)>(.+?)</\1>')

def add_runs(paragraph, text, size=12, bold=False, italic=False):
    for tok in TOKEN.split(text):
        if not tok:
            continue
        b, it, sub, sup = bold, italic, False, False
        if tok.startswith('**') and tok.endswith('**') and len(tok) > 4:
            tok, b = tok[2:-2], True
        elif tok.startswith('*') and tok.endswith('*') and len(tok) > 2:
            tok, it = tok[1:-1], True
        m = TAG.fullmatch(tok)
        if m and not (b or it):
            pass
        # split out any sub/sup tags inside the remaining token
        pos = 0
        for m in TAG.finditer(tok):
            if m.start() > pos:
                r = paragraph.add_run(tok[pos:m.start()])
                r.font.name, r.font.size = FONT, Pt(size)
                r.bold, r.italic = b, it
            r = paragraph.add_run(m.group(2))
            r.font.name, r.font.size = FONT, Pt(size)
            r.bold, r.italic = b, it
            if m.group(1) == 'sub':
                r.font.subscript = True
            else:
                r.font.superscript = True
            pos = m.end()
        if pos < len(tok):
            r = paragraph.add_run(tok[pos:])
            r.font.name, r.font.size = FONT, Pt(size)
            r.bold, r.italic = b, it

# --------------------------------------------------------------------------
# 4. Document + styles
# --------------------------------------------------------------------------
doc = Document()

# page setup: A4, 2.54 cm margins
for sec in doc.sections:
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.top_margin = sec.bottom_margin = Cm(2.54)
    sec.left_margin = sec.right_margin = Cm(2.54)
    sec.header_distance, sec.footer_distance = Cm(1.25), Cm(1.25)

def set_font(style_or_run_font, size=12, bold=None, italic=None, color=None):
    f = style_or_run_font
    f.name = FONT
    f.size = Pt(size)
    if bold is not None:
        f.bold = bold
    if italic is not None:
        f.italic = italic
    if color is not None:
        f.color.rgb = color

def force_eastasia(style, name=FONT):
    rpr = style.element.get_or_add_rPr()
    rf = rpr.find(qn('w:rFonts'))
    if rf is None:
        rf = OxmlElement('w:rFonts')
        rpr.append(rf)
    for a in ('w:ascii', 'w:hAnsi', 'w:cs', 'w:eastAsia'):
        rf.set(qn(a), name)

normal = doc.styles['Normal']
set_font(normal.font, 12, color=RGBColor(0, 0, 0))
force_eastasia(normal)
pf = normal.paragraph_format
pf.line_spacing = 2.0
pf.space_before = Pt(0)
pf.space_after = Pt(0)
pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
pf.first_line_indent = Cm(0)

# Title style
tstyle = doc.styles['Title']
set_font(tstyle.font, 12, bold=True, color=RGBColor(0, 0, 0))
force_eastasia(tstyle)
tp = tstyle.paragraph_format
tp.line_spacing = 2.0
tp.space_before = Pt(0)
tp.space_after = Pt(12)
tp.alignment = WD_ALIGN_PARAGRAPH.LEFT
# remove the default bottom border / big spacing of the built-in Title style
tppr = tstyle.element.get_or_add_pPr()
for tag in ('w:pBdr',):
    for el in tppr.findall(qn(tag)):
        tppr.remove(el)

# Heading 1
h1 = doc.styles['Heading 1']
set_font(h1.font, 12, bold=True, italic=False, color=RGBColor(0, 0, 0))
force_eastasia(h1)
h1p = h1.paragraph_format
h1p.line_spacing = 2.0
h1p.space_before = Pt(12)
h1p.space_after = Pt(6)
h1p.alignment = WD_ALIGN_PARAGRAPH.LEFT
h1p.keep_with_next = True

# Heading 2
h2 = doc.styles['Heading 2']
set_font(h2.font, 12, bold=True, italic=True, color=RGBColor(0, 0, 0))
force_eastasia(h2)
h2p = h2.paragraph_format
h2p.line_spacing = 2.0
h2p.space_before = Pt(6)
h2p.space_after = Pt(0)
h2p.alignment = WD_ALIGN_PARAGRAPH.LEFT
h2p.keep_with_next = True

for st in (h1, h2, tstyle):
    rpr = st.element.get_or_add_rPr()
    for el in rpr.findall(qn('w:color')):
        rpr.remove(el)
    c = OxmlElement('w:color')
    c.set(qn('w:val'), '000000')
    rpr.append(c)

# --------------------------------------------------------------------------
# 5. Footer with centred page number
# --------------------------------------------------------------------------
footer = doc.sections[0].footer
fp = footer.paragraphs[0]
fp.text = ''
fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
fp.paragraph_format.line_spacing = 1.0
fp.paragraph_format.space_before = Pt(0)
fp.paragraph_format.space_after = Pt(0)


def add_field(paragraph, instr):
    r = paragraph.add_run()
    fld_begin = OxmlElement('w:fldChar')
    fld_begin.set(qn('w:fldCharType'), 'begin')
    instr_el = OxmlElement('w:instrText')
    instr_el.set(qn('xml:space'), 'preserve')
    instr_el.text = instr
    fld_sep = OxmlElement('w:fldChar')
    fld_sep.set(qn('w:fldCharType'), 'separate')
    txt = OxmlElement('w:t')
    txt.text = '1'
    fld_end = OxmlElement('w:fldChar')
    fld_end.set(qn('w:fldCharType'), 'end')
    r._r.append(fld_begin)
    r._r.append(instr_el)
    r._r.append(fld_sep)
    r._r.append(txt)
    r._r.append(fld_end)
    r.font.name = FONT
    r.font.size = Pt(12)


add_field(fp, ' PAGE ')

# --------------------------------------------------------------------------
# 6. Emit body
# --------------------------------------------------------------------------
def add_table(rows):
    ncol = max(len(r) for r in rows)
    rows = [r + [''] * (ncol - len(r)) for r in rows]
    table = doc.add_table(rows=len(rows), cols=ncol)
    table.style = doc.styles['Table Grid']
    table.alignment = WD_TABLE_ALIGNMENT.LEFT
    table.autofit = False

    # ---- measure columns with real Times metrics (header bold, body regular) ----
    PT_PER_CM = 28.3464567
    CELL_MAR_CM = 0.09         # explicit left/right cell margin
    TBL_IND_CM = CELL_MAR_CM   # border sits on the text margin: ind == cell margin
    PAD_PT = CELL_MAR_CM * PT_PER_CM
    avail_pt = (15.92 - TBL_IND_CM) * PT_PER_CM

    def tw(s, bold):
        try:
            return _text_len(s, bold)
        except Exception:
            return len(s) * (5.6 if bold else 5.3)

    # A header may wrap between words (never inside a hyphenated one).  A data
    # value should stay on one line where the table allows it, but long prose
    # cells may wrap -- at spaces, or at their own hyphens, never mid-word.
    MAX_ONELINE_PT = 200.0

    def _hyph(t):
        return ('-' in t) or ('\u2013' in t) or ('\u2014' in t)

    hard, soft, mass, wrap_cost = [], [], [], []
    for ci in range(ncol):
        tok_max = 0.0
        one_line = 0.0
        mx = 0.0
        for ri, row in enumerate(rows):
            cell = row[ci].strip()
            bold = (ri == 0)
            toks = [t for t in cell.split(' ') if t]
            full = tw(cell, bold)
            tok_w = max([tw(t, bold) for t in toks] or [0.0])
            need = tok_w
            if bold:
                # Word fills a line before it pushes a hyphenated word to the
                # next one, so a header needs room for "prev + hyphenated word"
                for i, t in enumerate(toks):
                    if _hyph(t):
                        need = max(need, tw(t, bold))
                        if i > 0:
                            need = max(need, tw(toks[i - 1] + ' ' + t, bold))
            mx = max(mx, full)
            tok_max = max(tok_max, need)
            if bold:
                one_line = max(one_line, need)        # headers wrap at spaces
            else:
                one_line = max(one_line, min(full, MAX_ONELINE_PT))
        hard.append(tok_max + 2 * PAD_PT + 0.75)
        soft.append(max(one_line, tok_max) + 2 * PAD_PT + 0.75)
        mass.append(max(mx, tok_max) + 2 * PAD_PT)
        # wrapping a prose cell is cheap; splitting a statistic is not
        data = ' '.join(r[ci] for r in rows[1:]).replace(' ', '')
        alpha = sum(ch.isalpha() for ch in data) / max(1, len(data))
        wrap_cost.append(1.0 if alpha >= 0.30 else 0.15)

    def water_fill(floor, weight, total):
        """Share `total` in proportion to `weight`, never below `floor`."""
        out = [0.0] * len(floor)
        active = set(i for i in range(len(floor)))
        rem = total
        while active:
            s = sum(weight[i] for i in active)
            below = [i for i in active if rem * weight[i] / s < floor[i] - 1e-9]
            if not below:
                for i in active:
                    out[i] = rem * weight[i] / s
                return out
            for i in below:
                out[i] = floor[i]
                rem -= floor[i]
                active.discard(i)
            if not active:
                break
        if sum(floor) > 0:
            return [f * total / sum(floor) for f in floor]
        return [total / len(floor)] * len(floor)

    if sum(soft) <= avail_pt:
        widths = list(soft)
        slack = avail_pt - sum(soft)
        grow = [max(0.0, mass[i] - soft[i]) for i in range(ncol)]
        g = sum(grow)
        if g > 0:
            widths = [widths[i] + slack * grow[i] / g for i in range(ncol)]
        else:
            widths = [w + slack / ncol for w in widths]
    else:
        widths = water_fill(hard, mass, avail_pt)
        # statistics must not wrap: give numeric columns their one-line width
        # back, paid for by prose columns that still have slack above their floor
        for i in sorted(range(ncol), key=lambda k: wrap_cost[k]):
            if widths[i] < soft[i] - 0.5 and wrap_cost[i] < 0.5:
                need = soft[i] - widths[i]
                while need > 0.5:
                    # prose columns first, then anything still above its own target
                    donors = [j for j in range(ncol)
                              if j != i and wrap_cost[j] >= 0.5 and widths[j] - hard[j] > 0.5]
                    if not donors:
                        donors = [j for j in range(ncol)
                                  if j != i and widths[j] - max(hard[j], soft[j]) > 0.5]
                    if not donors:
                        break
                    spare = sum(widths[j] - max(hard[j], soft[j] if wrap_cost[j] < 0.5 else hard[j])
                                for j in donors)
                    if spare <= 0.5:
                        break
                    take = min(need, spare)
                    for j in donors:
                        base = max(hard[j], soft[j] if wrap_cost[j] < 0.5 else hard[j])
                        widths[j] -= take * (widths[j] - base) / spare
                    widths[i] += take
                    need -= take
    if __import__('os').environ.get('ARS_TABLE_DEBUG'):
        print('   avail=%.1f hard=%s soft=%s -> widths=%s cost=%s' % (
            avail_pt, [round(x, 1) for x in hard], [round(x, 1) for x in soft],
            [round(x, 1) for x in widths], wrap_cost))
    if sum(widths) <= 0:
        widths = [avail_pt / ncol] * ncol
    widths = [w * avail_pt / sum(widths) for w in widths]
    widths_cm = [w / PT_PER_CM for w in widths]

    # explicit table indent so the left border sits on the text margin
    tblPr = table._tbl.tblPr
    ind = OxmlElement('w:tblInd')
    ind.set(qn('w:w'), str(int(round(TBL_IND_CM * 567))))   # cm -> twips
    ind.set(qn('w:type'), 'dxa')
    tblPr.append(ind)
    layout = OxmlElement('w:tblLayout')
    layout.set(qn('w:type'), 'fixed')
    tblPr.append(layout)
    # explicit cell margins (left/right 0.10 cm; top/bottom 0)
    mar = OxmlElement('w:tblCellMar')
    for side, val in (('top', 0), ('left', 51), ('bottom', 0), ('right', 51)):
        el = OxmlElement('w:' + side)
        el.set(qn('w:w'), str(val))
        el.set(qn('w:type'), 'dxa')
        mar.append(el)
    tblPr.append(mar)
    # grid definition (authoritative under fixed layout)
    for gc, w in zip(table._tbl.tblGrid.gridCol_lst, widths_cm):
        gc.set(qn('w:w'), str(int(round(w * 567))))

    for ri, row in enumerate(rows):
        for ci, val in enumerate(row):
            cell = table.cell(ri, ci)
            cell.width = Cm(widths_cm[ci])
            p = cell.paragraphs[0]
            p.paragraph_format.line_spacing = 1.0
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.first_line_indent = Cm(0)
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            add_runs(p, val, size=10, bold=(ri == 0))
    # header row: repeat across pages, keep together
    trPr = table.rows[0]._tr.get_or_add_trPr()
    th = OxmlElement('w:tblHeader')
    th.set(qn('w:val'), 'true')
    trPr.append(th)
    for ci, w in enumerate(widths_cm):
        for cell in table.columns[ci].cells:
            cell.width = Cm(w)
    return table


first_block = True
prev_kind = None

# Title first (taken from the source H1; the H1 itself lives above the body slice)
p = doc.add_paragraph(style='Title')
tr = p.add_run(DOC_TITLE)
tr.font.name = FONT
tr.font.size = Pt(12)
tr.bold = True
tr.italic = False

for kind, payload in blocks:
    if kind == 'title':
        continue
    elif kind == 'h1':
        p = doc.add_paragraph(style='Heading 1')
        add_runs(p, payload, size=12, bold=True)
        if payload.strip() == 'Tables':
            # tables begin on a fresh page so no table is split where avoidable
            p.paragraph_format.page_break_before = True
    elif kind == 'h2':
        p = doc.add_paragraph(style='Heading 2')
        add_runs(p, payload, size=12, bold=True, italic=True)
    elif kind == 'table':
        add_table(payload)
        # spacer paragraph after table so following text is not glued to it
        sp = doc.add_paragraph()
        sp.paragraph_format.line_spacing = 1.0
        sp.paragraph_format.space_before = Pt(0)
        sp.paragraph_format.space_after = Pt(6)
        for r in sp.runs:
            r.font.size = Pt(6)
    elif kind == 'bullet':
        p = doc.add_paragraph()
        add_runs(p, payload, size=12)
    else:
        p = doc.add_paragraph()
        add_runs(p, payload, size=12)
        if payload.startswith('**Table 2.**'):
            # Tables 2 and 3 share the final page rather than leaving Table 3 alone
            p.paragraph_format.page_break_before = True
        first_block = False
    prev_kind = kind

# --------------------------------------------------------------------------
# 7. Post-pass: reference hanging indent
# --------------------------------------------------------------------------
ref_start = None
ref_end = None
paras = doc.paragraphs
for idx, p in enumerate(paras):
    if p.style.name == 'Heading 1' and p.text.strip() == 'References':
        ref_start = idx
    elif ref_start is not None and p.style.name == 'Heading 1' and p.text.strip() == 'Figure Legends':
        ref_end = idx
        break
for p in paras[ref_start + 1:ref_end]:
    if re.match(r'^\d+\.\s', p.text):
        p.paragraph_format.left_indent = Cm(0.5)
        p.paragraph_format.first_line_indent = Cm(-0.5)

# --------------------------------------------------------------------------
# 8. Remove theme-font references so nothing resolves to Calibri
# --------------------------------------------------------------------------
W_NS = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
for rf in doc.styles.element.iter(W_NS + 'rFonts'):
    for attr in ('asciiTheme', 'hAnsiTheme', 'eastAsiaTheme', 'cstheme'):
        if rf.get(W_NS + attr) is not None:
            del rf.attrib[W_NS + attr]
    for attr in ('ascii', 'hAnsi', 'cs', 'eastAsia'):
        rf.set(W_NS + attr, FONT)
for rf in doc.element.iter(W_NS + 'rFonts'):
    for attr in ('asciiTheme', 'hAnsiTheme', 'eastAsiaTheme', 'cstheme'):
        if rf.get(W_NS + attr) is not None:
            del rf.attrib[W_NS + attr]

doc.save(OUT)
print('saved:', OUT)
print('blocks:', len(blocks))
print('paragraphs:', len(doc.paragraphs))
print('tables:', len(doc.tables))
