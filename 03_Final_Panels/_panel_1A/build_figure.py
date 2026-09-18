# -*- coding: utf-8 -*-
"""Figure 1A - study overview, rebuilt as a publication-ready vector graphic.

Canvas 510.24 x 196 pt = 180 x 69.2 mm (full journal page width).
Text is kept as live <text> (Arial) so it stays editable/searchable; every
pictogram is an embedded <symbol> from an openly licensed scientific icon set:
  * Health Icons (MIT, https://healthicons.org)
  * Font Awesome Free 6 (icons: CC BY 4.0, https://fontawesome.com/license/free)
See ICON_CREDITS.md for the per-icon attribution line to paste in the legend.
"""
import math, os, re, sys

W, H = 510.24, 196.0

# ---------------------------------------------------------------------------
# Placed at the size it prints at, since 2026-09-15.
#
# The graphic was built 180 mm wide and the submitted Figure 1 page was 254 mm
# wide, so at the journal's column width its 6-8 pt type printed at ~5 pt.
# Figure 1 is now assembled like Figures 2-5, into a measured slot on the
# 171.10 mm page (03_Final_Panels/panel_rects_v2.csv, printed panel A), and
# this file is written at that slot's width. The drawing coordinates stay in
# the original 510.24 x 196 user units; only the declared size and the type
# change: every <text> is floored so that the smallest prints at MIN_PRINT_PT.
# The panel letter is no longer drawn here - the assembler draws every letter
# from letter_spec_v2.csv - so the "A" and its indent are gone.
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "00_Config"))
import slots                                    # noqa: E402
import panel_style_cns as _style                # noqa: E402
_style.apply()
#: ONE family, the one every other panel and the page's letters are set in,
#: resolved by panel_style_cns rather than a CSS stack: the assembler drops
#: a panel's root attributes when it nests the SVG, and cairo took the
#: generic fallback for "Arial, Helvetica, sans-serif" - the headings ran
#: out of their boxes and the italic n vanished on the assembled page
#: (2026-09-15). Set on a <g> round the drawing, which survives nesting.
F = _style.letter_font()[0]
SLOT_W_MM, SLOT_H_MM = slots.size_mm(1, "A")
if abs(SLOT_W_MM / SLOT_H_MM - W / H) > 0.01:
    sys.exit(f"the Figure 1 A slot ({SLOT_W_MM} x {SLOT_H_MM} mm) is not the "
             f"graphic's aspect ({W:.2f}:{H:.2f}); fix the slot in "
             f"build_grid_v2.REPAGED, do not letterbox the drawing")
PT_PER_MM = 72.0 / 25.4
UU_PER_PT = W / (SLOT_W_MM * PT_PER_MM)         # user units per printed point
MIN_PRINT_PT = 6.05                             # 6.0 exactly rounds to 5.99 in the PDF
MIN_UU = MIN_PRINT_PT * UU_PER_PT
INK = "#1a1a1a"
BLU_BD, BLU_HD, BLU_IC, BLU_TX = "#b7cbdd", "#eef3f8", "#8fb0c8", "#1b6b94"
PRE_F, PRE_S = "#cfe1ef", "#8fb2cd"
POR_F, POR_S = "#f9dcdc", "#dda9a9"
PNR_F, PNR_S = "#e2d0ea", "#b795c6"
SHIELD, CYCLE = "#3e627b", "#9d6ab5"
GRN_BD, GRN_HD, GRN_IC, TEAL_D = "#c6d7d1", "#f2f7f5", "#88ab9f", "#2b7f7a"
PER, ARROW = "#9db6ab", "#1c524c"
TAN_BD, TAN_HD, TAN_IC = "#e2c69d", "#fdf8f2", "#c9a878"
SUB_BD, SUB_F, BROWN, ORANGE = "#e7cca7", "#fffdf9", "#8a4a33", "#c2762e"
STO, LIV, OVA, LN, PBMC = "#d98d8d", "#b8705e", "#cb8a9c", "#7fa36a", "#a63a35"

# ---------- icon symbols from the licensed source files ----------
ICONS = {
    "brain":  "hi_body_neurology.svg",
    "people": "hi_people_people.svg",
    "stomach": "hi_body_stomach.svg",
    "liver":  "hi_body_liver.svg",
    "ovary":  "hi_body_female-reproductive_system.svg",
    "ln":     "hi_body_lymph-nodes.svg",
    "dna":    "hi_body_dna.svg",
    "cells":  "hi_body_cell-nuclei.svg",
    "bloodtube": "hi_devices_medical-sample.svg",
    "scope":  "hi_devices_microscope-with_specimen.svg",
    "vial":   "fa_vial.svg",
    "man":    "fa_person.svg",
    "woman":  "fa_person-dress.svg",
    "shield": "fa_shield.svg",
    "lock":   "fa_lock.svg",
    "cycle":  "fa_arrows-rotate.svg",
    "grid":   "fa_table-cells.svg",
    "check":  "fa_clipboard-check.svg",
}

def load_symbol(sid, fname):
    s = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "icons_src", fname), encoding="utf-8").read()
    vb = re.search(r'viewBox="([^"]+)"', s).group(1)
    inner = s[s.index(">", s.index("<svg")) + 1: s.rindex("</svg>")]
    inner = re.sub(r"<!--.*?-->", "", inner, flags=re.S)
    inner = re.sub(r"<title>.*?</title>", "", inner, flags=re.S)
    inner = re.sub(r"\s+", " ", inner).strip()
    return f'<symbol id="ic-{sid}" viewBox="{vb}">{inner}</symbol>'

o = []
def add(s): o.append(s)
def esc(t): return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
def T(x, y, s, size=6.5, w="normal", anchor="middle", fill=INK, raw=False):
    size = max(size, MIN_UU)                    # floor, in printed points
    b = f' font-weight="{w}"' if w != "normal" else ""
    body = s if raw else esc(s)
    add(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}"{b} text-anchor="{anchor}" fill="{fill}">{body}</text>')

# AACR: italic lowercase n. Set as three <text> runs laid out from the
# face's own advance widths (cnsfig.rich._advance_pt), not as a <tspan>:
# the assembler's cairo renderer printed the tspan run on top of the text
# before it and lost the rest (2026-09-15).
from cnsfig.rich import _advance_pt              # noqa: E402


def TN(x, y, prefix, v, suffix="", size=6.5, anchor="middle", fill=INK):
    """prefix + italic n + ' = v' + suffix, centred (or anchored) at x."""
    size = max(size, MIN_UU)
    parts = [(prefix, False), ("n", True), (f" = {v}{suffix}", False)]
    widths = [_advance_pt(t, size, it) if t else 0.0 for t, it in parts]
    total = sum(widths)
    x0 = x - total / 2 if anchor == "middle" else (x - total if anchor == "end" else x)
    for (t, it), w in zip(parts, widths):
        if not t:
            continue
        # Leading/trailing spaces carry their advance but not their ink;
        # the PDF drops a leading space, so the run is written stripped and
        # the space's width is kept in x0.
        lead = len(t) - len(t.lstrip(" "))
        if lead:
            x0 += _advance_pt(" " * lead, size, it)
        face = ' font-style="italic"' if it else ''
        add(f'<text x="{x0:.2f}" y="{y:.2f}" font-size="{size}"{face} '
            f'text-anchor="start" fill="{fill}">{esc(t.strip(" "))}</text>')
        x0 += w - (_advance_pt(" " * lead, size, it) if lead else 0.0)
def rrect(x, y, w, h, r, fill, stroke="none", sw=0.6):
    add(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{r}" ry="{r}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
def line(x1, y1, x2, y2, stroke=INK, sw=0.7, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" '
        f'stroke-width="{sw}" stroke-linecap="round"{d}/>')
def circ(cx, cy, r, fill, stroke="none", sw=0.6):
    add(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
def use(sid, x, y, w, h, color):
    add(f'<use href="#ic-{sid}" xlink:href="#ic-{sid}" x="{x:.2f}" y="{y:.2f}" '
        f'width="{w:.2f}" height="{h:.2f}" fill="{color}" color="{color}"/>')
def g_open(**kw): add("<g " + " ".join(f'{k.replace("_","-")}="{v}"' for k, v in kw.items()) + ">")
def g_close(): add("</g>")
def head(x, y, ang, L=3.0, wd=2.0, fill=INK):
    bx, by = x - L*math.cos(ang), y - L*math.sin(ang)
    px, py = -math.sin(ang)*wd/2, math.cos(ang)*wd/2
    add(f'<path d="M{x:.2f},{y:.2f} L{bx+px:.2f},{by+py:.2f} L{bx-px:.2f},{by-py:.2f} Z" fill="{fill}"/>')
def arrow(x1, y1, x2, y2, stroke=INK, sw=0.7, dash=None, both=False, L=3.0, wd=2.0):
    a = math.atan2(y2-y1, x2-x1)
    sx, sy = (x1 + 0.8*L*math.cos(a), y1 + 0.8*L*math.sin(a)) if both else (x1, y1)
    line(sx, sy, x2 - 0.8*L*math.cos(a), y2 - 0.8*L*math.sin(a), stroke, sw, dash)
    head(x2, y2, a, L, wd, stroke)
    if both: head(x1, y1, a + math.pi, L, wd, stroke)

# ================= document =================
add(f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
    f'width="{SLOT_W_MM:.2f}mm" height="{SLOT_H_MM:.2f}mm" viewBox="0 0 {W:.2f} {H:.2f}" font-family="{F}">')
add('<title>Figure 1A. Study overview: conceptual framework, design, sample collection, validation cohorts</title>')
add('<!-- Pictograms: Health Icons (MIT, healthicons.org) and Font Awesome Free 6 '
    '(icons CC BY 4.0, fontawesome.com/license/free). See ICON_CREDITS.md. -->')
add("<defs>")
for sid, fn in ICONS.items():
    add(load_symbol(sid, fn))
add("</defs>")
add(f'<rect width="{W:.2f}" height="{H:.2f}" fill="#ffffff"/>')
add(f'<g font-family="{F}">')                    # closed before </svg>

# The panel letter is the assembler's (letter_spec_v2.csv), not this file's.
line(5, 20.5, 498, 20.5, INK, 1.1)
head(506, 20.5, 0, 7.0, 5.2, INK)

PX, PW, PY0, PY1 = [5, 130, 255, 380], 120.0, 25.0, 192.0
def shell(i, bd, hd, ic, title, glyph):
    x = PX[i]
    add(f'<clipPath id="cp{i}"><rect x="{x}" y="{PY0}" width="{PW}" height="{PY1-PY0}" rx="4" ry="4"/></clipPath>')
    rrect(x, PY0, PW, PY1-PY0, 4, "#ffffff", bd, 0.8)
    add(f'<g clip-path="url(#cp{i})">'); rrect(x, PY0, PW, 22, 0, hd); g_close()
    line(x, PY0+22, x+PW, PY0+22, bd, 0.6)
    circ(x+13, 36, 7.5, ic)
    use(glyph, x+13-5.4, 36-5.4, 10.8, 10.8, "#ffffff")
    T(x+24, 39, title, 8, "bold", "start")

# ---- Panel 1: Conceptual Framework ----
# The two arms are cross-sectional response contrasts at each treatment
# timepoint. Only one patient contributed both timepoints, so the schematic
# must not imply a paired pre-to-post trajectory.
g_open(id="panel-conceptual-framework")
shell(0, BLU_BD, BLU_HD, BLU_IC, "Conceptual Framework", "brain")
cxL, cxR = 33.0, 97.0
T(cxL, 58, "Responders", 6.8, "bold", fill=BLU_TX)
T(cxR, 58, "Nonresponders", 6.8, "bold", fill=BLU_TX)
for cx, lab in ((cxL, "Pre-R"), (cxR, "Pre-NR")):
    rrect(cx-19, 64, 38, 14, 2.5, PRE_F, PRE_S, 0.7); T(cx, 73.4, lab, 7, "bold")
for cx, lab, f_, s_ in ((cxL, "Post-R", POR_F, POR_S), (cxR, "Post-NR", PNR_F, PNR_S)):
    rrect(cx-19, 122, 38, 14, 2.5, f_, s_, 0.7); T(cx, 131.4, lab, 7, "bold")

arrow(53, 71, 77, 71, SHIELD, 1.0, both=True, L=3.0, wd=2.2)
line(65, 79, 65, 84, SHIELD, 0.7, dash="0.6 1.7")
use("shield", 59, 85, 12, 12, SHIELD)
use("lock", 61.9, 88.6, 6.2, 7.1, "#ffffff")
T(65, 105, "Pre-treatment non-response", 6.4, "bold", fill=SHIELD)
T(65, 112.5, "Pre-NR vs Pre-R", 6.0)

arrow(53, 129, 77, 129, CYCLE, 1.0, both=True, L=3.0, wd=2.2)
T(65, 153, "Post-treatment non-response", 6.4, "bold", fill=CYCLE)
T(65, 162, "Post-NR vs Post-R", 6.0)
g_close()

# ---- Panel 2: Study Design ----
g_open(id="panel-study-design")
shell(1, GRN_BD, GRN_HD, GRN_IC, "Study Design", "people")
T(190, 88, "Anti\u2013PD-1 +", 6.6, "bold"); T(190, 96.5, "chemotherapy", 6.6, "bold")
for cx in (156.0, 222.0):
    use("man", cx-19.5, 102, 21.25, 34, PER)
    use("woman", cx-1.75, 102, 21.25, 34, PER)
line(179, 120, 197.5, 120, ARROW, 2.4); head(201, 120, 0, 5.0, 6.4, ARROW)
T(156, 156, "Pre-Treatment", 6.6, "bold");  T(156, 165.5, "(N = 24)", 6.6, "bold")
T(222, 156, "Post-Treatment", 6.6, "bold"); T(222, 165.5, "(N = 12)", 6.6, "bold")
g_close()

# ---- Panel 3: Sample Collection ----
g_open(id="panel-sample-collection")
shell(2, TAN_BD, TAN_HD, TAN_IC, "Sample Collection", "vial")
rrect(259, 52, 52, 118, 3, SUB_F, SUB_BD, 0.7)
T(285, 65, "Primary tumor", 7.2, "bold", fill=BROWN)
use("stomach", 260, 72, 50, 66, STO)
TN(285, 155, "Stomach, ", 32, size=6.1)
rrect(315, 52, 56, 56, 3, SUB_F, SUB_BD, 0.7)
T(343, 64, "Metastases", 7, "bold", fill=ORANGE)
line(343, 69, 343, 104, ORANGE, 0.7, dash="1.2 1.6")
use("liver", 321, 69, 18, 18, LIV); use("ovary", 349, 69, 18, 18, OVA)
T(329, 95, "Liver", 6.4);  TN(329, 103, "", 12, size=6.4)
T(357, 95, "Ovary", 6.4);  TN(357, 103, "", 3, size=6.4)
rrect(315, 112, 56, 58, 3, SUB_F, SUB_BD, 0.7)
T(343, 124, "Additional", 7, "bold", fill=ORANGE)
line(343, 129, 343, 165, ORANGE, 0.7, dash="1.2 1.6")
use("ln", 321, 130, 18, 18, LN); use("bloodtube", 349, 130, 18, 18, PBMC)
T(329, 156, "LN", 6.4);   TN(329, 164, "", 12, size=6.4)
T(357, 156, "PBMC", 6.4); TN(357, 164, "", 11, size=6.4)
T(315, 180, "All samples were profiled", 6.5)
T(315, 188, "by scRNA-seq", 6.5)
g_close()

# ---- Panel 4: Validation Cohorts ----
g_open(id="panel-validation-cohorts")
shell(3, GRN_BD, GRN_HD, GRN_IC, "Validation Cohorts", "check")
rows = [(69,  "scope", "In-house Experiment", ["\u2022 H&E staining", "\u2022 Immunohistochemistry"]),
        (103, "dna",   "Bulk RNA-seq",        ["PRJEB25780", "TCGA-STAD"]),
        (137, "cells", "External scRNA-seq",  ["GSE239676", "GSE183904"]),
        (171, "grid",  "Spatial transcriptomics", ["GSE251950"])]
for k, (cy, gl, title, lines) in enumerate(rows):
    use(gl, 388, cy-10, 20, 20, TEAL_D)
    dy = -6.5 if len(lines) == 2 else -4.5
    T(412, cy+dy, title, 6.5, "bold", "start", TEAL_D)
    for j, ln in enumerate(lines):
        T(412, cy+dy+9.5+j*8.2, ln, 6.3, anchor="start")
    if k < 3:
        line(386, cy+17, 494, cy+17, GRN_BD, 0.7, dash="1.2 1.6")
g_close()
add("</g>")                                       # the font-family group
add("</svg>")
HERE = os.path.dirname(os.path.abspath(__file__))
open(os.path.join(HERE, "figure1A.svg"), "w", encoding="utf-8").write("\n".join(o) + "\n")
print(f"figure1A.svg: {SLOT_W_MM:.2f} x {SLOT_H_MM:.2f} mm; type floored at "
      f"{MIN_PRINT_PT} pt printed ({MIN_UU:.2f} user units)")

credits = """# Icon credits - Figure 1A

All pictograms are openly licensed vector icons, embedded as SVG `<symbol>`s
(no bitmaps, no BioRender-derived artwork).

| Element | Icon | Source | License |
|---|---|---|---|
| Conceptual Framework (header) | neurology / brain | Health Icons | MIT |
| Study Design, Validation Cohorts (headers) | people | Health Icons | MIT |
| Sample Collection (header) | vial | Font Awesome Free 6 | CC BY 4.0 |
| Pre-treatment non-response | shield + lock | Font Awesome Free 6 | CC BY 4.0 |
| Patients (pre/post) | person, person-dress | Font Awesome Free 6 | CC BY 4.0 |
| Primary tumour | stomach | Health Icons | MIT |
| Liver metastasis | liver | Health Icons | MIT |
| Ovarian metastasis | female reproductive system | Health Icons | MIT |
| Lymph node | lymph-nodes | Health Icons | MIT |
| PBMC | medical sample (blood tube) | Health Icons | MIT |
| In-house experiment | microscope with specimen | Health Icons | MIT |
| Bulk RNA-seq | dna | Health Icons | MIT |
| External scRNA-seq | cell-nuclei | Health Icons | MIT |
| Spatial transcriptomics | table-cells | Font Awesome Free 6 | CC BY 4.0 |

Suggested figure-legend line (CC BY 4.0 requires attribution; MIT does not,
but the courtesy credit is kept):

> Icons from Health Icons (healthicons.org, MIT) and Font Awesome Free 6.7
> (fontawesome.com, icons licensed CC BY 4.0); icons were recoloured and
> rescaled for this figure.

Health Icons: https://github.com/resolvetosavelives/healthicons (MIT)
Font Awesome Free: https://fontawesome.com/license/free (icons CC BY 4.0)
"""
open(os.path.join(HERE, "ICON_CREDITS.md"), "w", encoding="utf-8").write(credits)
print("figure1A.svg", os.path.getsize("figure1A.svg"), "bytes")
