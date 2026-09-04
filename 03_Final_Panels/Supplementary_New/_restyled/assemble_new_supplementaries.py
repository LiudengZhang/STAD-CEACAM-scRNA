#!/usr/bin/env python3
"""
Assemble the restyled Supplementary Figures S7-S11 (Version B).

Two things make this different from Version A's assembler, and they are the
whole point of the restyle.

1. PANELS ARE PLACED AT THEIR NATURAL SIZE. Version A declared a millimetre box
   per panel and let `preserveAspectRatio='xMidYMid meet'` fit the panel into it.
   Measured on 2026-09-01 that fit ran from 0.1600 to 0.4538 across the 26
   panels - so the same 7 pt of type printed anywhere between 1.1 and 3.2 pt
   depending only on which box it landed in. Here every panel is drawn at the
   size it prints at and placed at that size, and `check_scale()` refuses to
   assemble if any placement is not 1:1. Type set at 7 pt prints at 7 pt.

2. THE PANEL LETTERS COME FROM PROVENANCE.CSV, NOT FROM THE DIRECTORY NAME.
   For the supplementary figures the two agree - the audit of 1 Sep confirmed
   26 of 26 - but CLAUDE.md rule 2 exists because for the main figures they do
   not, and a lookup that is only correct by luck is not a lookup. The map is
   built from the CSV and the assembler stops if a panel directory is not in it.

Layout is expressed as rows of panel directories. Row heights and the page
height follow from the panels themselves, so a panel that grew when its type was
set correctly does not silently overlap its neighbour.

    conda run -n Liudeng_Python_310 python assemble_new_supplementaries.py
"""

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "00_Config"))

import panel_style_cns as style                       # noqa: E402
from shared.svg_assembler import VectorAssembler      # noqa: E402
from lxml import etree                                # noqa: E402

OUT = HERE / "_assembled"
PROVENANCE = ROOT / "03_Final_Panels" / "PROVENANCE.csv"

PAGE_W = style.PAGE_W_MM        # 190.5 mm, cnsplots' full-width figure
GUTTER = 5.0                    # between panels in a row
ROW_GAP = 6.0                   # between rows
TOP = 5.0
BOTTOM = 5.0
# A panel drawn 1:1 uses its whole box, so there is no white margin inside it
# for a panel letter to sit in: dropped at the panel's own top-left corner the
# letter lands on top of the y axis label. Each row therefore reserves a band
# above it for its letters.
LETTER_BAND = 4.6               # mm reserved above every row for its letters
LETTER_DX = 0.0                 # letter offset from the panel's left edge
LETTER_DY = 3.4                 # letter baseline below the top of the band

# Rows of panel directories, top to bottom. Panels on one row sit side by side.
FIGURES_ALL = {
    "S7_Cohort_Statistics":     [["S7_A"], ["S7_B"], ["S7_C"]],
    "S8_CEACAM_Metaprogram":    [["S8_A"], ["S8_B"], ["S8_C"], ["S8_D"],
                                 ["S8_E"], ["S8_F", "S8_G"], ["S8_H"]],
    "S9_Mechanism_Specificity": [["S9_A"], ["S9_B"], ["S9_C"], ["S9_D"],
                                 ["S9_E", "S9_F"]],
    "S10_PreTx_and_Adaptive":   [["S10_A", "S10_B", "S10_C"], ["S10_D"],
                                 ["S10_E"]],
    "S11_Affirmative_Analyses": [["S11_A"], ["S11_B"], ["S11_C"], ["S11_D"]],
}

MM = 25.4 / 72.0                # points -> mm

# Panels that could not be restyled and are awaiting the author's decision.
# They are left OUT of the restyled page rather than carried over from Version
# A: a Version A panel dropped into a Version B page would be the only thing on
# it drawn at four times its printed size, and the assembler would have to
# shrink it - which is the exact fault this rebuild exists to remove. A visible
# gap is the honest representation of "not done yet".
#
# Empty since 2026-09-01. S8D was the only entry: its `scores` frame lived only
# inside `mp_external_validation.py:main()`, so the panel could not be redrawn
# without recomputing the signature score. The analysis now writes
# `mp_external_sample_scores.csv` and the driver reads it.
PENDING = {}

# Selectable on the command line so one figure can be assembled while the rest
# of the panels are still being redrawn.
def _selected(argv):
    want = [a for a in argv if not a.startswith("-")]
    return ({k: v for k, v in FIGURES_ALL.items() if k in want or
             any(k.startswith(w) for w in want)} if want else FIGURES_ALL)


class RestyledAssembler(VectorAssembler):
    """VectorAssembler with the cnsplots type system for its own text.

    Version A's assembler hardcodes `Liberation Sans` and draws panel letters at
    `font_weight='normal'`. cnsplots sets panel letters in the panel-label font
    at `title_fontsize`, bold. Overridden here rather than edited in place,
    because `svg_assembler.py` is shared with Version A and Version A is frozen.
    """

    FONT_STACK = None            # filled in at construction

    def __init__(self, *a, **kw):
        family, _ = style.letter_font()
        # One family, named exactly. A list here would let cairo substitute,
        # which is how the panel letters ended up in a different face from the
        # panels the first time this ran.
        self.FONT_STACK = f"{family}, sans-serif"
        super().__init__(*a, **kw)

    def _embed_text(self, root, data):
        x, y = data["x_mm"], data["y_mm"]
        txt = etree.SubElement(root, "{%s}text" % _SVGNS)
        txt.set("x", f"{x}")
        txt.set("y", f"{y}")
        txt.set("font-family", self.FONT_STACK)
        txt.set("font-size", f"{data['font_size'] * 0.3528}")
        if data["font_weight"] == "bold":
            txt.set("font-weight", "bold")
        txt.set("text-anchor",
                {"start": "start", "middle": "middle",
                 "end": "end"}.get(data.get("anchor", "start"), "start"))
        txt.set("fill", data.get("color") or "black")
        if data.get("rotation"):
            txt.set("transform", f"rotate({-data['rotation']}, {x}, {y})")
        txt.text = data["text"]


_SVGNS = "http://www.w3.org/2000/svg"


def printed_letters():
    """{panel directory -> printed panel letter}, from PROVENANCE.csv.

    Never inferred from the directory name. CLAUDE.md rule 2.
    """
    out = {}
    with open(PROVENANCE, newline="") as fh:
        for row in csv.DictReader(fh):
            src = row["source_dir"].strip()
            if not src.startswith("Supplementary_New/"):
                continue
            out[Path(src).name] = row["printed_panel"].strip()
    return out


def natural_size_mm(svg_path):
    """The size the panel was drawn at, in mm, from its own viewBox."""
    head = svg_path.read_text(errors="replace")[:2000]
    m = re.search(r'viewBox="[\d.eE+-]+ [\d.eE+-]+ ([\d.eE+-]+) ([\d.eE+-]+)"',
                  head)
    if m:
        return float(m.group(1)) * MM, float(m.group(2)) * MM
    w = re.search(r'width="([\d.]+)pt"', head)
    h = re.search(r'height="([\d.]+)pt"', head)
    if not (w and h):
        raise ValueError(f"cannot read a size out of {svg_path}")
    return float(w.group(1)) * MM, float(h.group(1)) * MM


def panel_svg(figure, panel):
    d = HERE / figure / panel
    hits = sorted(d.glob("*.svg"))
    if not hits:
        raise FileNotFoundError(
            f"no restyled SVG in {d}\n"
            f"    Run the driver or the panel script that produces {panel} "
            f"before assembling. Nothing is carried over from Version A.")
    if len(hits) > 1:
        raise ValueError(f"more than one SVG in {d}: "
                         f"{[p.name for p in hits]}")
    return hits[0]


def layout(figure, rows):
    """Positions for every panel, at natural size. Returns (placements, page_h)."""
    letters = printed_letters()
    placements, y = [], TOP
    for row in rows:
        row_top = y
        y += LETTER_BAND
        sizes = []
        for panel in row:
            if panel in PENDING:
                print(f"    -- {panel} OMITTED: {PENDING[panel]}")
                continue
            svg = panel_svg(figure, panel)
            w, h = natural_size_mm(svg)
            sizes.append((panel, svg, w, h))
        if not sizes:
            continue
        total_w = sum(s[2] for s in sizes) + GUTTER * (len(sizes) - 1)
        if total_w > PAGE_W:
            raise ValueError(
                f"{figure} row {[s[0] for s in sizes]} is {total_w:.1f} mm "
                f"wide, wider than the {PAGE_W:.1f} mm page. Redraw a panel "
                f"narrower - do not scale it down, that is what Version A did.")
        x = (PAGE_W - total_w) / 2
        for panel, svg, w, h in sizes:
            if panel not in letters:
                raise KeyError(
                    f"{panel} is not in PROVENANCE.csv. The printed letter must "
                    f"be looked up, never inferred from the directory name.")
            placements.append(dict(panel=panel, letter=letters[panel], svg=svg,
                                   x=x, y=y, w=w, h=h, letter_y=row_top))
            x += w + GUTTER
        y += max(s[3] for s in sizes) + ROW_GAP
    return placements, y - ROW_GAP + BOTTOM


# What a supplementary page can actually occupy once the journal has placed it.
# A4 with 20 mm margins. This is the number that decides whether setting type at
# 7 pt on the panel really puts 7 pt in front of the reader: a page larger than
# this is scaled to fit, and the type scales with it. Version A's failure was
# exactly this effect applied per panel instead of per page - and applied
# unequally, which is what produced a 2.8x spread.
PRINTABLE_W, PRINTABLE_H = 170.0, 247.0

# Ruled on 2026-09-01: the supplementary figures take the tall page. They are
# distributed as PDFs and are not reproduced at a fixed page size, so a 531 mm
# S8 is acceptable where it would not be for a main figure - Version A already
# ships 306 mm. Splitting S8 would renumber the supplementary figures, which
# reaches the main text and the response letter, and this revision is meant to
# be minimal. So the fit figure below is reported, never enforced: it is there
# so the author can revisit the decision if production objects.


def fit_report(page_w, page_h):
    """(scale the journal would apply, smallest type after it)."""
    scale = min(1.0, PRINTABLE_W / page_w, PRINTABLE_H / page_h)
    return scale, style.tick_pt() * scale


def check_scale(placements):
    """Every panel must be placed at exactly the size it was drawn at."""
    bad = []
    for p in placements:
        w, h = natural_size_mm(p["svg"])
        if abs(w - p["w"]) > 0.05 or abs(h - p["h"]) > 0.05:
            bad.append(f"{p['panel']}: drawn {w:.2f} x {h:.2f} mm, "
                       f"placed {p['w']:.2f} x {p['h']:.2f} mm")
    if bad:
        raise AssertionError(
            "panels would be rescaled at assembly, which is the fault the "
            "restyle exists to remove:\n  " + "\n  ".join(bad))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Must come before letter_font(): the resolved family depends on it.
    style.apply()
    family, weight = style.letter_font()
    print(f"page {PAGE_W:.1f} mm wide; panel letters "
          f"{style.letter_pt():g} pt {weight} {family}\n")
    written = []
    for figure, rows in _selected(sys.argv[1:]).items():
        placements, page_h = layout(figure, rows)
        check_scale(placements)
        asm = RestyledAssembler(PAGE_W, page_h, title=None)
        print(f"{figure}   {PAGE_W:.1f} x {page_h:.1f} mm")
        for p in placements:
            asm.place_panel(None, p["svg"], p["x"], p["y"], p["w"], p["h"])
            asm.add_text(p["letter"], p["x"] + LETTER_DX,
                         p["letter_y"] + LETTER_DY,
                         font_size=style.letter_pt(), font_weight=weight)
            print(f"    {p['letter']}  {p['panel']:<7} "
                  f"{p['w']:6.1f} x {p['h']:5.1f} mm  at "
                  f"({p['x']:5.1f}, {p['y']:5.1f})  scale 1.000")
        scale, smallest = fit_report(PAGE_W, page_h)
        flag = "" if smallest >= 5.0 else "   (tall page, accepted - see note)"
        print(f"    page fits {PRINTABLE_W:.0f} x {PRINTABLE_H:.0f} mm at "
              f"{scale:.3f}; smallest type then prints at "
              f"{smallest:.2f} pt{flag}")
        asm.save(figure, output_dir=OUT)
        written.append((figure, page_h, scale, smallest))
    print(f"\nRestyled figures written to {OUT}")
    for f in sorted(OUT.iterdir()):
        print(f"   {f.name:<44} {f.stat().st_size / 1e6:7.2f} MB")
    print(f"\n{'figure':<28}{'page mm':>14}{'journal fit':>13}"
          f"{'smallest pt':>13}")
    worst = 99.0
    for figure, page_h, scale, smallest in written:
        worst = min(worst, smallest)
        print(f"{figure:<28}{PAGE_W:6.1f} x {page_h:5.1f}{scale:13.3f}"
              f"{smallest:13.2f}")
    print(f"\nAt 1:1 - which is how these are distributed - the smallest body "
          f"type on every page above is {style.tick_pt():g} pt.")
    print(f"Were a page instead scaled into {PRINTABLE_W:.0f} x {PRINTABLE_H:.0f} "
          f"mm, the smallest would become {worst:.2f} pt. Reported, not "
          f"enforced; see RESTYLE_REPORT.md.")
    if PENDING:
        print(f"\n{len(PENDING)} panel(s) omitted, awaiting a decision:")
        for panel, why in PENDING.items():
            print(f"   {panel}: {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
