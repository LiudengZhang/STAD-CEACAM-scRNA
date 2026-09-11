#!/usr/bin/env python3
"""
Assemble the supplementary figures S7-S9.

Two things govern this assembler, and they are the reason it replaced the one
that declared a millimetre box per panel.

1. PANELS ARE PLACED AT THEIR NATURAL SIZE. The earlier assembler declared a
   box per panel and let `preserveAspectRatio='xMidYMid meet'` fit the panel
   into it. Across the supplementary set that fit ran from 0.16 to 0.45, so the
   same 7 pt of type printed anywhere between 1.1 and 3.2 pt depending only on
   which box it landed in. Here every panel is drawn at the size it prints at
   and placed at that size, and `check_scale()` refuses to assemble if any
   placement is not 1:1. Type set at 7 pt prints at 7 pt.

2. THE PANEL LETTERS COME FROM PROVENANCE.CSV, NOT FROM THE DIRECTORY NAME.
   For the supplementary figures the two agree, but a lookup that is only
   correct by luck is not a lookup - for the main figures they do not agree,
   which is why the rule exists. The map is built from the manifest and the
   assembler stops if a panel directory is not in it.

Layout is expressed as rows of panel directories. Row heights and the page
height follow from the panels themselves, so a panel that grew when its type
was set correctly cannot silently overlap its neighbour.

    python assemble_new_supplementaries.py
"""

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "00_Config"))

import panel_style_cns as style                       # noqa: E402
from shared.svg_assembler import VectorAssembler      # noqa: E402
from lxml import etree                                # noqa: E402

OUT = HERE / "_assembled"
PROVENANCE = ROOT / "03_Final_Panels" / "PROVENANCE.csv"
# Rows for panels that are new in this pass are held here until they are merged
# into the manifest. Same columns; a row whose build_path is "delete" names a
# retired panel and is not a lookup.
PROVENANCE_DELTA = ROOT / "03_Final_Panels" / "PROVENANCE_delta_S.csv"

PAGE_W = 183.0                  # double-column width, mm
MARGIN = 6.0                    # outer margin, mm
GUTTER = 5.0                    # between panels in a row
ROW_GAP = 6.0                   # between rows
TOP = 5.0
BOTTOM = 5.0
# A panel drawn 1:1 uses its whole box, so there is no white margin inside it
# for a panel letter to sit in: dropped at the panel's own top-left corner the
# letter lands on top of the y axis label. Each row therefore reserves a band
# above it for its letters.
LETTER_BAND = 5.0               # mm reserved above every row for its letters
LETTER_PT = 10.0                # panel letters, matching the main figures
LETTER_DX = 0.0                 # letter offset from the panel's left edge
LETTER_DY = 3.7                 # letter baseline below the top of the band

# The figure title is the name the caption gives it; it is not drawn on the
# page, which is how this supplementary set has always been set.
FIGURES = {
    "S7_CEACAM5_vs_CEACAM6": dict(
        title="CEACAM5 versus CEACAM6",
        rows=[["S7_A"], ["S7_B"], ["S7_C"]],
    ),
    "S8_Spatial_Confounders": dict(
        title="Spatial confounders and bulk validation",
        rows=[["S8_A"], ["S8_B"]],
    ),
    "S9_MoMac_Identity_NFkB": dict(
        title="MoMac identity and NF-kB specificity",
        rows=[["S9_A"], ["S9_B"], ["S9_C", "S9_D"], ["S9_E"]],
    ),
}

MM = 25.4 / 72.0                # points -> mm


# Selectable on the command line so one figure can be assembled while the rest
# of the panels are still being drawn.
def _selected(argv):
    want = [a for a in argv if not a.startswith("-")]
    return ({k: v for k, v in FIGURES.items() if k in want or
             any(k.startswith(w) for w in want)} if want else FIGURES)


class PanelAssembler(VectorAssembler):
    """VectorAssembler with the figure set's own type system for its own text.

    The base assembler hardcodes a font family and draws panel letters at
    `font_weight='normal'`. Overridden here rather than edited in place,
    because `svg_assembler.py` is shared with the main figures.
    """

    FONT_STACK = None            # filled in at construction

    def __init__(self, *a, **kw):
        family, _ = style.letter_font()
        # One family, named exactly. A list here would let the renderer
        # substitute, which is how panel letters end up in a different face
        # from the panels they label.
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


def _rows_of(path):
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def printed_letters():
    """{panel path within its tree -> printed panel letter}, from the manifest.

    Never inferred from the directory name. The pending delta is read after the
    manifest and overrides it, so a panel added in this pass can be looked up
    before the delta has been merged; a delta row marked for deletion names a
    retired panel and is skipped.
    """
    out = {}
    for path in (PROVENANCE, PROVENANCE_DELTA):
        for row in _rows_of(path):
            if row.get("build_path", "").strip() == "delete":
                continue
            src = row["source_dir"].strip()
            if not src.startswith(("Supplementary_New/",
                                   "Supplementary_Fixes/")):
                continue
            key = src.split("/", 1)[1]
            out[key] = row["printed_panel"].strip()
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
            f"no panel SVG in {d}\n"
            f"    Run the driver that produces {panel} before assembling. "
            f"Nothing is carried over from an earlier figure.")
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
            svg = panel_svg(figure, panel)
            w, h = natural_size_mm(svg)
            sizes.append((panel, svg, w, h))
        total_w = sum(s[2] for s in sizes) + GUTTER * (len(sizes) - 1)
        # A width read back out of a viewBox in points carries a rounding
        # error of a few parts in 10^15, so a row that exactly fills the page
        # measures a hair wider than it. The tolerance is the same 0.05 mm
        # check_scale() calls 1:1.
        live_w = PAGE_W - 2 * MARGIN
        if total_w > live_w + 0.05:
            raise ValueError(
                f"{figure} row {[s[0] for s in sizes]} is {total_w:.1f} mm "
                f"wide, wider than the {live_w:.1f} mm the page has between "
                f"its margins. Redraw a panel narrower - do not scale it down, "
                f"which is the fault this assembler exists to remove.")
        # Rows are left-aligned on the margin rather than centred, so every
        # panel letter on the page sits in the same column.
        x = MARGIN
        for panel, svg, w, h in sizes:
            key = f"{figure}/{panel}"
            if key not in letters:
                raise KeyError(
                    f"{key} is not in PROVENANCE.csv or the pending delta. "
                    f"The printed letter must be looked up, never inferred "
                    f"from the directory name.")
            placements.append(dict(panel=panel, letter=letters[key], svg=svg,
                                   x=x, y=y, w=w, h=h, letter_y=row_top))
            x += w + GUTTER
        y += max(s[3] for s in sizes) + ROW_GAP
    return placements, y - ROW_GAP + BOTTOM


# What a supplementary page can occupy once the journal has placed it: A4 with
# 20 mm margins. This is the number that decides whether setting type at 7 pt on
# the panel really puts 7 pt in front of the reader, because a page larger than
# this is scaled to fit and the type scales with it. Reported, never enforced:
# these are distributed as PDFs and are not reproduced at a fixed page size, and
# splitting a page to reach it would renumber the supplementary figures, which
# reaches the main text and the response letter.
PRINTABLE_W, PRINTABLE_H = 170.0, 247.0


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
            "panels would be rescaled at assembly, which is the fault this "
            "assembler exists to remove:\n  " + "\n  ".join(bad))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Must come before letter_font(): the resolved family depends on it.
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    family, weight = style.letter_font()
    print(f"page {PAGE_W:.1f} mm wide; panel letters "
          f"{LETTER_PT:g} pt {weight} {family}\n")
    written = []
    for figure, spec in _selected(sys.argv[1:]).items():
        placements, page_h = layout(figure, spec["rows"])
        check_scale(placements)
        asm = PanelAssembler(PAGE_W, page_h, title=None)
        print(f"{figure}   {PAGE_W:.1f} x {page_h:.1f} mm   {spec['title']}")
        for p in placements:
            asm.place_panel(None, p["svg"], p["x"], p["y"], p["w"], p["h"])
            asm.add_text(p["letter"], p["x"] + LETTER_DX,
                         p["letter_y"] + LETTER_DY,
                         font_size=LETTER_PT, font_weight=weight)
            print(f"    {p['letter']}  {p['panel']:<7} "
                  f"{p['w']:6.1f} x {p['h']:5.1f} mm  at "
                  f"({p['x']:5.1f}, {p['y']:5.1f})  scale 1.000")
        scale, smallest = fit_report(PAGE_W, page_h)
        print(f"    page fits {PRINTABLE_W:.0f} x {PRINTABLE_H:.0f} mm at "
              f"{scale:.3f}; smallest type then prints at {smallest:.2f} pt")
        asm.save(figure, output_dir=OUT)
        written.append((figure, page_h, scale, smallest))
    print(f"\nAssembled figures written to {OUT}")
    print(f"\n{'figure':<28}{'page mm':>14}{'journal fit':>13}"
          f"{'smallest pt':>13}")
    worst = 99.0
    for figure, page_h, scale, smallest in written:
        worst = min(worst, smallest)
        print(f"{figure:<28}{PAGE_W:6.1f} x {page_h:5.1f}{scale:13.3f}"
              f"{smallest:13.2f}")
    print(f"\nAt 1:1 - which is how these are distributed - the smallest body "
          f"type on every page above is {style.tick_pt():g} pt.")
    print(f"Were a page instead scaled into {PRINTABLE_W:.0f} x "
          f"{PRINTABLE_H:.0f} mm, the smallest would become {worst:.2f} pt.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
