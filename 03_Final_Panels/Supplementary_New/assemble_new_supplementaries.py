#!/usr/bin/env python3
"""
Assemble the supplementary figures S1-S10.

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
import slots                                          # noqa: E402
from shared.svg_assembler import VectorAssembler      # noqa: E402
from lxml import etree                                # noqa: E402

OUT = HERE / "_assembled"
PROVENANCE = ROOT / "03_Final_Panels" / "PROVENANCE.csv"
# Rows for panels that are new in this pass are held here until they are merged
# into the manifest. Same columns; a row whose build_path is "delete" names a
# retired panel and is not a lookup.
PROVENANCE_DELTA = ROOT / "03_Final_Panels" / "PROVENANCE_delta_S.csv"

# THE PAGE IS THE MAIN FIGURES' PAGE WIDTH, since 2026-09-15. At 183 mm with
# 6 mm margins the journal scaled the page to its 170 mm printable width and
# the 6 pt type with it, to 5.6 pt; every figure is published at one width, so
# the supplementary pages are set at the main figures' 171.10 mm with the
# panels flush to the edges - the journal's own margins surround the page.
PAGE_W = 171.10                 # the main figures' page width, mm
MARGIN = 0.0                    # the page IS the printable area
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
# page, which is how this supplementary set has always been set. The one
# string the assembler sets besides the letters is the "continued" note at
# the top right of the second and later pages of a figure that has more than
# one (none since 2026-09-16; S1 for a day), at the body size, so a reader of the merged supplement knows the
# page belongs to the figure before it.
CONTINUED_NOTE = "Figure {n}, continued"
CONTINUED_PT = 7.0
CONTINUED_DY = 3.7              # baseline below the page top, in the TOP band

# The tallest a supplementary page may be: the main figures' page
# (00_Config/slots.page_mm), since the evening of 2026-09-15. AACR holds
# supplementary data to the main article's standard of presentation and its
# type to 6 pt at print size; a page taller than the main figures' would be
# scaled to fit and its 6 pt type with it. Enforced here and in
# 12_Figure_Refactor/sweep_pages.py; a figure that cannot fit takes a second
# page (`pages` in FIGURES) rather than a smaller type.
PAGE_H_MAX = slots.page_mm(2)[1]
FIGURES = {
    # S1-S6 since 2026-09-15 (evening; the author's fourth reading: "the type
    # in the supplementary figures still looks off"). Until then the six
    # submitted pages were carried over with their panels untouched, at
    # 2.4-4.9 pt body type; every panel is now redrawn at 1:1 into
    # Supplementary_New/S<n>_*/S<n>_<letter>/ and assembled here like S7-S9.
    #
    # TEN FIGURES SINCE 2026-09-16 (the author's fifth reading). The submitted
    # S1 held QC and annotation together: seven rows, 431 mm at 6 pt. On the
    # evening of 2026-09-15 it went onto two pages (`pages`, CONTINUED_NOTE);
    # the next morning the author wanted it on one page and, when it would
    # not fit ("B and C on one row" needs 188 mm of width for 70 rotated
    # sample labels at 6 pt), ruled to SPLIT it: S1 = QC (A-E), S2 = cell-type
    # annotation (the former F, G, H as A, B, C), and the former S2-S9 are
    # S3-S10. The directories were renamed with the printed numbers, so this
    # table's keys and the panel names are the printed names again (rule 2).
    # The two-page mechanism stays in the code, exercised by nothing today.
    "S1_QC_Annotation": dict(
        title="Quality control",
        # A 38, B 44, C 46, D/E 40 mm tall since 2026-09-16; ~216 mm.
        rows=[["S1_A"], ["S1_B"], ["S1_C"], ["S1_D", "S1_E"]],
    ),
    "S2_Cell_Annotation": dict(
        title="Cell-type annotation",
        # A and B share a row and one key: the two dot plots are scaled alike
        # (standard_scale='var', the same largest dot) and B's key at the
        # right of the row decodes both (_dotplot.draw(shared_key="S2 B")).
        rows=[["S2_A", "S2_B"], ["S2_C"]],
    ),
    "S3_CEACAM_Metaprogram_Validation": dict(
        title="CEACAM5/6 epithelial state and metaprogram validation",
        # A across the row and G beside H since 2026-09-16 ("G sits oddly and
        # breaks the letter order"): the letters read A-H row by row.
        rows=[["S3_A"], ["S3_B", "S3_C", "S3_D", "S3_E"], ["S3_F"],
              ["S3_G", "S3_H"]],
    ),
    "S4_CD8_TCells": dict(
        title="CD8+ T-cell states and adaptive immune composition",
        rows=[["S4_A"], ["S4_B", "S4_C", "S4_D"], ["S4_E"]],
    ),
    "S5_Spatial_Validation": dict(
        title="Spatial validation (GSE251950)",
        # F beside E since 2026-09-16, a 3 x 2 page.
        rows=[["S5_A", "S5_B"], ["S5_C", "S5_D"], ["S5_E", "S5_F"]],
    ),
    "S6_Immune_Modules": dict(
        title="Immune-module proportions by treatment phase and response",
        rows=[["S6_A", "S6_B", "S6_C"], ["S6_D", "S6_E"]],
    ),
    "S7_CD274_Remaining": dict(
        title="PD-L1 (CD274) expression in the remaining cell types",
        rows=[["S7_A", "S7_B", "S7_C", "S7_D", "S7_E"],
              ["S7_F", "S7_G", "S7_H", "S7_I"]],
    ),
    "S8_CEACAM5_vs_CEACAM6": dict(
        title="CEACAM5 versus CEACAM6",
        rows=[["S8_A"], ["S8_B"], ["S8_C"]],
    ),
    "S9_Spatial_Confounders": dict(
        title=("Density-stratified spatial analyses and independent bulk "
               "validation"),
        rows=[["S9_A", "S9_B"], ["S9_C"], ["S9_D", "S9_E"]],
    ),
    "S10_MoMac_Identity_NFkB": dict(
        title="MoMac identity and NF-kB specificity",
        # Until 2026-09-16 A and E shared the top row and E printed before B;
        # the author's fifth reading: "E appears before B". A takes the row
        # (161 mm, its labels get room), B the next, and C, D, E the last
        # (66 + 48 + 37 mm), so the letters read in order.
        rows=[["S10_A"], ["S10_B"], ["S10_C", "S10_D", "S10_E"]],
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
            if not src.startswith("Supplementary_New/"):
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
    page_h = y - ROW_GAP + BOTTOM
    if page_h > PAGE_H_MAX + 0.05:
        raise ValueError(
            f"{figure}: {page_h:.1f} mm tall, taller than the {PAGE_H_MAX:.2f} mm "
            f"main-figure page a supplementary page is held to. Give the "
            f"figure a second page (`pages` in FIGURES) - do not scale it.")
    return placements, page_h


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
        pages = spec.get("pages") or [spec["rows"]]
        n_fig = re.match(r"S\d+", figure).group(0)
        page_pdfs = []
        for k, rows in enumerate(pages, start=1):
            placements, page_h = layout(figure, rows)
            check_scale(placements)
            asm = PanelAssembler(PAGE_W, page_h, title=None)
            tag = f" (page {k} of {len(pages)})" if len(pages) > 1 else ""
            print(f"{figure}{tag}   {PAGE_W:.1f} x {page_h:.1f} mm   {spec['title']}")
            for p in placements:
                asm.place_panel(None, p["svg"], p["x"], p["y"], p["w"], p["h"])
                asm.add_text(p["letter"], p["x"] + LETTER_DX,
                             p["letter_y"] + LETTER_DY,
                             font_size=LETTER_PT, font_weight=weight)
                print(f"    {p['letter']}  {p['panel']:<7} "
                      f"{p['w']:6.1f} x {p['h']:5.1f} mm  at "
                      f"({p['x']:5.1f}, {p['y']:5.1f})  scale 1.000")
            if k > 1:
                asm.add_text(CONTINUED_NOTE.format(n=n_fig), PAGE_W - MARGIN,
                             CONTINUED_DY, font_size=CONTINUED_PT,
                             font_weight="normal", anchor="end")
            scale, smallest = fit_report(PAGE_W, page_h)
            print(f"    page fits {PRINTABLE_W:.0f} x {PRINTABLE_H:.0f} mm at "
                  f"{scale:.3f}; smallest type then prints at {smallest:.2f} pt")
            stem = figure if len(pages) == 1 else f"{figure}_p{k}"
            asm.save(stem, output_dir=OUT)
            page_pdfs.append(OUT / f"{stem}.pdf")
            written.append((figure + tag, page_h, scale, smallest))
        if len(pages) > 1:
            # One PDF per figure, its pages in order; the per-page SVG and
            # PNG stay under their _p<k> names and a stale single-page
            # SVG/PNG of the same figure goes, so nothing on disk can be
            # mistaken for the whole figure.
            import fitz
            out = fitz.open()
            for f in page_pdfs:
                with fitz.open(f) as d:
                    out.insert_pdf(d)
            out.save(OUT / f"{figure}.pdf", garbage=4, deflate=True)
            out.close()
            for f in page_pdfs:
                f.unlink()
            for ext in (".svg", ".png"):
                (OUT / f"{figure}{ext}").unlink(missing_ok=True)
            print(f"    {figure}.pdf: {len(pages)} pages")
    print(f"\nAssembled figures written to {OUT}")
    print(f"\n{'figure':<44}{'page mm':>14}{'journal fit':>13}"
          f"{'smallest pt':>13}")
    worst = 99.0
    for figure, page_h, scale, smallest in written:
        worst = min(worst, smallest)
        print(f"{figure:<44}{PAGE_W:6.1f} x {page_h:5.1f}{scale:13.3f}"
              f"{smallest:13.2f}")
    print(f"\nAt 1:1 - which is how these are distributed - the smallest body "
          f"type on every page above is {style.tick_pt():g} pt.")
    print(f"Were a page instead scaled into {PRINTABLE_W:.0f} x "
          f"{PRINTABLE_H:.0f} mm, the smallest would become {worst:.2f} pt.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
