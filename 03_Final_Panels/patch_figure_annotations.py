"""
Put exact two-sided P values into the printed main figures.

Reviewer 1 point R1.3c asked for two-sided tests and exact P values. The tests
were converted and the panel scripts now print the value, but the figures in the
paper cannot be rebuilt from those scripts: the panel code was re-lettered after
submission - it now emits Figure 2 as A-Q with split N1/N2 and O1/O2 panels
where the paper prints A-N, and has no panel L at all - and the pre-submission
assembler is neither on disk nor in version control. Re-assembling would
silently renumber the figures and invalidate every panel reference in the text.

So the annotations are edited in place instead, in the submitted PDFs, and
nothing else about the figures changes. Two mechanisms, chosen by what each
figure actually contains:

  Figures 3 and 5   the annotations are live text, so the old glyphs are
                    redacted and the value is drawn at the same anchor
  Figure 2          the annotations are vector outlines with no text layer, so
                    the glyph cluster is covered and the value drawn over it
  Figure 4          two label corrections rather than P values, same mechanism:
                    the gene label "CD16" is not an HGNC symbol (the gene is
                    FCGR3A), and the C3 state is renamed Mac -> MoMac because
                    the lineage analysis added in this revision places it on the
                    monocyte-macrophage continuum. The replacements are longer
                    than what they replace, so each is redrawn on the anchor its
                    neighbours share - right edge for tick labels, midpoint for
                    the UMAP annotation and the rotated axis titles.

Figure 1 is different in kind: panel A is a schematic, not a measurement, and
the author replaced it with a redrawn vector version (_panel_1A/). That panel
is taller relative to its width than the one it replaces, so the page grows
and panels B and C move down unchanged - nothing inside them is rescaled or
re-lettered. The pictograms are openly licensed and the CC BY 4.0 ones need a
credit line, which edits.py adds to the Figure 1 legend.

Every value below traces to twosided_sweep.csv, except the three panels the
sweep did not cover, which are marked and traced to their panel scripts.

Run: python patch_figure_annotations.py
"""

from pathlib import Path
import sys

import fitz

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "00_Config"))
from paths import REVIEWER_MATERIALS  # noqa: E402

SRC = REVIEWER_MATERIALS / "figures_submitted"
OUT = Path(__file__).parent / "Main_Figures" / "_patched"
PANEL_1A = Path(__file__).parent / "Main_Figures" / "_panel_1A" / "figure1A.pdf"

# Two measured panels are redrawn rather than annotated, because the numbers
# behind them were wrong rather than merely unlabelled.
#
#   Figure 2D  the submitted panel reports rho = 0.93 over "n = 49,696
#              pre-treatment epithelial cells". Both come from the damaged .X
#              of Epithelial.h5ad (00_Data_Audit/FINDINGS.md, sections 1 and 7):
#              49,696 is exactly the number of rows of that matrix that are not
#              NaN, out of 106,653, and 0.93 is the correlation among them. Read
#              from .raw the same cells give rho = 0.44. The replacement shows
#              kNN metacells of 10 cells, and prints the per-cell value beside
#              the metacell value.
#   Figure 5H  the submitted radar was built from the same corrupted
#              differential-expression lineage. It disagrees with the clean
#              recompute in the two places the text now makes claims about:
#              B cells after treatment (+1.01 there, -0.99 in the recompute) and
#              monocytes/macrophages before treatment (-0.98 there, +1.06).
#
# Each entry is (figure, panel letter, slot in mm, letter origin in mm, source
# PDF). The slot is bounded by the white gutters measured in the submitted file,
# so no neighbouring panel is touched; LETTER_INSET keeps the panel letter clear
# of the artwork placed beside it.
LETTER_INSET = 3.0      # mm

PANEL_SWAPS = {
    "Figure 2": [
        ("D", (79.8, 5.5, 118.3, 30.5), (80.50, 9.49),
         Path(__file__).parent / "Main_Figures" / "02_Figure_2" / "02_D"
         / "ceacam_metacell_correlation.pdf"),
    ],
    "Figure 5": [
        ("H", (118.8, 45.5, 171.1, 81.0), (119.54, 50.65),
         Path(__file__).parent / "Main_Figures" / "05_Figure_5" / "05_F"
         / "nfkb_radar_celltype_enrichment.pdf"),
    ],
}

MM = 72 / 25.4          # millimetres to PDF points

# Geometry of the submitted Figure 1, measured from the file itself: panel A
# occupies x 13.3-237.7 mm and ends at y 56.9 mm; the first drawing of the
# B/C row starts at y 69.6 mm, and their panel letters sit at y 67.6 mm.
FIG1_A_BAND = 62.0      # mm; everything above this is panel A
FIG1_A_LEFT = 13.3      # mm; left edge of panel A, shared with B
FIG1_A_WIDTH = 224.4    # mm; width of panel A, shared with the B/C row
FIG1_A_TOP = 13.0       # mm; top of the panel-A block, including its letter
FIG1_A_GAP = 10.7       # mm; gap between panel A and the B/C letters
# The submitted page runs to 177.8 mm but its ink stops at 123.4 mm, so the
# copied block is trimmed rather than carrying 54 mm of trailing white into a
# page that is already taller.
FIG1_KEEP_BOTTOM = 126.0  # mm

FONT = "helv"          # metrically close to the ArialMT the figures use
WHITE = (1, 1, 1)

# ---------------------------------------------------------------- Figures 3, 5
# (old text, search rectangle, new text). The rectangle disambiguates repeated
# glyphs; it is the span's own bounding box, widened by a point.
TEXT_PATCHES = {
    "Figure 3": [
        ("**", (244.2, 168.3, 252.6, 177.9), "P = 0.014", "H, epithelial density"),
        ("*", (331.2, 163.8, 337.4, 173.5), "P = 0.084", "I, distance to stroma"),
        ("*", (417.7, 163.2, 423.9, 172.9), "P = 0.064", "J, distance to immune"),
    ],
    "Figure 5": [
        ("ns", (49.2, 159.8, 56.3, 167.5), "P = 0.924", "D, Kruskal-Wallis among others"),
        ("*", (72.0, 152.8, 77.4, 161.3), "P = 0.035", "D, BACH1 regulon"),
        ("ns", (49.2, 254.8, 56.3, 262.5), "P = 0.898", "E, Kruskal-Wallis among others"),
        ("*", (71.9, 246.7, 77.3, 255.2), "P = 0.038", "E, NFKB1 regulon"),
        ("*", (38.8, 334.0, 44.3, 343.4), "P = 0.052", "J, monocytes/macrophages"),
        ("*", (98.0, 334.0, 103.5, 343.4), "P = 0.082", "J, epithelial cells"),
        ("*", (157.1, 334.0, 162.6, 343.4), "P = 0.055", "J, fibroblasts"),
        ("*", (220.1, 338.2, 225.5, 347.6), "P = 0.032", "J, dendritic cells"),
        ("*", (357.3, 338.2, 362.8, 347.6), "P = 0.030", "L, IL-6/JAK/STAT3 in CD4+ T"),
    ],
}

# --------------------------------------------------------------------- Figure 2
# (glyph-cluster box, baseline y, new text). The boxes were read off the page's
# own drawing list, so they are the exact extent of what is being replaced.
DRAW_PATCHES = {
    "Figure 2": [
        ((277.6, 100.8, 298.9, 104.7), 104.7, "P = 0.057", "E, CEACAM5/6+ cluster proportion"),
        ((325.3, 188.3, 327.7, 192.7), 192.7, "P = 0.083", "H, S-MP4"),
        ((415.0, 188.3, 438.0, 192.7), 192.7, "P = 0.188", "I, S-MP5"),
        ((370.1, 301.6, 371.9, 304.8), 304.8, "P = 0.057", "K left, CEACAM6"),
        ((424.1, 301.6, 440.9, 304.8), 304.8, "P = 0.114", "K right, CEACAM5"),
        ((368.3, 372.2, 370.1, 374.4), 374.4, "P = 0.056", "L left, CEACAM6 PRJEB25780"),
        ((435.6, 372.2, 437.4, 374.4), 374.4, "P = 0.059", "L right, CEACAM5 PRJEB25780"),
        ((420.1, 455.1, 421.9, 459.1), 459.1, "P = 0.057", "N, IHC summed"),
    ],
}

# --------------------------------------------------------------------- Figure 4
# (old text, search rectangle, new text, size, font, alignment, note). The
# rectangle is the span's own bounding box widened by a point, as above.
# Alignment: "right" for the tick labels, whose right edges line up with their
# neighbours; "centre" for the UMAP annotation and the two rotated axis titles.
# "rot-top" is a rotated label whose top edge lines up with its neighbours.
LABEL_PATCHES = {
    "Figure 4": [
        ("CD16", (306.6, 381.2, 314.7, 398.6), "FCGR3A", 6.04, "heit", "rot-top",
         "E, gene label, CD16 is not an HGNC symbol"),
        ("Mac_IL1B", (150.1, 340.3, 179.3, 348.3), "MoMac_IL1B", 6.04, "helv",
         "right", "E, state label"),
        ("Mac_IL1B", (288.3, 433.1, 308.3, 439.1), "MoMac_IL1B", 4.0, "helv",
         "right", "G, state label"),
        ("Mac_IL1B", (53.6, 127.8, 70.5, 133.2), "MoMac_IL1B", 3.32, "helv",
         "right", "A, state label"),
        ("Mac_IL1B", (59.9, 467.2, 89.9, 475.2), "MoMac_IL1B", 6.0, "hebo",
         "centre", "F, UMAP annotation"),
        ("Mac_IL1B Proportion", (395.3, 429.5, 401.8, 473.8),
         "MoMac_IL1B Proportion", 4.5, "helv", "rot-centre", "H, axis title"),
        ("Mac_IL1B Proportion", (393.2, 501.5, 399.7, 545.8),
         "MoMac_IL1B Proportion", 4.5, "helv", "rot-centre", "H, axis title"),
    ],
}

SIZES = {"Figure 2": 6.0, "Figure 3": 5.0, "Figure 5": 3.5}


def centred(page, box, baseline, text, size):
    """Draw `text` centred on the horizontal midpoint of `box`."""
    width = fitz.get_text_length(text, fontname=FONT, fontsize=size)
    mid = (box[0] + box[2]) / 2
    page.insert_text((mid - width / 2, baseline), text,
                     fontname=FONT, fontsize=size, color=(0, 0, 0))


def patch_text(page, name, size):
    """Redact the old glyphs, then draw the value at the same anchor."""
    done = []
    for old, box, new, where in TEXT_PATCHES[name]:
        rect = fitz.Rect(*box)
        found = [r for r in page.search_for(old, quads=False) if rect.contains(r)]
        if len(found) != 1:
            raise SystemExit(
                f"{name}: {old!r} at {box} matched {len(found)} spans, expected 1")
        page.add_redact_annot(found[0])
        done.append((found[0], new, where))
    # Text only: line art and images are what the panels are drawn with.
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                          graphics=fitz.PDF_REDACT_LINE_ART_NONE)
    for rect, new, where in done:
        centred(page, (rect.x0, rect.y0, rect.x1, rect.y1), rect.y1, new, size)
    return [w for _, _, w in done]


def patch_drawings(page, name, size):
    """Cover the vector glyph cluster and draw the value over it."""
    done = []
    for box, baseline, new, where in DRAW_PATCHES[name]:
        # Pad sideways and upward only. The significance bracket sits a point or
        # two below the glyphs, and padding downward erases the middle of it.
        page.draw_rect(fitz.Rect(box[0] - 1.5, box[1] - 1.5, box[2] + 1.5, box[3] + 0.2),
                       color=None, fill=WHITE)
        centred(page, box, baseline, new, size)
        done.append(where)
    return done


def patch_labels(page, name):
    """Redact a label and redraw the longer replacement on the shared anchor."""
    done = []
    for old, box, new, size, font, align, where in LABEL_PATCHES[name]:
        rect = fitz.Rect(*box)
        found = [r for r in page.search_for(old, quads=False) if rect.contains(r)]
        if len(found) != 1:
            raise SystemExit(
                f"{name}: {old!r} at {box} matched {len(found)} spans, expected 1")
        page.add_redact_annot(found[0])
        done.append((found[0], new, size, font, align, where))
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                          graphics=fitz.PDF_REDACT_LINE_ART_NONE)
    for rect, new, size, font, align, where in done:
        length = fitz.get_text_length(new, fontname=font, fontsize=size)
        if align.startswith("rot"):
            # Rotated 90 degrees, reading upward: the baseline runs along the
            # right edge of the glyph box and the run grows downward as it gets
            # longer, so the anchor is the top edge or the midpoint in y.
            x = rect.x1 - 0.11 * size
            y = rect.y0 + length if align == "rot-top" else \
                (rect.y0 + rect.y1) / 2 + length / 2
            page.insert_text((x, y), new, fontname=font, fontsize=size,
                             rotate=90, color=(0, 0, 0))
        else:
            x = rect.x1 - length if align == "right" else \
                (rect.x0 + rect.x1) / 2 - length / 2
            page.insert_text((x, rect.y1), new, fontname=font, fontsize=size,
                             color=(0, 0, 0))
    return [w for *_, w in done]


def replace_panels(page, name):
    """
    Redraw whole panels whose numbers changed, in place, inside their own slot.

    The slot is cleared and the panel letter is redrawn at the origin the
    submitted file used, so the lettering of the figure is untouched. The
    replacement keeps its own aspect ratio and is centred in what is left of
    the slot.
    """
    done = []
    for letter, (x0, y0, x1, y1), (lx, ly), src in PANEL_SWAPS.get(name, []):
        if not src.exists():
            sys.exit(f"{src} not found - run its panel script first")
        # Redaction, not a white rectangle. Painting over the slot leaves the
        # superseded panel in the content stream: the first attempt at this left
        # "rho = 0.93 (***)" extractable, and a PDF text search still found it
        # under the replacement. apply_redactions removes the objects. The slot
        # is bounded by the white gutters measured in the submitted file, so
        # nothing outside it intersects.
        slot = fitz.Rect(x0 * MM, y0 * MM, x1 * MM, y1 * MM)
        page.add_redact_annot(slot)
        page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_REMOVE,
                              graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_TOUCHED)
        page.draw_rect(slot, color=WHITE, fill=WHITE, width=0)
        page.insert_text((lx * MM, ly * MM), letter, fontname="hebo",
                         fontsize=10, color=(0, 0, 0))

        art = fitz.open(src)
        ar = art[0].rect
        ax0, ay0 = lx + LETTER_INSET, y0 + 0.5
        aw, ah = x1 - ax0 - 0.3, y1 - ay0 - 0.3
        scale = min(aw / (ar.width / MM), ah / (ar.height / MM))
        w, h = ar.width / MM * scale, ar.height / MM * scale
        ox, oy = ax0 + (aw - w) / 2, ay0 + (ah - h) / 2
        page.show_pdf_page(
            fitz.Rect(ox * MM, oy * MM, (ox + w) * MM, (oy + h) * MM), art, 0)
        art.close()
        done.append(f"{letter}, panel redrawn from {src.parent.name}/{src.name}")
    return done


def replace_figure_1a(doc):
    """
    Swap panel A of Figure 1 for the redrawn vector panel and let the page grow.

    Panels B and C are copied across as a single clipped block, so they keep
    their own scale, position relative to each other and lettering; only their
    vertical offset changes. Returns the new document.
    """
    if not PANEL_1A.exists():
        sys.exit(f"{PANEL_1A} not found")
    src = doc[0]
    a = fitz.open(PANEL_1A)
    ar = a[0].rect
    scale = (FIG1_A_WIDTH * MM) / ar.width
    a_h = ar.height * scale / MM                      # mm

    keep_top = FIG1_A_BAND * MM                       # source y where B/C starts
    keep_bot = min(FIG1_KEEP_BOTTOM * MM, src.rect.height)
    keep_h = (keep_bot - keep_top) / MM               # mm
    # The gap the submitted figure leaves between panel A and the B/C letters is
    # preserved, and the clip starts 5.6 mm above those letters.
    dest_top = FIG1_A_TOP + a_h + FIG1_A_GAP - (67.6 - FIG1_A_BAND)

    out = fitz.open()
    page = out.new_page(width=src.rect.width,
                        height=(dest_top + keep_h) * MM)
    page.show_pdf_page(
        fitz.Rect(FIG1_A_LEFT * MM, FIG1_A_TOP * MM,
                  (FIG1_A_LEFT + FIG1_A_WIDTH) * MM, (FIG1_A_TOP + a_h) * MM),
        a, 0)
    page.show_pdf_page(
        fitz.Rect(0, dest_top * MM, src.rect.width, (dest_top + keep_h) * MM),
        doc, 0, clip=fitz.Rect(0, keep_top, src.rect.width, keep_bot))
    a.close()
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    total = 0
    for i in range(1, 7):
        name = f"Figure {i}"
        src = SRC / f"{name}.pdf"
        doc = fitz.open(src)
        page = doc[0]
        changed = []
        if name in TEXT_PATCHES:
            changed += patch_text(page, name, SIZES[name])
        if name in DRAW_PATCHES:
            changed += patch_drawings(page, name, SIZES[name])
        if name in LABEL_PATCHES:
            changed += patch_labels(page, name)
        if name in PANEL_SWAPS:
            changed += replace_panels(page, name)
        if name == "Figure 1":
            doc = replace_figure_1a(doc)
            changed.append("A, schematic replaced with the redrawn vector panel")
        dst = OUT / f"Figure_{i}.pdf"
        if dst.exists():
            dst.chmod(0o644)
        doc.save(dst, garbage=3, deflate=True)
        doc.close()
        dst.chmod(0o644)
        total += len(changed)
        print(f"{name}: {len(changed)} annotation(s)")
        for w in changed:
            print(f"    {w}")
    print(f"\n{total} annotations replaced; written to {OUT}")


if __name__ == "__main__":
    main()
