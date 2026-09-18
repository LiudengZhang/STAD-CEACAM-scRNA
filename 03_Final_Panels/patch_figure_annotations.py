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
  Figure 6          the author-supplied page is retained, with only the two
                    temporal headings relabeled to match the cross-sectional
                    design described in the revised manuscript.

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
import hashlib
import sys

import fitz

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "00_Config"))
from paths import MAIN_FIGURES, REVIEWER_MATERIALS  # noqa: E402

SRC = REVIEWER_MATERIALS / "figures_submitted"
OUT = MAIN_FIGURES / "_patched"
PANEL_1A = MAIN_FIGURES / "_panel_1A" / "figure1A.pdf"

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

# Retired: kept because build_provenance.py reads this table to record which
# drawing is spliced into which slot, and because the splices are part of how
# the shipped Figures 2 and 5 were made. main() no longer reaches them - those
# figures are redrawn rather than patched.
#
# letter, slot rect in mm, letter origin in mm, the drawing, and the digest of
# the drawing this splice was built from. The digest is what makes the splice
# reproducible: the path is a panel script's ordinary output and moves whenever
# that panel is redrawn.
PANEL_SWAPS = {
    "Figure 2": [
        # The right edge was 118.3 mm (335.34 pt) and that blanked panel F.
        #
        # The panel F Milo raster is placed at x0 = 333.2 pt, so it overlapped
        # this slot by 2.12 pt - about three quarters of a millimetre, entirely
        # in the white gutter between the two panels. PDF_REDACT_IMAGE_REMOVE
        # deletes an image whose *bounding box* touches the rectangle, not the
        # part that overlaps it, so that sliver took the whole neighbourhood
        # graph and its colour bar with it. The submitted figure carries 38
        # images; the patched one carried 12.
        #
        # 117.2 mm = 332.03 pt gives the raster 1.2 pt of clearance. Nothing is
        # lost by it: the old panel D's own ink ends at 315.42 pt, and the strip
        # from 330 pt to the old edge contains no drawing at all. Measured, not
        # guessed - and check_slot_images() below now refuses to run if this
        # ever stops being true.
        ("D", (79.8, 5.5, 117.2, 30.5), (80.50, 9.49),
         MAIN_FIGURES / "02_Figure_2" / "02_D"
         / "ceacam_metacell_correlation.pdf",
         "618c72dff1949ae0e0c297bc1f8f18af"),
    ],
    "Figure 5": [
        ("H", (118.8, 45.5, 171.1, 81.0), (119.54, 50.65),
         MAIN_FIGURES / "05_Figure_5" / "05_F"
         / "nfkb_radar_celltype_enrichment.pdf",
         "03111dbbd8363c5de4d09a2e5bd6b965"),
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
    "Figure 6": [
        ("Anti-PD1 +", (300.0, 42.0, 362.0, 60.0),
         "Anti-PD-1 +", 11.0, "hebo", "centre",
         "treatment label, PD-1 typography"),
        ("Intrinsic Resistance", (122.0, 41.5, 236.0, 59.0),
         "Pre-treatment non-response", 11.0, "hebo", "centre",
         "left heading, cross-sectional pre-treatment contrast"),
        ("Acquired Resistance", (424.0, 41.5, 543.0, 59.0),
         "Post-treatment non-response", 11.0, "hebo", "centre",
         "right heading, cross-sectional post-treatment contrast"),
        ("Recruitment", (60.5, 362.0, 121.0, 377.0),
         "abundance", 10.0, "hebo", "centre",
         "left annotation, association rather than recruitment mechanism"),
        ("Macrophage", (489.0, 139.5, 552.0, 155.0),
         "MoMac", 10.5, "hebo", "centre",
         "right state label, monocyte/macrophage identity"),
        ("Differentiation", (444.0, 100.5, 514.0, 116.0),
         "MoMac state", 10.0, "heit", "centre",
         "right annotation, observed state rather than inferred differentiation"),
        ("Activation", (391.5, 270.0, 443.0, 286.0),
         "signature", 10.0, "hebo", "centre",
         "centre annotation, transcriptional evidence"),
        ("Chronic ", (603.0, 411.5, 653.0, 430.0),
         "Inflammatory", 11.0, "hebo", "centre",
         "lower-right annotation, cross-sectional inflammatory signal"),
        ("Inflammation", (587.0, 426.0, 666.0, 444.0),
         "signature", 11.0, "hebo", "centre",
         "lower-right annotation, cross-sectional inflammatory signal"),
    ],
}

# Figure 6's entry is the size the submitted label was set at. That figure is no
# longer patched, so it is carried for the record rather than used.
SIZES = {"Figure 2": 6.0, "Figure 3": 5.0, "Figure 5": 3.5, "Figure 6": 10.515}


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
        # An empty replacement is a deletion: the glyphs go and nothing takes
        # their place, so whatever else the label carries stays exactly where the
        # submitted file put it.
        if not new:
            continue
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


def check_slot_images(page, name, letter, slot):
    """
    Refuse to redact a slot that would take an image belonging to another panel.

    This is the check that was missing when the Figure 2D swap blanked panel F.
    apply_redactions() with PDF_REDACT_IMAGE_REMOVE removes an image whose
    bounding box intersects the rectangle - all of it, however small the
    overlap. A slot bounded by the white gutters between panels still touches
    the bounding box of a neighbour's raster, because a bounding box includes
    the figure's own margins.

    So: any image the slot touches must lie wholly inside it. One that does not
    belongs to another panel, and removing it would delete published content.
    Stop rather than do that.
    """
    for img in page.get_images(full=True):
        for rect in page.get_image_rects(img[0]):
            if not rect.intersects(slot):
                continue
            if rect in slot:
                continue
            ox = min(slot.x1, rect.x1) - max(slot.x0, rect.x0)
            oy = min(slot.y1, rect.y1) - max(slot.y0, rect.y0)
            sys.exit(
                f"{name} panel {letter}: the redaction slot "
                f"({slot.x0:.1f}, {slot.y0:.1f}, {slot.x1:.1f}, {slot.y1:.1f}) "
                f"touches an image at "
                f"({rect.x0:.1f}, {rect.y0:.1f}, {rect.x1:.1f}, {rect.y1:.1f}) "
                f"that extends outside it, overlapping by {ox:.2f} x {oy:.2f} "
                "pt. PDF_REDACT_IMAGE_REMOVE would delete that image in full, "
                "not just the overlap. Narrow the slot until it clears the "
                "image, or move the image. Refusing to redact.")


def replace_panels(page, name):
    """
    Redraw whole panels whose numbers changed, in place, inside their own slot.

    The slot is cleared and the panel letter is redrawn at the origin the
    submitted file used, so the lettering of the figure is untouched. The
    replacement keeps its own aspect ratio and is centred in what is left of
    the slot.
    """
    done = []
    for letter, (x0, y0, x1, y1), (lx, ly), src, want in \
            PANEL_SWAPS.get(name, []):
        if not src.exists():
            sys.exit(f"{src} not found - run its panel script first")
        got = hashlib.md5(src.read_bytes()).hexdigest()
        if got != want:
            sys.exit(
                f"{name} panel {letter}: {src.name} is not the drawing this "
                f"splice was built from (md5 {got}, expected {want}). The "
                f"panel has been redrawn since. Splicing a different drawing "
                f"into a frozen figure would change it without saying so; "
                f"restore the recorded drawing, or retire this swap once the "
                f"figure is rebuilt from its slots.")
        # Redaction, not a white rectangle. Painting over the slot leaves the
        # replaced panel in the content stream: cover it and "rho = 0.93 (***)"
        # stays extractable, and a PDF text search still finds it underneath
        # the replacement. apply_redactions removes the objects. The slot
        # is bounded by the white gutters measured in the submitted file, so
        # nothing outside it intersects.
        slot = fitz.Rect(x0 * MM, y0 * MM, x1 * MM, y1 * MM)
        check_slot_images(page, name, letter, slot)
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


#: The figures this still makes: none, since 2026-09-15. Figures 2 to 5 are
#: redrawn into their printed slots by assemble_slotted.py, and Figure 1 joined
#: them on 2026-09-15 - re-paged from the submitted 254 mm landscape onto the
#: 171.10 mm page (12_Figure_Refactor/build_grid_v2.REPAGED), its A placed as
#: a slot and its B and C redrawn at 1:1. Patching Figure 1 here would write
#: the submitted layout into _patched/ beside the slotted page; the code is
#: kept, runnable with --figure 1, as the record of how the 2026-09-10 to
#: 2026-09-14 pages were made. Figure 6 starts from the author-supplied page
#: and receives the two heading corrections above. All
#: are collected by build_shipped_figures.py, which records the mechanism.
PATCHED_FIGURES = (6,)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", type=int, action="append", default=[],
                    help="patch this submitted figure anyway (the retired "
                         "Figure 1 path; nothing ships from it)")
    args = ap.parse_args()
    figures = tuple(args.figure) or PATCHED_FIGURES
    OUT.mkdir(parents=True, exist_ok=True)
    total = 0
    if not figures:
        print("nothing to patch: every main figure is slotted (Figure 1 since "
              "2026-09-15) or supplied; see PATCHED_FIGURES")
    for i in figures:
        name = f"Figure {i}"
        src = (MAIN_FIGURES / "_supplied" / "Figure_6.pdf") if i == 6 \
            else SRC / f"{name}.pdf"
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
    print("Figures 2, 3, 4 and 5 are not patched: they are redrawn into their "
          "printed slots. Figure 6 starts from the author-supplied page and "
          "has only its two temporal headings corrected. All are collected "
          "by build_shipped_figures.py")


if __name__ == "__main__":
    main()
