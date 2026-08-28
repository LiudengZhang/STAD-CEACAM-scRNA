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
