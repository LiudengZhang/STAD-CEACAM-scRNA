#!/usr/bin/env python3
"""
Figure 2 panel M - representative immunohistochemistry, one row of six.

Non-responder: H&E, CEACAM5, CEACAM6.  Responder: H&E, CEACAM5, CEACAM6.
Each image is auto-cropped to its largest tissue fragment by connected
component analysis.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. The margins are millimetres of paper; they
are not fitted to the ink, because six axes holding images of equal aspect are
resized by matplotlib at draw time and the subplot parameters do not describe
where their ink ends up.

  printed panel  Figure 2 M       (PROVENANCE.csv; NOT inferred from "02_L")

  The output stem is ihc_representative_1x6, which is what this script has
  always written and what the printed panel is: one row of six images, aspect
  3.7. The name ihc_representative_2x3 belongs to a different, two-row drawing
  that is aspect 1.4 and is not the panel on the page.

MARK
    The earlier drawing used a canvas four times the printed size and set the
    body type it actually draws - the six column titles and the two group
    labels - at 7 * SCALE. MARK carries the non-type point sizes across to the
    1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The title pad and the group-label labelpad, both lengths in points rather
    than type, take it. The spine width the earlier drawing set is not carried
    over: a spine is axes furniture and cnsplots sets axes.linewidth.

THE IMAGE BORDERS
    cnsplots turns the top and right spines off, which is the right default for
    a data axes and the wrong one for a photograph: it would leave each of the
    six micrographs framed on two sides only. The earlier drawing's own loop
    runs over all four spines and colours them, so its border is four-sided.
    `set_visible(True)` inside that loop restores it.

Every image file, every crop, every rotation, every margin fraction inside the
image processing and every string is the earlier drawing's. The drawing code is
the same code.
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from PIL import Image
from scipy import ndimage
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import IHC_THUMBNAILS                          # noqa: E402
import panel_style_cns as style                           # noqa: E402
import slots                                              # noqa: E402

SCALE = 4                           # the earlier canvas multiplier
SMALL_PT = 7.0                      # the earlier smallest drawn body type

PANEL_LETTER = "M"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(2, PANEL_LETTER)
MARGIN = dict(left=0.6, right=0.6, top=5.3, bottom=5.3, wspace=0.08)

# Paths
THUMB_DIR = IHC_THUMBNAILS
OUTPUT_DIR = Path(__file__).parent

# Sample files
SAMPLES = {
    "NR": {
        "label": "Non-responder (P01)",
        "HE": THUMB_DIR / "HE" / "P01_HE.png",
        "CEACAM5": THUMB_DIR / "CEACAM5" / "P01_CEACAM5.png",
        "CEACAM6": THUMB_DIR / "CEACAM6" / "P01_CEACAM6.png",
    },
    "R": {
        "label": "Responder (P22)",
        "HE": THUMB_DIR / "HE" / "P22_HE.png",
        "CEACAM5": THUMB_DIR / "CEACAM5" / "P22_CEACAM5.png",
        "CEACAM6": THUMB_DIR / "CEACAM6" / "P22_CEACAM6.png",
    },
}


def crop_to_largest_fragment(img_array, pad_frac=0.10):
    """Find the largest connected tissue fragment and crop to it."""
    gray = np.mean(img_array[:, :, :3], axis=2)
    if gray.max() <= 1.0:
        gray = (gray * 255).astype(np.uint8)
    else:
        gray = gray.astype(np.uint8)

    tissue_mask = gray < 215
    tissue_mask = ndimage.binary_closing(tissue_mask, iterations=5)
    tissue_mask = ndimage.binary_opening(tissue_mask, iterations=3)

    labeled, n_features = ndimage.label(tissue_mask)
    if n_features == 0:
        return img_array

    component_sizes = ndimage.sum(tissue_mask, labeled, range(1, n_features + 1))
    largest_id = np.argmax(component_sizes) + 1

    largest_mask = labeled == largest_id
    rows = np.where(np.any(largest_mask, axis=1))[0]
    cols = np.where(np.any(largest_mask, axis=0))[0]

    rmin, rmax = rows[0], rows[-1]
    cmin, cmax = cols[0], cols[-1]

    h, w = img_array.shape[:2]
    pad_r = int((rmax - rmin) * pad_frac)
    pad_c = int((cmax - cmin) * pad_frac)

    rmin = max(0, rmin - pad_r)
    rmax = min(h, rmax + pad_r)
    cmin = max(0, cmin - pad_c)
    cmax = min(w, cmax + pad_c)

    print(f"  Fragment: {cmax-cmin}x{rmax-rmin} px "
          f"(from {w}x{h}, {n_features} fragments)")

    return img_array[rmin:rmax, cmin:cmax]


def make_square_crop(img_array):
    """Crop to square from center."""
    h, w = img_array.shape[:2]
    size = min(h, w)
    top = (h - size) // 2
    left = (w - size) // 2
    return img_array[top:top + size, left:left + size]


#: Millimetres between the foot of an image and its grouping rule.
RULE_DROP_MM = 1.4


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, 1, 6)

    # Layout: 1 row x 6 cols - NR(H&E, CEACAM5, CEACAM6), R(H&E, CEACAM5, CEACAM6)
    layout = [
        ("HE", "NR", axes[0]),
        ("CEACAM5", "NR", axes[1]),
        ("CEACAM6", "NR", axes[2]),
        ("HE", "R", axes[3]),
        ("CEACAM5", "R", axes[4]),
        ("CEACAM6", "R", axes[5]),
    ]

    for marker, response, ax in layout:
        img_path = SAMPLES[response][marker]
        print(f"{marker} {response}: {img_path.name}")
        if img_path.exists():
            img = np.array(Image.open(img_path).convert("RGB"))
            # Zoom out more for R images so tissue fits entirely in box
            pad = 0.20 if response == "R" else 0.12
            img_cropped = crop_to_largest_fragment(img, pad_frac=pad)
            img_cropped = make_square_crop(img_cropped)
            if marker == "HE":
                # Rotate 45 deg CCW for tissue orientation
                from scipy.ndimage import rotate as ndi_rotate
                img_cropped = ndi_rotate(img_cropped, 45, reshape=True,
                                         order=1, cval=255)
                h, w = img_cropped.shape[:2]
                margin = int(min(h, w) * 0.18)
                img_cropped = img_cropped[margin:h-margin, margin:w-margin]
                img_cropped = make_square_crop(img_cropped)
            if marker == "CEACAM6":
                # Rot 180 deg to match CEACAM5 orientation
                img_cropped = np.rot90(img_cropped, k=2)
                if response == "R":
                    # Additional 45 deg CCW rotation + zoom in
                    from scipy.ndimage import rotate as ndi_rotate
                    img_cropped = ndi_rotate(img_cropped, 45, reshape=True,
                                             order=1, cval=255)
                    # Crop center to remove white corners
                    h, w = img_cropped.shape[:2]
                    margin_h = int(h * 0.18)
                    margin_w = int(w * 0.18)
                    img_cropped = img_cropped[margin_h:h-margin_h, margin_w:w-margin_w]
                    img_cropped = make_square_crop(img_cropped)
            ax.imshow(img_cropped)
        else:
            ax.text(0.5, 0.5, "Image not found", transform=ax.transAxes,
                    ha='center', va='center', fontsize=style.tick_pt())

        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            # The loop colours all four spines, so the border round each
            # photograph is four-sided; cnsplots turns the top and right
            # spines off, which is right for a data axes and wrong for an
            # image frame. See THE IMAGE BORDERS in the header.
            spine.set_visible(True)
            spine.set_color('black')

    # Column titles - marker names above each image
    marker_titles = ["H&E", "CEACAM5", "CEACAM6", "H&E", "CEACAM5", "CEACAM6"]
    # Upright, all six (2026-09-14, evening): these are immunohistochemical
    # stains, so the names are the PROTEINS' - a gene symbol is italic, its
    # protein is not, and the submitted page set them upright. The page
    # gate holds this panel to the regular face (check_restyled_panel.py:
    # UPRIGHT_PROTEIN).
    marker_styles = ["normal", "normal", "normal", "normal", "normal", "normal"]
    for i, (title, style_) in enumerate(zip(marker_titles, marker_styles)):
        axes[i].set_title(title, fontstyle=style_, pad=8 * SCALE * MARK)

    # Group labels below each triplet, on the middle column of each, as the
    # shipped page prints them. They hang under the grouping rules drawn after
    # the margins are set; see THE TWO GROUPING RULES.
    axes[1].set_xlabel("Pre-NR", labelpad=8 * SCALE * MARK)
    axes[4].set_xlabel("Pre-R", labelpad=8 * SCALE * MARK)

    style.margins_mm(fig, **MARGIN)

    # THE TWO GROUPING RULES, RESTORED  (2026-09-11)
    #   The published page draws a rule under each triplet of images and hangs
    #   the group name off it, which is what tells the reader that the first
    #   three images are one patient group and the last three another. The
    #   redraw kept the two names and dropped both rules, so six images sat in
    #   a row with two words under them and nothing saying where one group
    #   ended. Neither rule carries a value; each spans exactly the images it
    #   groups, read from their own placed positions rather than typed.
    for first, last in ((0, 2), (3, 5)):
        a = axes[first].get_position()
        b = axes[last].get_position()
        y = a.y0 - RULE_DROP_MM / PANEL_H_MM
        fig.add_artist(mlines.Line2D([a.x0, b.x1], [y, y],
                                     transform=fig.transFigure,
                                     color='black', linewidth=style.RULE_PT,
                                     solid_capstyle='butt'))

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, OUTPUT_DIR / "ihc_representative_1x6")
    print(f"\nSaved: {OUTPUT_DIR / 'ihc_representative_1x6'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
