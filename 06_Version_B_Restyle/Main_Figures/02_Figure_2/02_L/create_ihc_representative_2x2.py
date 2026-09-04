#!/usr/bin/env python3
"""
Figure 2, printed panel M, RESTYLED (Version B) - 1x6 representative IHC
images.

NR: H&E, CEACAM5, CEACAM6  |  R: H&E, CEACAM5, CEACAM6
Auto-crops to largest tissue fragment via connected component analysis.

Version A is
`03_Revised_Panels/Main_Figures/02_Figure_2/02_L/create_ihc_representative_2x2.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every image file, every crop, every rotation, every margin fraction inside the
image processing and every string is Version A's. The drawing code is the same
code.

  printed panel  Figure 2 M     (PROVENANCE.csv; NOT inferred from "02_L")
  printed rect   123.2 x 30.3 mm    (panel_rects.csv)
  Version B box  160.0 x 37.0 mm

  NOTE on the file name, carried over unchanged and flagged, not fixed:
  Version A's script is called `create_ihc_representative_2x2.py`, draws a 1x6
  layout, and writes `ihc_representative_1x6.*`, while PROVENANCE.csv records
  the panel's source file as `02_L/ihc_representative_2x3.svg` and
  `assemble_figure_2.py` reads that name. Three older layouts (2x2, 2x3, 3x2)
  still sit beside it in the Version A directory. The output stem here is
  Version A's own - `ihc_representative_1x6` - as PANEL_SPEC requires; which
  file the Version B assembler should place is a question for the assembler,
  not for this panel.

MARK
    Version A drew at SCALE = 4 (12.0 x 2.5 cm x 4 = 480 x 100 mm) and the body
    type it actually draws is the six column titles and the two group labels,
    all at `7 * SCALE`. So

        SCALE = 4, SMALL_PT = 7
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 28 = 0.25

    (The only smaller size in the file, `6 * SCALE`, is on the "Image not
    found" placeholder, which draws nothing when the images are present. It is
    left at MARK-scaled type for the same reason as everything else.)

    The title `pad` and the group-label `labelpad`, both `8 * SCALE`, are
    lengths in points, not type, so they take MARK: 32 pt becomes 8 pt, which
    is what they printed at. Version A's `spine.set_linewidth(0.5)` on the six
    image borders is NOT carried over: a spine is axes furniture and cnsplots
    sets axes.linewidth.

JUDGEMENT CALL - the image borders
    cnsplots turns the top and right spines off, which is the right default for
    a data axes and the wrong one for a photograph: it would leave each of the
    six micrographs framed on two sides only. Version A's own loop runs over all
    four spines and colours them, so its border is four-sided. `set_visible(True)`
    is therefore added inside that loop. It restores Version A's drawing rather
    than introducing anything, and it is the one place in this panel where the
    library's default is overridden.

    The panel grew from 123.2 x 30.3 mm to 160 x 37 mm so that the six column
    titles ("CEACAM5", "CEACAM6") fit at 8 pt over 24.5 mm-wide images.
"""

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import IHC_THUMBNAILS                          # noqa: E402
import panel_style_cns as style                           # noqa: E402

SCALE = 4                           # Version A's canvas multiplier
SMALL_PT = 7.0                      # Version A's smallest drawn body type

PRINTED_MM = (123.2, 30.3)          # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 160.0, 37.0
MARGIN = dict(left=1.5, right=1.5, top=6.0, bottom=6.5, wspace=0.08)

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


def main():
    family = style.apply()
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
            # set_visible(True) is the one addition to Version A's loop, and it
            # is here to keep Version A's drawing rather than to change it: the
            # loop runs over all four spines and colours them, so the border it
            # draws round each photograph is four-sided. cnsplots turns the top
            # and right spines off, which is right for a data axes and wrong for
            # an image frame - it would leave every micrograph in an L-shaped
            # rule. See JUDGEMENT CALL in the header.
            spine.set_visible(True)
            spine.set_color('black')

    # Column titles - marker names above each image
    marker_titles = ["H&E", "CEACAM5", "CEACAM6", "H&E", "CEACAM5", "CEACAM6"]
    marker_styles = ["normal", "italic", "italic", "normal", "italic", "italic"]
    for i, (title, style_) in enumerate(zip(marker_titles, marker_styles)):
        axes[i].set_title(title, fontstyle=style_, pad=8 * SCALE * MARK)

    # Group labels - NR / R below each triplet
    # Use xlabel on the middle column of each group
    axes[1].set_xlabel("NR", labelpad=8 * SCALE * MARK)
    axes[4].set_xlabel("R", labelpad=8 * SCALE * MARK)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    style.save_panel(fig, OUTPUT_DIR / "ihc_representative_1x6")
    print(f"\nSaved: {OUTPUT_DIR / 'ihc_representative_1x6'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
