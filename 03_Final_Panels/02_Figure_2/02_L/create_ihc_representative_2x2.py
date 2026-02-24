#!/usr/bin/env python3
"""
Panel 2L: 1x6 representative IHC images
NR (P01): H&E, CEACAM5, CEACAM6  |  R (P22): H&E, CEACAM5, CEACAM6
Auto-crops to largest tissue fragment via connected component analysis.
"""

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import IHC_THUMBNAILS

# 4x scaling
SCALE = 4
DPI = 300
CM_TO_INCH = 1 / 2.54

# Panel: 1 row x 6 cols — wide and short
PANEL_W_CM = 12.0 * SCALE
PANEL_H_CM = 2.5 * SCALE

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
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    fig, axes = plt.subplots(
        1, 6,
        figsize=(PANEL_W_CM * CM_TO_INCH, PANEL_H_CM * CM_TO_INCH),
    )

    # Layout: 1 row x 6 cols — NR(H&E, CEACAM5, CEACAM6), R(H&E, CEACAM5, CEACAM6)
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
                # Rotate 45° CCW for tissue orientation
                from scipy.ndimage import rotate as ndi_rotate
                img_cropped = ndi_rotate(img_cropped, 45, reshape=True,
                                         order=1, cval=255)
                h, w = img_cropped.shape[:2]
                margin = int(min(h, w) * 0.18)
                img_cropped = img_cropped[margin:h-margin, margin:w-margin]
                img_cropped = make_square_crop(img_cropped)
            if marker == "CEACAM6":
                # Rot 180° to match CEACAM5 orientation
                img_cropped = np.rot90(img_cropped, k=2)
                if response == "R":
                    # Additional 45° CCW rotation + zoom in
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
                    ha='center', va='center', fontsize=6 * SCALE)

        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
            spine.set_color('black')

    # Column titles — marker names above each image
    marker_titles = ["H&E", "CEACAM5", "CEACAM6", "H&E", "CEACAM5", "CEACAM6"]
    marker_styles = ["normal", "italic", "italic", "normal", "italic", "italic"]
    for i, (title, style) in enumerate(zip(marker_titles, marker_styles)):
        axes[i].set_title(title, fontsize=7 * SCALE, fontweight='regular',
                          fontstyle=style, pad=8 * SCALE)

    # Group labels — NR / R below each triplet
    # Use xlabel on the middle column of each group
    axes[1].set_xlabel("NR", fontsize=7 * SCALE, labelpad=8 * SCALE)
    axes[4].set_xlabel("R", fontsize=7 * SCALE, labelpad=8 * SCALE)

    plt.subplots_adjust(wspace=0.08)

    # Save PNG and SVG
    for ext in ['png', 'svg']:
        out_path = OUTPUT_DIR / f"ihc_representative_1x6.{ext}"
        fig.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\nSaved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
