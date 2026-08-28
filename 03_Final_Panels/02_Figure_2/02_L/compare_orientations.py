#!/usr/bin/env python3
"""
Generate comparison grid: CEACAM5 as reference, 8 orientations of CEACAM6.
Row 1: CEACAM5 NR (reference) | CEACAM5 R (reference)
Rows 2-5: 8 transformations of CEACAM6 (4 rotations × with/without horizontal mirror)
"""

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import IHC_THUMBNAILS  # noqa: E402

THUMB_DIR = IHC_THUMBNAILS
OUTPUT_DIR = Path(__file__).parent

FILES = {
    "CEACAM5_NR": THUMB_DIR / "CEACAM5" / "P01_CEACAM5.png",
    "CEACAM5_R":  THUMB_DIR / "CEACAM5" / "P22_CEACAM5.png",
    "CEACAM6_NR": THUMB_DIR / "CEACAM6" / "P01_CEACAM6.png",
    "CEACAM6_R":  THUMB_DIR / "CEACAM6" / "P22_CEACAM6.png",
}

TRANSFORMS = [
    ("1: Original",           lambda img: img),
    ("2: Rot 90° CCW",        lambda img: np.rot90(img, k=1)),
    ("3: Rot 180°",           lambda img: np.rot90(img, k=2)),
    ("4: Rot 90° CW",         lambda img: np.rot90(img, k=3)),
    ("5: Mirror LR",          lambda img: np.fliplr(img)),
    ("6: Mirror UD",          lambda img: np.flipud(img)),
    ("7: Rot 90° CCW + Mirror LR", lambda img: np.fliplr(np.rot90(img, k=1))),
    ("8: Rot 90° CW + Mirror LR",  lambda img: np.fliplr(np.rot90(img, k=3))),
]


def crop_to_largest_fragment(img_array, pad_frac=0.12):
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
    rmin, rmax = max(0, rmin - pad_r), min(h, rmax + pad_r)
    cmin, cmax = max(0, cmin - pad_c), min(w, cmax + pad_c)
    return img_array[rmin:rmax, cmin:cmax]


def make_square_crop(img_array):
    h, w = img_array.shape[:2]
    size = min(h, w)
    top = (h - size) // 2
    left = (w - size) // 2
    return img_array[top:top + size, left:left + size]


def load_and_crop(path):
    img = np.array(Image.open(path).convert("RGB"))
    img = crop_to_largest_fragment(img)
    img = make_square_crop(img)
    return img


def main():
    # Load all images
    c5_nr = load_and_crop(FILES["CEACAM5_NR"])
    c5_r  = load_and_crop(FILES["CEACAM5_R"])
    c6_nr = load_and_crop(FILES["CEACAM6_NR"])
    c6_r  = load_and_crop(FILES["CEACAM6_R"])

    # Grid: 9 rows × 2 cols (row 0 = CEACAM5 reference, rows 1-8 = transforms)
    n_rows = 1 + len(TRANSFORMS)
    fig, axes = plt.subplots(n_rows, 2, figsize=(8, 4 * n_rows))

    # Row 0: CEACAM5 reference
    axes[0, 0].imshow(c5_nr)
    axes[0, 0].set_title("CEACAM5 NR (REFERENCE)", fontsize=10, fontweight='bold', color='green')
    axes[0, 1].imshow(c5_r)
    axes[0, 1].set_title("CEACAM5 R (REFERENCE)", fontsize=10, fontweight='bold', color='green')

    # Rows 1-8: CEACAM6 transforms
    for i, (label, transform) in enumerate(TRANSFORMS):
        row = i + 1
        axes[row, 0].imshow(transform(c6_nr))
        axes[row, 0].set_title(f"CEACAM6 NR — {label}", fontsize=9)
        axes[row, 1].imshow(transform(c6_r))
        axes[row, 1].set_title(f"CEACAM6 R — {label}", fontsize=9)

    for ax in axes.flat:
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    out = OUTPUT_DIR / "orientation_comparison.png"
    fig.savefig(out, dpi=150, facecolor='white')
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
