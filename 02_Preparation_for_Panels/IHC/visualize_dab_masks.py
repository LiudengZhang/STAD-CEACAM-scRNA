#!/usr/bin/env python3
"""
QC visualization: show original IHC image + DAB detection overlay for each sample.
8 rows (patients) x 4 columns (CEACAM5 original, CEACAM5 mask, CEACAM6 original, CEACAM6 mask)
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from skimage.color import rgb2hed, rgb2gray
from scipy.ndimage import binary_closing, binary_opening
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from paths import IHC_THUMBNAILS

THUMB_DIR = IHC_THUMBNAILS
OUTPUT_DIR = Path(__file__).parent

THRESHOLDS = [0.02, 0.05, 0.07]

SAMPLES = {
    'P01 (NR)': ('P01_CEACAM5.png',  'P01_CEACAM6.png'),
    'P02 (NR)': ('P02_CEACAM5.png',  'P02_CEACAM6.png'),
    'P25 (NR)': ('P25_CEACAM5.png',  'P25_CEACAM6.png'),
    'P26 (NR)': ('P26_CEACAM5.png',  'P26_CEACAM6.png'),
    'P03 (R)':  ('P03_CEACAM5.png',  'P03_CEACAM6.png'),
    'P04 (R)':  ('P04_CEACAM5.png',  'P04_CEACAM6.png'),
    'P21 (R)':  ('P21_CEACAM5.png',  'P21_CEACAM6.png'),
    'P22 (R)':  ('P22_CEACAM5.png',  'P22_CEACAM6.png'),
}


def get_tissue_and_dab(img_array):
    hed = rgb2hed(img_array)
    dab = hed[:, :, 2]
    gray = (rgb2gray(img_array) * 255).astype(np.uint8)
    tissue = gray < 230
    tissue = binary_closing(tissue, iterations=3)
    tissue = binary_opening(tissue, iterations=2)
    return dab, tissue


def make_mask_viz(img_array, mask):
    out = img_array.copy()
    out[mask, 0] = 255
    out[mask, 1] = 50
    out[mask, 2] = 50
    return out


# For each sample: 1 original col + 3 threshold cols, for C5 and C6
# Layout: 8 rows x 8 cols (C5 orig, C5@0.02, C5@0.05, C5@0.07, C6 orig, C6@0.02, C6@0.05, C6@0.07)
n_cols = 2 + 2 * len(THRESHOLDS)  # 8
fig, axes = plt.subplots(8, n_cols, figsize=(4 * n_cols, 5 * 8))
fig.suptitle('IHC DAB Detection QC — threshold comparison (0.02 / 0.05 / 0.07)',
             fontsize=18, fontweight='bold', y=0.998)

for row, (label, (c5_file, c6_file)) in enumerate(SAMPLES.items()):
    for m_idx, (marker, fname) in enumerate([('CEACAM5', c5_file), ('CEACAM6', c6_file)]):
        img_path = THUMB_DIR / marker / fname
        img = np.array(Image.open(img_path).convert('RGB'))
        dab_ch, tissue = get_tissue_and_dab(img)
        tissue_px = np.sum(tissue)

        # Original column
        col_base = m_idx * (1 + len(THRESHOLDS))
        ax = axes[row, col_base]
        ax.imshow(img)
        ax.set_title(f'{label}\n{marker}', fontsize=9, fontweight='bold')
        ax.axis('off')

        # Threshold columns
        for t_idx, thresh in enumerate(THRESHOLDS):
            dab_pos = (dab_ch > thresh) & tissue
            pct = np.sum(dab_pos) / tissue_px * 100 if tissue_px > 0 else 0
            viz = make_mask_viz(img, dab_pos)

            ax = axes[row, col_base + 1 + t_idx]
            ax.imshow(viz)
            ax.set_title(f't={thresh} → {pct:.1f}%', fontsize=9,
                         color='red' if pct > 5 else 'black')
            ax.axis('off')

plt.tight_layout()
out_path = OUTPUT_DIR / 'dab_detection_qc.png'
fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out_path}")
