#!/usr/bin/env python3
"""
S1_E: Doublet Score Distribution
Histogram of Scrublet doublet scores with predicted doublet overlay.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import FULL_DATASET_H5AD

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

# =============================================================================
# Configuration
# =============================================================================
SCALE = 4
FONT_SIZE = 6 * SCALE
DPI = 300
FIG_WIDTH_CM = 8 * SCALE
FIG_HEIGHT_CM = 6 * SCALE
FIG_WIDTH = FIG_WIDTH_CM / 2.54
FIG_HEIGHT = FIG_HEIGHT_CM / 2.54

OUTPUT_DIR = Path(__file__).parent
OUTPUT_FILE = OUTPUT_DIR / "panel_S1_E.png"

# =============================================================================
# Load data
# =============================================================================
print("Loading data...")
adata = sc.read_h5ad(FULL_DATASET_H5AD, backed='r')

doublet_score = adata.obs['doublet_score'].values.astype(float)
predicted_doublet = adata.obs['predicted_doublet'].values.astype(bool)

n_singlet = (~predicted_doublet).sum()
n_doublet = predicted_doublet.sum()
print(f"Singlets: {n_singlet}, Doublets: {n_doublet} ({100*n_doublet/len(predicted_doublet):.1f}%)")

# =============================================================================
# Create figure
# =============================================================================
plt.rcParams.update({
    'font.size': FONT_SIZE,
    'axes.labelsize': FONT_SIZE,
    'axes.titlesize': FONT_SIZE,
    'xtick.labelsize': FONT_SIZE * 0.8,
    'ytick.labelsize': FONT_SIZE * 0.8,
    'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

# Histogram of singlets vs doublets
bins = np.linspace(0, 0.3, 61)
ax.hist(doublet_score[~predicted_doublet], bins=bins, alpha=0.7, color='#4daf4a',
        label=f'Singlets (n={n_singlet:,})', density=True)
ax.hist(doublet_score[predicted_doublet], bins=bins, alpha=0.7, color='#e41a1c',
        label=f'Doublets (n={n_doublet:,})', density=True)

ax.set_xlabel('Doublet Score')
ax.set_ylabel('Density')
ax.set_title('Doublet Score Distribution', fontsize=FONT_SIZE * 0.9)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for spine in ['bottom', 'left']:
    ax.spines[spine].set_linewidth(0.5 * SCALE)
ax.tick_params(width=0.5 * SCALE, length=3 * SCALE)

ax.set_xlim(0.0, 0.3)
ax.legend(fontsize=FONT_SIZE * 0.6, frameon=False, loc='upper right')

plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_FILE.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

print(f"Saved: {OUTPUT_FILE}")
print("Done!")
