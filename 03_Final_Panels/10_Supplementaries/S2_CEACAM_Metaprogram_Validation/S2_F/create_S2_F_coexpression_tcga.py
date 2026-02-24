#!/usr/bin/env python3
"""
S2_F: CEACAM5/CEACAM6 Co-expression Scatter Plot — TCGA-STAD
Uses BayesPrism deconvolved epithelial expression from tumor-only samples.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import TCGA_BAYESPRISM_EPI

SCALE = 4
FONT_SIZE = 6 * SCALE
LINEWIDTH = 0.5 * SCALE
DPI = 300

PRINT_WIDTH_MM = 60
PRINT_HEIGHT_MM = 60
MM_TO_INCH = 1 / 25.4
FIG_WIDTH = PRINT_WIDTH_MM * MM_TO_INCH * SCALE
FIG_HEIGHT = PRINT_HEIGHT_MM * MM_TO_INCH * SCALE

OUTPUT_DIR = Path(__file__).parent
OUTPUT_FILE = OUTPUT_DIR / "panel_S2_F.png"

print("Loading TCGA BayesPrism epithelial expression...")
epi = pd.read_csv(TCGA_BAYESPRISM_EPI, sep='\t', index_col=0)
if 'CEACAM5' in epi.index:
    epi = epi.T
print(f"Shape: {epi.shape}")

ceacam5 = np.log2(epi['CEACAM5'].values + 1)
ceacam6 = np.log2(epi['CEACAM6'].values + 1)

r, p = stats.spearmanr(ceacam5, ceacam6)
print(f"Spearman correlation: ρ = {r:.3f}, p = {p:.2e}")

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'font.size': FONT_SIZE,
    'axes.linewidth': LINEWIDTH,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

ax.scatter(ceacam5, ceacam6, c='#666666', s=20 * SCALE,
           alpha=0.5, edgecolors='white', linewidths=0.3 * SCALE, zorder=3)

slope, intercept, _, _, _ = stats.linregress(ceacam5, ceacam6)
x_line = np.linspace(ceacam5.min(), ceacam5.max(), 100)
ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=LINEWIDTH * 1.5, alpha=0.7)

sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
ax.text(0.05, 0.95, f'ρ = {r:.2f} ({sig})\nn = {len(ceacam5)}',
        transform=ax.transAxes, fontsize=FONT_SIZE * 0.8,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax.set_xlabel('CEACAM5 (log$_2$+1)', fontsize=FONT_SIZE)
ax.set_ylabel('CEACAM6 (log$_2$+1)', fontsize=FONT_SIZE)
ax.set_title('CEACAM5/6 Co-expression\n(TCGA-STAD)', fontsize=FONT_SIZE * 0.9)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for spine in ['bottom', 'left']:
    ax.spines[spine].set_linewidth(LINEWIDTH)
ax.tick_params(width=0.5 * SCALE, length=3 * SCALE, labelsize=5 * SCALE)

plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_FILE.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {OUTPUT_FILE}")
