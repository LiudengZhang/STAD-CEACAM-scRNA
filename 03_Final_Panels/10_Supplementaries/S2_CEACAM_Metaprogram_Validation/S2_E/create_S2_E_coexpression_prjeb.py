#!/usr/bin/env python3
"""
S2_E: CEACAM5/CEACAM6 Co-expression Scatter Plot — PRJEB25780 (TIGER)
Uses BayesPrism deconvolved epithelial expression. Colored by R vs NR.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *

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
OUTPUT_FILE = OUTPUT_DIR / "panel_S2_E.png"

# Data paths (from paths.py: TIGER_BAYESPRISM_EPI, TIGER_META)
TIGER_EPI_EXPR = TIGER_BAYESPRISM_EPI

RESPONSE_COLORS = {'R': '#2166AC', 'NR': '#B2182B'}

# Load data
print("Loading TIGER BayesPrism epithelial expression...")
epi = pd.read_csv(TIGER_EPI_EXPR, sep='\t', index_col=0)
if 'CEACAM5' in epi.index:
    epi = epi.T
print(f"  Shape: {epi.shape}")

print("Loading metadata...")
meta = pd.read_csv(TIGER_META, sep='\t')
resp_map = meta.set_index('sample_id')['response_NR'].to_dict()
resp_map = {k: ('R' if v == 'R' else 'NR') for k, v in resp_map.items() if pd.notna(v)}

# Build dataframe
common = [s for s in epi.index if s in resp_map]
df = pd.DataFrame({
    'CEACAM5': np.log2(epi.loc[common, 'CEACAM5'].values + 1),
    'CEACAM6': np.log2(epi.loc[common, 'CEACAM6'].values + 1),
    'Response': [resp_map[s] for s in common],
}, index=common)
df = df[df['Response'].isin(['R', 'NR'])]
print(f"  Samples: {len(df)} (R={sum(df.Response=='R')}, NR={sum(df.Response=='NR')})")

# Correlation
r, p = stats.spearmanr(df['CEACAM5'], df['CEACAM6'])
print(f"  Spearman ρ = {r:.3f}, p = {p:.2e}")

# Plot
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

for response in ['NR', 'R']:
    subset = df[df['Response'] == response]
    ax.scatter(subset['CEACAM5'], subset['CEACAM6'],
               c=RESPONSE_COLORS[response], s=30 * SCALE,
               alpha=0.7, edgecolors='white', linewidths=0.3 * SCALE,
               label=f'{response} (n={len(subset)})', zorder=3)

slope, intercept, _, _, _ = stats.linregress(df['CEACAM5'], df['CEACAM6'])
x_line = np.linspace(df['CEACAM5'].min(), df['CEACAM5'].max(), 100)
ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=LINEWIDTH * 1.5, alpha=0.7)

sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
ax.text(0.05, 0.95, f'ρ = {r:.2f} ({sig})\nn = {len(df)}',
        transform=ax.transAxes, fontsize=FONT_SIZE * 0.8,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax.set_xlabel('Epi. CEACAM5 (log$_2$+1)', fontsize=FONT_SIZE)
ax.set_ylabel('Epi. CEACAM6 (log$_2$+1)', fontsize=FONT_SIZE)
ax.set_title('CEACAM5/6 Co-expression\n(PRJEB25780)', fontsize=FONT_SIZE * 0.9)

ax.legend(loc='lower right', fontsize=FONT_SIZE * 0.7, frameon=True)
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
