#!/usr/bin/env python3
"""
S2_D: TCGA-STAD Kaplan-Meier Survival by CEACAM5/6 Expression
Uses BayesPrism deconvolved epithelial expression (median split).
Proper KM curves with lifelines + log-rank test.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import TCGA_BAYESPRISM_EPI, TCGA_CLINICAL

# ==============================================================================
# Configuration
# ==============================================================================
SCALE = 4
FONT_SIZE = 6 * SCALE
DPI = 300

PRINT_WIDTH_MM = 100
PRINT_HEIGHT_MM = 55
MM_TO_INCH = 1 / 25.4
FIG_WIDTH = PRINT_WIDTH_MM * MM_TO_INCH * SCALE
FIG_HEIGHT = PRINT_HEIGHT_MM * MM_TO_INCH * SCALE

SCRIPT_DIR = Path(__file__).parent
CLINICAL_FILE = TCGA_CLINICAL
OUTPUT_FILE = SCRIPT_DIR / "panel_S2_D.png"

COLORS = {'High': '#B2182B', 'Low': '#2166AC'}

# ==============================================================================
# Load Data
# ==============================================================================
print("Loading BayesPrism epithelial expression...")
epi = pd.read_csv(TCGA_BAYESPRISM_EPI, sep='\t', index_col=0)
if 'CEACAM5' in epi.index:
    epi = epi.T
print(f"  Epithelial expression: {epi.shape}")

print("Loading clinical data...")
clin = pd.read_csv(CLINICAL_FILE, sep='\t')
print(f"  Clinical records: {len(clin)}")

# Build OS time and event
clin['os_time'] = clin['days_to_death'].fillna(clin['days_to_last_follow_up'])
clin['os_event'] = (clin['vital_status'] == 'Dead').astype(int)
clin = clin.dropna(subset=['os_time'])
clin = clin[clin['os_time'] > 0]
clin = clin.set_index('submitter_id')
print(f"  After filtering: {len(clin)} patients with OS data")

# Align
common = sorted(set(epi.index) & set(clin.index))
print(f"  Common samples: {len(common)}")

# ==============================================================================
# Create Figure
# ==============================================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'font.size': FONT_SIZE,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

fig, axes = plt.subplots(1, 2, figsize=(FIG_WIDTH, FIG_HEIGHT))

for idx, gene in enumerate(['CEACAM5', 'CEACAM6']):
    ax = axes[idx]

    expr = np.log2(epi.loc[common, gene].values + 1)
    median_val = np.median(expr)
    high_mask = expr >= median_val
    low_mask = ~high_mask

    high_ids = [common[i] for i in range(len(common)) if high_mask[i]]
    low_ids = [common[i] for i in range(len(common)) if low_mask[i]]

    t_high = clin.loc[high_ids, 'os_time'].values / 30.44  # days to months
    e_high = clin.loc[high_ids, 'os_event'].values
    t_low = clin.loc[low_ids, 'os_time'].values / 30.44
    e_low = clin.loc[low_ids, 'os_event'].values

    # KM fitting
    kmf_high = KaplanMeierFitter()
    kmf_high.fit(t_high, event_observed=e_high, label=f'High (n={len(high_ids)})')
    kmf_high.plot_survival_function(ax=ax, color=COLORS['High'],
                                     linewidth=0.5 * SCALE, ci_show=False)

    kmf_low = KaplanMeierFitter()
    kmf_low.fit(t_low, event_observed=e_low, label=f'Low (n={len(low_ids)})')
    kmf_low.plot_survival_function(ax=ax, color=COLORS['Low'],
                                    linewidth=0.5 * SCALE, ci_show=False)

    # Log-rank test
    result = logrank_test(t_high, t_low, event_observed_A=e_high, event_observed_B=e_low)
    p_val = result.p_value

    p_str = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
    ax.text(0.95, 0.95, f'Log-rank\n{p_str}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=FONT_SIZE * 0.75,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray', linewidth=0.5))

    ax.set_title(f'{gene} (Epi. deconvolved)', fontsize=FONT_SIZE * 0.85)
    ax.set_xlabel('Time (months)', fontsize=FONT_SIZE * 0.85)
    ax.set_ylabel('Overall Survival' if idx == 0 else '', fontsize=FONT_SIZE * 0.85)
    ax.set_ylim(0, 1.05)
    ax.legend(loc='lower left', fontsize=FONT_SIZE * 0.65, frameon=True)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(0.5 * SCALE)
    ax.tick_params(width=0.5 * SCALE, length=3 * SCALE, labelsize=5 * SCALE)

    print(f"\n  {gene}: High={len(high_ids)}, Low={len(low_ids)}, P={p_val:.4f}")

fig.suptitle('TCGA-STAD: Overall Survival by Epithelial Expression (Median Split)',
             fontsize=FONT_SIZE * 0.8, y=1.02)
plt.tight_layout()

plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_FILE.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f"\nSaved: {OUTPUT_FILE}")
