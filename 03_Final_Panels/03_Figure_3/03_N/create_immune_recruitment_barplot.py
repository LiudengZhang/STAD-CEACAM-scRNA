#!/usr/bin/env python3
"""
Panel K: CEACAM vs Immune Infiltration — Region-Level Bar Chart
Moved from old 03_G. Scaling violations fixed. Wider panel (2-col in assembly).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_REGION_COMPARISON

DPI = 300
PANEL_WIDTH_CM = 7.0 * 4   # 28 cm electronic → 7 cm print (wider for 2-col)
PANEL_HEIGHT_CM = 5.5 * 4   # 22 cm electronic → 5.5 cm print
CM_TO_INCH = 1 / 2.54
SCALE = 4

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("=" * 60)
    print("Panel K: Immune Recruitment by CEACAM Region")
    print("=" * 60)

    print("\n[1/3] Loading region-level results...")
    results_df = pd.read_csv(SPATIAL_REGION_COMPARISON)
    results_df = results_df.sort_values('difference', ascending=False)
    print(f"  Loaded results for {len(results_df)} immune cell types")

    print("\n[2/3] Creating visualization...")
    fig_width = PANEL_WIDTH_CM * CM_TO_INCH
    fig_height = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    SHORT_NAMES = {
        'Monocytes/Macrophages': 'Mo/Mac',
        'CD4+ T cells': 'CD4+ T',
        'CD8+ T cells': 'CD8+ T',
        'Plasma cells': 'Plasma',
    }

    df_sorted = results_df.sort_values('difference', ascending=False)
    df_sorted['label'] = df_sorted['immune_cell'].map(lambda x: SHORT_NAMES.get(x, x))
    x_pos = np.arange(len(df_sorted))
    colors = ['#d62728' if d > 0 else '#1f77b4' for d in df_sorted['difference']]

    bars = ax.bar(x_pos, df_sorted['difference'], color=colors, edgecolor='black',
                  alpha=0.8, width=0.7, linewidth=1.0)

    ax.axhline(y=0, color='black', linestyle='-', linewidth=1.0)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df_sorted['label'], fontsize=7 * SCALE, rotation=45, ha='right')
    ax.set_ylabel('Difference in Proportion\n(CEACAM-high - CEACAM-low)', fontsize=8 * SCALE)
    ax.set_title('Immune Recruitment by CEACAM Region', fontsize=9 * SCALE, fontweight='normal')
    ax.tick_params(axis='both', labelsize=7 * SCALE, width=1.0, length=4)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    y_min = df_sorted['difference'].min() - 0.01
    y_max = df_sorted['difference'].max() + 0.02
    ax.set_ylim(y_min, y_max)

    for i, (_, row) in enumerate(df_sorted.iterrows()):
        pval = row['pvalue']
        sig = ''
        if pval < 0.001:
            sig = '***'
        elif pval < 0.01:
            sig = '**'
        elif pval < 0.05:
            sig = '*'

        diff = row['difference']
        if diff >= 0:
            y_pos_text = diff + 0.002
            va = 'bottom'
        else:
            y_pos_text = diff - 0.002
            va = 'top'
        ax.text(i, y_pos_text, sig, ha='center', va=va, fontsize=8 * SCALE, fontweight='bold')

    plt.tight_layout()

    print("\n[3/3] Saving figure...")
    out = os.path.join(BASE_DIR, 'ceacam_immune_region_summary.png')
    plt.savefig(out, dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.savefig(out.replace('.png', '.svg'), format='svg', facecolor='white', bbox_inches='tight')
    plt.savefig(out.replace('.png', '.pdf'), dpi=DPI, facecolor='white', bbox_inches='tight')
    print(f"  Saved: {out}")
    plt.close()

    print("\nDone!")


if __name__ == '__main__':
    main()
