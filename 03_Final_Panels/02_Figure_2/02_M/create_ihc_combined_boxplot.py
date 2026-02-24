#!/usr/bin/env python3
"""
Panel 2M: Combined CEACAM5+6 IHC staining boxplot (R vs NR)
One-tailed Mann-Whitney U test (NR > R), n=4 R, n=4 NR (all pre-treatment)
Data from color deconvolution (Ruifrok & Johnston 2001) via skimage rgb2hed.
DAB OD threshold = 0.02. Source: 02_Preparation_for_Panels/IHC/quantify_ceacam_ihc.py
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import IHC_COLOR_DECONV_CSV

# 4x scaling (scale canvas + fonts only, NOT linewidths/markers)
SCALE = 4
DPI = 300
CM_TO_INCH = 1 / 2.54

# Panel: 3.2cm x 3.0cm print (match J/K sizing)
PANEL_WIDTH_CM = 3.2 * SCALE
PANEL_HEIGHT_CM = 3.0 * SCALE

OUTPUT_DIR = Path(__file__).parent

# Colors matching Figure 2 panels J/K
COLOR_R = '#0072B2'
COLOR_NR = '#D55E00'
MEDIAN_R = '#005689'
MEDIAN_NR = '#A34700'

# Read IHC data from color deconvolution results
df = pd.read_csv(IHC_COLOR_DECONV_CSV)

# Pivot: combined CEACAM5 + CEACAM6 per patient
pivot = df.pivot_table(index=['patient', 'group'], columns='marker',
                       values='staining_pct').reset_index()
pivot['combined'] = pivot['CEACAM5'] + pivot['CEACAM6']

R_VALS = pivot[pivot['group'] == 'R']['combined'].values
NR_VALS = pivot[pivot['group'] == 'NR']['combined'].values


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    fig, ax = plt.subplots(
        figsize=(PANEL_WIDTH_CM * CM_TO_INCH, PANEL_HEIGHT_CM * CM_TO_INCH)
    )

    # One-tailed Mann-Whitney U (NR > R)
    stat, pval = stats.mannwhitneyu(NR_VALS, R_VALS, alternative='greater')
    print(f"R  (n={len(R_VALS)}): {R_VALS} -> mean={R_VALS.mean():.2f}%")
    print(f"NR (n={len(NR_VALS)}): {NR_VALS} -> mean={NR_VALS.mean():.2f}%")
    print(f"Fold: {NR_VALS.mean() / R_VALS.mean():.2f}x")
    print(f"Mann-Whitney U (one-tailed, NR > R): U={stat:.1f}, P={pval:.4f}")

    # Boxplot
    bp = ax.boxplot(
        [R_VALS, NR_VALS],
        positions=[1, 2],
        widths=0.6,
        patch_artist=True,
        boxprops=dict(linewidth=0.5),
        whiskerprops=dict(color='black', linewidth=0.5),
        capprops=dict(color='black', linewidth=0.5),
        flierprops=dict(markersize=0),  # hide outliers, show strip points instead
    )

    bp['boxes'][0].set_facecolor(COLOR_R)
    bp['boxes'][0].set_edgecolor('black')
    bp['boxes'][1].set_facecolor(COLOR_NR)
    bp['boxes'][1].set_edgecolor('black')
    bp['medians'][0].set_color(MEDIAN_R)
    bp['medians'][0].set_linewidth(0.8)
    bp['medians'][1].set_color(MEDIAN_NR)
    bp['medians'][1].set_linewidth(0.8)

    # Strip points (essential with small n)
    rng = np.random.default_rng(42)
    jitter_r = rng.uniform(-0.12, 0.12, len(R_VALS))
    jitter_nr = rng.uniform(-0.12, 0.12, len(NR_VALS))

    ax.scatter(1 + jitter_r, R_VALS, c=COLOR_R, s=20, zorder=5,
               edgecolors='black', linewidths=0.3)
    ax.scatter(2 + jitter_nr, NR_VALS, c=COLOR_NR, s=20, zorder=5,
               edgecolors='black', linewidths=0.3)

    # Significance bracket
    y_max = max(np.max(R_VALS), np.max(NR_VALS))
    y_bracket = y_max * 1.15

    ax.plot([1, 1, 2, 2],
            [y_bracket, y_bracket * 1.05, y_bracket * 1.05, y_bracket],
            'k-', linewidth=0.5)

    p_text = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'

    ax.text(1.5, y_bracket * 1.08, p_text, ha='center', va='bottom',
            fontsize=6 * SCALE)

    # Labels
    ax.set_title('IHC (CEACAM5+6)', fontsize=7 * SCALE, fontweight='normal')
    ax.set_ylabel('Staining area (%)', fontsize=6 * SCALE)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f'R\n(n={len(R_VALS)})', f'NR\n(n={len(NR_VALS)})'],
                        fontsize=6 * SCALE)
    ax.set_ylim(0, y_max * 1.40)

    # Spines (no scaling on linewidths!)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.tick_params(axis='both', labelsize=6 * SCALE, width=0.5, length=4)

    fig.subplots_adjust(left=0.30, right=0.95, top=0.88, bottom=0.18)

    # Save
    for ext in ['png', 'svg']:
        out_path = OUTPUT_DIR / f"ihc_combined_boxplot.{ext}"
        fig.savefig(out_path, dpi=DPI, facecolor='white')
    plt.close()
    print(f"\nSaved to {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
