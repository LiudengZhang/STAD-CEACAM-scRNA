#!/usr/bin/env python3
"""
Panel 2L: CEACAM6 boxplot from PRJEB25780 (BayesPrism epithelial-deconvolved)
Target-size approach: panel created at exact assembly slot dimensions.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TIGER_BAYESPRISM_EPI, TIGER_META

# Target slot from assembly layout (mm)
# Row 3 stacked L/M: w = 180 * 0.15 = 27mm, h = (50 - 1) / 2 = 24.5mm
PANEL_W_MM = 27.0
PANEL_H_MM = 24.5
MM_TO_INCH = 1 / 25.4
DPI = 300

OUTPUT_DIR = Path(__file__).parent

# Standard R vs NR colors
COLOR_R = '#2166AC'
COLOR_NR = '#B2182B'


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6,
        'axes.labelsize': 6,
        'xtick.labelsize': 5,
        'ytick.labelsize': 5,
        'axes.linewidth': 0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 2,
        'ytick.major.size': 2,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    # Load BayesPrism deconvolved epithelial expression
    print("Loading BayesPrism epithelial expression...")
    epi_expr = pd.read_csv(TIGER_BAYESPRISM_EPI, sep='\t', index_col=0)
    meta_df = pd.read_csv(TIGER_META, sep='\t')
    meta_df = meta_df[meta_df['Treatment'] != 'Normal']

    # Match samples
    common = [s for s in meta_df['sample_id'] if s in epi_expr.index]
    meta_df = meta_df[meta_df['sample_id'].isin(common)].set_index('sample_id')
    epi_expr = epi_expr.loc[common]

    gene = 'CEACAM6'
    r_vals = np.log2(epi_expr.loc[meta_df['response_NR'] == 'R', gene].values + 1)
    nr_vals = np.log2(epi_expr.loc[meta_df['response_NR'] == 'N', gene].values + 1)

    print(f"R (n={len(r_vals)}): mean={np.mean(r_vals):.3f}")
    print(f"NR (n={len(nr_vals)}): mean={np.mean(nr_vals):.3f}")

    stat, pval = stats.mannwhitneyu(nr_vals, r_vals, alternative='greater')
    print(f"Mann-Whitney U (NR > R): U={stat:.0f}, P={pval:.4f}")

    # Plot
    fig, ax = plt.subplots(figsize=(PANEL_W_MM * MM_TO_INCH, PANEL_H_MM * MM_TO_INCH))

    data = [r_vals, nr_vals]
    positions = [1, 2]
    colors = [COLOR_R, COLOR_NR]

    # Violin
    vp = ax.violinplot(data, positions=positions, widths=0.7,
                       showmeans=False, showmedians=False, showextrema=False)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors[i])
        body.set_edgecolor('black')
        body.set_linewidth(0.5)
        body.set_alpha(0.7)

    # Boxplot
    bp = ax.boxplot(data, positions=positions, widths=0.15, patch_artist=True,
                    showfliers=False,
                    boxprops=dict(facecolor='white', linewidth=0.5),
                    whiskerprops=dict(color='black', linewidth=0.5),
                    capprops=dict(color='black', linewidth=0.5),
                    medianprops=dict(color='black', linewidth=0.8))

    # Jittered points
    np.random.seed(42)
    for i, (pos, vals, color) in enumerate(zip(positions, data, colors)):
        jitter = np.random.uniform(-0.08, 0.08, len(vals))
        ax.scatter(pos + jitter, vals, c=color, s=6, alpha=0.8,
                   edgecolors='white', linewidths=0.3, zorder=3)

    # Significance bracket
    y_max = max(np.max(r_vals), np.max(nr_vals))
    y_bracket = y_max * 1.1
    ax.plot([1, 1, 2, 2],
            [y_bracket, y_bracket * 1.03, y_bracket * 1.03, y_bracket],
            'k-', linewidth=0.5)
    p_text = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
    ax.text(1.5, y_bracket * 1.05, p_text, ha='center', va='bottom', fontsize=4)

    # Labels
    ax.set_title(r'$\it{CEACAM6}$ (Epi)', fontsize=5, pad=2)
    ax.set_ylabel('Expression (log2)', fontsize=5)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['R', 'NR'], fontsize=5)
    ax.set_ylim(0, y_max * 1.35)

    plt.tight_layout(pad=0.3)

    # Save
    for ext in ['svg', 'png', 'pdf']:
        out = OUTPUT_DIR / f"ceacam6_prjeb25780_boxplot.{ext}"
        fig.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
        print(f"Saved: {out}")
    plt.close()


if __name__ == '__main__':
    main()
