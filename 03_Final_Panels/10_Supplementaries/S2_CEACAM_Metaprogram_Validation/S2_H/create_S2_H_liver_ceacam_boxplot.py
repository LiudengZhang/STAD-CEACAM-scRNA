#!/usr/bin/env python3
"""
S2_H: CEACAM5/6 expression in liver epithelial cells, R vs NR.
Two side-by-side boxplots. Style matches Figure 2 02_K (no jitter, thin lines).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from pathlib import Path

SCALE = 4
DPI = 300
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 5.5 * 4   # two subplots side by side
PANEL_HEIGHT_CM = 3.0 * 4  # match Fig 2 height

# Colors matching Fig 2: Blue for R, Orange for NR
COLORS = {'R': '#0072B2', 'NR': '#D55E00'}
MEDIAN_COLORS = {'R': '#005689', 'NR': '#A34700'}

DATA_DIR = Path(__file__).parent.parent / 'data_prep'
OUTPUT_DIR = Path(__file__).parent


def main():
    print("S2_H: Liver CEACAM5/6 Expression R vs NR")

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    df = pd.read_csv(DATA_DIR / 'liver_stomach_ceacam_per_sample.csv')
    liver = df[df['site'] == 'Liver'].copy()

    fig, axes = plt.subplots(1, 2, figsize=(
        PANEL_WIDTH_CM * CM_TO_INCH,
        PANEL_HEIGHT_CM * CM_TO_INCH
    ))

    for idx, gene in enumerate(['CEACAM5', 'CEACAM6']):
        ax = axes[idx]
        r_vals = liver[liver['response'] == 'R'][gene].dropna().values
        nr_vals = liver[liver['response'] == 'NR'][gene].dropna().values

        # Boxplot — unified structural widths (0.5 * SCALE)
        data = [r_vals, nr_vals]
        bp = ax.boxplot(data, positions=[1, 2], widths=0.6, patch_artist=True,
                        boxprops=dict(linewidth=0.5 * SCALE),
                        whiskerprops=dict(color='black', linewidth=0.5 * SCALE),
                        capprops=dict(color='black', linewidth=0.5 * SCALE),
                        flierprops=dict(marker='o', markerfacecolor='white', markersize=4,
                                       markeredgecolor='black', markeredgewidth=0.5 * SCALE),
                        medianprops=dict(linewidth=0.8 * SCALE))

        bp['boxes'][0].set_facecolor(COLORS['R'])
        bp['boxes'][1].set_facecolor(COLORS['NR'])
        bp['medians'][0].set_color(MEDIAN_COLORS['R'])
        bp['medians'][1].set_color(MEDIAN_COLORS['NR'])

        # Significance bracket
        all_vals = np.concatenate([r_vals, nr_vals])
        y_max = np.nanmax(all_vals)
        y_bracket = y_max * 1.15

        if len(r_vals) >= 2 and len(nr_vals) >= 2:
            stat, p_val = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        else:
            p_val = np.nan

        ax.plot([1, 1, 2, 2], [y_bracket, y_bracket * 1.05, y_bracket * 1.05, y_bracket],
                'k-', linewidth=0.5 * SCALE)

        if np.isnan(p_val):
            p_text = 'N/A'
        else:
            p_text = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
        ax.text(1.5, y_bracket * 1.08, p_text, ha='center', va='bottom', fontsize=6 * SCALE)

        ax.set_title(f'{gene}', fontsize=6.5 * SCALE, fontweight='normal', style='italic')
        ax.set_ylabel('Expression' if idx == 0 else '', fontsize=6 * SCALE)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(['R', 'NR'], fontsize=6 * SCALE)
        ax.set_ylim(0, y_max * 1.40)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.5 * SCALE)
        ax.spines['bottom'].set_linewidth(0.5 * SCALE)
        ax.tick_params(axis='both', labelsize=5 * SCALE, width=0.5 * SCALE, length=3 * SCALE)

    fig.suptitle('Liver metastasis epithelial cells', fontsize=6.5 * SCALE,
                 fontweight='normal', y=0.98)
    fig.subplots_adjust(left=0.18, right=0.95, top=0.85, bottom=0.15, wspace=0.35)

    out = OUTPUT_DIR / 'panel_S2_H'
    fig.savefig(str(out) + '.png', dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(str(out) + '.svg', format='svg', bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}.png / .svg")
    plt.close()


if __name__ == '__main__':
    main()
