#!/usr/bin/env python3
"""
Panel G: Distance to Immune Boxplot (CEACAM-high vs CEACAM-low)
Moved from old 03_F. Scaling violations fixed.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA

DPI = 300
PANEL_WIDTH_CM = 5.0 * 4
PANEL_HEIGHT_CM = 5.5 * 4
CM_TO_INCH = 1 / 2.54
SCALE = 4

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

CEACAM_LOW_COLOR = '#1f77b4'
CEACAM_HIGH_COLOR = '#d62728'


def create_paired_boxplot(sample_data, col, title, ylabel, output_path, y_max_limit=None):
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 9 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    fig_width = PANEL_WIDTH_CM * CM_TO_INCH
    fig_height = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    low_data = sample_data[sample_data['group'] == 'CEACAM-low'][col].values
    high_data = sample_data[sample_data['group'] == 'CEACAM-high'][col].values

    stat, pval = stats.wilcoxon(high_data, low_data, alternative='greater')

    bp = ax.boxplot([low_data, high_data], positions=[1, 2], widths=0.5, patch_artist=True,
                    boxprops=dict(linewidth=1.0),
                    whiskerprops=dict(linewidth=1.0),
                    capprops=dict(linewidth=1.0),
                    medianprops=dict(linewidth=2.0))
    bp['boxes'][0].set_facecolor(CEACAM_LOW_COLOR)
    bp['boxes'][1].set_facecolor(CEACAM_HIGH_COLOR)
    for box in bp['boxes']:
        box.set_alpha(0.7)
        box.set_edgecolor('black')
    for median in bp['medians']:
        median.set_color('black')

    for i in range(len(low_data)):
        ax.plot([1, 2], [low_data[i], high_data[i]], 'k-', alpha=0.3, linewidth=0.8)

    ax.scatter([1]*len(low_data), low_data, color=CEACAM_LOW_COLOR, s=40, zorder=3,
               edgecolor='black', linewidth=0.5)
    ax.scatter([2]*len(high_data), high_data, color=CEACAM_HIGH_COLOR, s=40, zorder=3,
               edgecolor='black', linewidth=0.5)

    data_y_max = max(max(low_data), max(high_data))
    y_range = data_y_max - min(min(low_data), min(high_data))
    bracket_y = data_y_max + 0.08 * y_range

    ax.plot([1, 1, 2, 2], [bracket_y, bracket_y + 0.02*y_range, bracket_y + 0.02*y_range, bracket_y],
            color='black', linewidth=1)

    if pval <= 0.0001:
        pval_text = '****'
    elif pval <= 0.001:
        pval_text = '***'
    elif pval <= 0.01:
        pval_text = '**'
    elif pval <= 0.05:
        pval_text = '*'
    else:
        pval_text = 'ns'
    ax.text(1.5, bracket_y + 0.04*y_range, pval_text, ha='center', va='bottom', fontsize=9*SCALE)

    ax.set_title(title, fontsize=9*SCALE, fontweight='normal')
    ax.set_ylabel(ylabel, fontsize=8*SCALE)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['CEACAM-\nlow', 'CEACAM-\nhigh'], fontsize=9*SCALE)
    ax.tick_params(axis='both', labelsize=8*SCALE, width=1.0, length=4)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    if y_max_limit is not None:
        ax.set_ylim(0, y_max_limit)

    n_pairs = len(low_data)
    ax.text(0.98, 0.02, f'n = {n_pairs} samples\nWilcoxon signed-rank (one-tailed)',
            transform=ax.transAxes, fontsize=7*SCALE, va='bottom', ha='right', color='#555555')

    plt.tight_layout()

    plt.savefig(output_path, dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.svg'), format='svg', facecolor='white', bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def main():
    print("=" * 60)
    print("Panel G: Distance to Immune Boxplot")
    print("=" * 60)

    print("\n[1/3] Loading spot data...")
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    print(f"  Total spots (all samples): {len(df):,}")

    print("\n[2/3] Computing sample-level means...")
    sample_metrics = []
    for sample in df['sample'].unique():
        sample_df = df[df['sample'] == sample]
        for group in ['CEACAM-high', 'CEACAM-low']:
            group_df = sample_df[sample_df['CEACAM_group'] == group]
            if len(group_df) > 0:
                sample_metrics.append({
                    'sample': sample,
                    'group': group,
                    'distance_to_immune': group_df['distance_to_immune'].mean(),
                })

    sample_data = pd.DataFrame(sample_metrics)
    valid_samples = sample_data.groupby('sample').filter(lambda x: len(x) == 2)['sample'].unique()
    sample_data = sample_data[sample_data['sample'].isin(valid_samples)]
    sample_data = sample_data.sort_values(['sample', 'group'])
    print(f"  Samples with paired data: {len(valid_samples)}")

    print("\n[3/3] Creating boxplot...")
    create_paired_boxplot(sample_data, 'distance_to_immune', 'Distance to Immune',
                         'Distance (a.u.)',
                         os.path.join(BASE_DIR, 'ceacam_spatial_boxplot.png'), y_max_limit=1200)

    print("\nDone!")


if __name__ == '__main__':
    main()
