#!/usr/bin/env python3
"""
Panel C: IM-MoMac Post-Treatment Response Boxplot
Fixed version: Blue for R (left), Orange for NR (right)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

# Configuration - 4x scaling for Nature Cancer
# Target print size: ~3.0 x 3.2 cm (smaller to fit Row 2 height)
SCALE = 4
FIGURE_SIZE = (3.0 * SCALE / 2.54, 3.2 * SCALE / 2.54)  # 4x print size
DPI = 300
FONT_SIZE = 5 * SCALE  # 5pt at print

# Standard colors
COLOR_R = '#0072B2'    # Blue for Responder
COLOR_NR = '#D55E00'   # Vermillion for Non-responder
PALETTE = [COLOR_R, COLOR_NR]  # R first (left), NR second (right)

# Darker shades for median lines
MEDIAN_COLOR_R = '#005689'   # Darker blue
MEDIAN_COLOR_NR = '#A34700'  # Darker vermillion

def main():
    print("=" * 60)
    print("Panel C: IM-MoMac Response Boxplot")
    print("=" * 60)

    script_dir = Path(__file__).parent

    # Load data from central config paths
    print("\nLoading data...")
    module_props = pd.read_csv(FIG4_MODULE_PROPORTIONS, index_col=0)
    clinical_data = pd.read_csv(FIG4_CLINICAL_METADATA)

    if 'Sample ID' in clinical_data.columns:
        clinical_data = clinical_data.set_index('Sample ID')

    merged = module_props.join(clinical_data, how='inner')

    # Filter for post-treatment stomach samples
    filtered = merged[
        (merged['Sample site'] == 'Stomach') &
        (merged['Treatment phase'] == 'Post') &
        (merged['stomach_post_grouping'].isin(['No-response', 'Responsed']))
    ].copy()

    print(f"  Responsed (R): {len(filtered[filtered['stomach_post_grouping'] == 'Responsed'])}")
    print(f"  No-response (NR): {len(filtered[filtered['stomach_post_grouping'] == 'No-response'])}")

    # Statistical test (Mann-Whitney U)
    group_r = filtered[filtered['stomach_post_grouping'] == 'Responsed']['Module_2'].values
    group_nr = filtered[filtered['stomach_post_grouping'] == 'No-response']['Module_2'].values

    stat, p_value = mannwhitneyu(group_r, group_nr, alternative='two-sided')
    print(f"  Mann-Whitney p-value: {p_value:.4f}")

    # Create figure
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': FONT_SIZE,
        'axes.linewidth': 1.0 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Create boxplot with R on left, NR on right
    sns.boxplot(
        data=filtered,
        x='stomach_post_grouping',
        y='Module_2',
        order=['Responsed', 'No-response'],  # R first (left)
        palette=PALETTE,
        ax=ax,
        boxprops=dict(edgecolor='black', linewidth=1.0 * SCALE),
        medianprops=dict(color='black', linewidth=1.5 * SCALE),  # Will be overridden
        whiskerprops=dict(linewidth=1.0 * SCALE),
        capprops=dict(linewidth=1.0 * SCALE),
        flierprops=dict(marker='o', markerfacecolor='white', markeredgecolor='black',
                        markersize=4 * SCALE, markeredgewidth=1 * SCALE),
        width=0.6
    )

    # Color median lines to match box colors (darker shades)
    for i, artist in enumerate(ax.patches):
        if i == 0:
            ax.lines[4].set_color(MEDIAN_COLOR_R)  # First median
        elif i == 1:
            ax.lines[9].set_color(MEDIAN_COLOR_NR)  # Second median

    # Update x-tick labels
    ax.set_xticklabels(['R', 'NR'], fontsize=FONT_SIZE)

    # Add significance bracket
    y_max = filtered['Module_2'].max()
    y_range = filtered['Module_2'].max() - filtered['Module_2'].min()
    bracket_height = y_max + y_range * 0.1
    bracket_top = bracket_height + y_range * 0.05

    ax.plot([0, 0, 1, 1], [bracket_height, bracket_top, bracket_top, bracket_height],
            lw=0.8 * SCALE, c='black')
    p_text = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
    ax.text(0.5, bracket_top + y_range * 0.02, p_text,
            ha='center', va='bottom', fontsize=FONT_SIZE * 0.85)

    # Styling
    ax.set_ylabel('Module Proportion', fontsize=FONT_SIZE)
    ax.set_xlabel('')
    ax.tick_params(labelsize=FONT_SIZE, width=1.0 * SCALE, length=4 * SCALE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0 * SCALE)
    ax.spines['bottom'].set_linewidth(1.0 * SCALE)

    ax.set_ylim(ax.get_ylim()[0], bracket_top + y_range * 0.15)

    plt.tight_layout()

    # Save
    output_path = script_dir / 'module2_post_response_boxplot_k5.png'
    fig.savefig(output_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(output_path.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_path}")

    plt.close()

if __name__ == "__main__":
    main()
