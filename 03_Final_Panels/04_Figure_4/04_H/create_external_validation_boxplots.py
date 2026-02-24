#!/usr/bin/env python3
"""
Panel G: External Validation Boxplots
- Linghua: Primary_tumor vs Adjacent_normal only (exclude Peritoneal_carcinomatosis)
- Kumar: Primary_tumor vs Normal
Standard colors: Tumor (pink) vs Normal (teal)
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

# 4× scaling method
SCALE = 4
DPI = 300

# Target print size per boxplot: ~3.0 x 3.0 cm (taller, similar to Panel C)
FIGURE_SIZE = (3.0 * SCALE / 2.54, 3.0 * SCALE / 2.54)

# Font sizes
LABEL_FONTSIZE = 5 * SCALE
TICK_FONTSIZE = 5 * SCALE
TITLE_FONTSIZE = 5 * SCALE

# Colors for Tumor vs Normal comparison
COLOR_TUMOR = '#E377C2'   # Pink for Tumor
COLOR_NORMAL = '#17BECF'  # Teal for Normal

# Darker shades for median lines
MEDIAN_COLOR_TUMOR = '#B85A9A'   # Darker pink
MEDIAN_COLOR_NORMAL = '#0F8A99'  # Darker teal

def create_boxplot(data, title, output_path, y_label='Mac_3 Proportion'):
    """Create a single boxplot for T vs N comparison."""

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': TICK_FONTSIZE,
        'axes.linewidth': 1.0 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Create boxplot - Tumor on left, Normal on right
    sns.boxplot(
        data=data,
        x='Group',
        y='Mac3_Proportion',
        order=['T', 'N'],
        palette=[COLOR_TUMOR, COLOR_NORMAL],
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
            ax.lines[4].set_color(MEDIAN_COLOR_TUMOR)  # First median
        elif i == 1:
            ax.lines[9].set_color(MEDIAN_COLOR_NORMAL)  # Second median

    # Statistical test
    tumor_vals = data[data['Group'] == 'T']['Mac3_Proportion'].values
    normal_vals = data[data['Group'] == 'N']['Mac3_Proportion'].values

    stat, p_value = mannwhitneyu(tumor_vals, normal_vals, alternative='two-sided')

    # Add significance bracket
    y_max = data['Mac3_Proportion'].max()
    y_range = data['Mac3_Proportion'].max() - data['Mac3_Proportion'].min()
    if y_range == 0:
        y_range = 0.1
    bracket_height = y_max + y_range * 0.1
    bracket_top = bracket_height + y_range * 0.08

    ax.plot([0, 0, 1, 1], [bracket_height, bracket_top, bracket_top, bracket_height],
            lw=0.8 * SCALE, c='black')

    p_text = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
    ax.text(0.5, bracket_top + y_range * 0.02, p_text,
            ha='center', va='bottom', fontsize=TICK_FONTSIZE * 0.85)

    # Styling
    ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight='normal')
    ax.set_ylabel(y_label, fontsize=LABEL_FONTSIZE)
    ax.set_xlabel('')
    ax.tick_params(labelsize=TICK_FONTSIZE, width=1.0 * SCALE, length=4 * SCALE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0 * SCALE)
    ax.spines['bottom'].set_linewidth(1.0 * SCALE)

    ax.set_ylim(ax.get_ylim()[0], bracket_top + y_range * 0.25)

    plt.tight_layout()
    fig.savefig(output_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(output_path.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_path.name}")
    print(f"    T: n={len(tumor_vals)}, mean={tumor_vals.mean():.4f}")
    print(f"    N: n={len(normal_vals)}, mean={normal_vals.mean():.4f}")
    print(f"    p-value: {p_value:.4e}")

    plt.close()

def main():
    print("=" * 60)
    print("Panel G: External Validation Boxplots")
    print("=" * 60)

    script_dir = Path(__file__).parent

    # === Linghua Dataset ===
    print("\n[1/2] Linghua (GSE239676)...")
    df_linghua = pd.read_csv(FIG4_EXTERNAL_LINGHUA)

    # Filter to Primary_tumor vs Adjacent_normal only (exclude Peritoneal_carcinomatosis)
    df_linghua_filtered = df_linghua[df_linghua['Tissue'].isin(['Primary_tumor', 'Adjacent_normal'])].copy()
    df_linghua_filtered['Group'] = df_linghua_filtered['Tissue'].map({
        'Primary_tumor': 'T',
        'Adjacent_normal': 'N'
    })

    print(f"  Original samples: {len(df_linghua)}")
    print(f"  After filtering (T vs N only): {len(df_linghua_filtered)}")

    output_linghua = script_dir / 'mac3_boxplot_linghua.png'
    create_boxplot(df_linghua_filtered, 'GSE239676', output_linghua)

    # === Kumar Dataset ===
    print("\n[2/2] Kumar (GSE183904)...")
    df_kumar = pd.read_csv(FIG4_EXTERNAL_KUMAR)

    # Already has Normal vs Primary_tumor
    df_kumar['Group'] = df_kumar['Tissue'].map({
        'Primary_tumor': 'T',
        'Normal': 'N'
    })

    output_kumar = script_dir / 'mac3_boxplot_kumar.png'
    create_boxplot(df_kumar, 'GSE183904', output_kumar)

    print("\n" + "=" * 60)
    print("Panel G complete")
    print("=" * 60)

if __name__ == "__main__":
    main()
