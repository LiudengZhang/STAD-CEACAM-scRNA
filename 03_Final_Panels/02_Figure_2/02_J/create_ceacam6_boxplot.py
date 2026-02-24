#!/usr/bin/env python3
"""
Create Panel 2B: CEACAM6 boxplot comparing Responders vs Non-Responders
Sample-level aggregation (merged Pre+Post) with Mann-Whitney U test
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import os

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

# Nature Cancer specifications - 4× scaling method
# Electronic size = Print target × 4, then scale down for crisp text
DPI = 300
PANEL_WIDTH_CM = 3.2 * 4   # 12.8 cm electronic → 3.2 cm print
PANEL_HEIGHT_CM = 3.0 * 4  # 12 cm electronic → 3.0 cm print (increased for title)
CM_TO_INCH = 1 / 2.54
SCALE = 4  # Scaling factor for fonts and line widths

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# DATA_PATH - now using EPITHELIAL_H5AD from central config
OUTPUT_DIR = BASE_DIR
USE_RAW = True  # Use raw counts instead of normalized

# Colors: Blue for R, Red for NR
COLORS = {
    'Responsed': '#0072B2',      # Blue for responders
    'No-response': '#D55E00',    # Vermillion/red for non-responders
}
MEDIAN_COLORS = {
    'Responsed': '#005689',      # Darker blue median line
    'No-response': '#A34700',    # Darker red median line
}

def main():
    # Set up matplotlib - 4× scaled fonts
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,  # 28pt electronic → 7pt print
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    # Load data
    print("Loading data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    print(f"Loaded {adata.n_obs} cells")

    # Filter for PRE-TREATMENT only with valid response groups
    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
    adata_filtered = adata[pre_mask & valid_mask].copy()
    print(f"Pre-treatment cells with response status: {adata_filtered.n_obs}")

    # Get CEACAM6 expression - prefer raw counts for better statistical power
    if USE_RAW and adata_filtered.raw is not None and 'CEACAM6' in adata_filtered.raw.var_names:
        print("Using raw counts for CEACAM6 expression")
        ceacam6_idx = adata_filtered.raw.var_names.get_loc('CEACAM6')
        if hasattr(adata_filtered.raw.X, 'toarray'):
            ceacam6_expr = adata_filtered.raw.X[:, ceacam6_idx].toarray().flatten()
        else:
            ceacam6_expr = adata_filtered.raw.X[:, ceacam6_idx].flatten()
    elif 'CEACAM6' in adata_filtered.var_names:
        print("Using normalized data for CEACAM6 expression")
        ceacam6_idx = adata_filtered.var_names.get_loc('CEACAM6')
        if hasattr(adata_filtered.X, 'toarray'):
            ceacam6_expr = adata_filtered.X[:, ceacam6_idx].toarray().flatten()
        else:
            ceacam6_expr = adata_filtered.X[:, ceacam6_idx].flatten()
    else:
        print("Error: CEACAM6 not found")
        return

    # Add to obs for aggregation
    adata_filtered.obs['CEACAM6'] = ceacam6_expr

    # Calculate sample-level mean expression
    sample_means = adata_filtered.obs.groupby('sample', observed=True).agg({
        'CEACAM6': 'mean',
        'stomach_pre_grouping': 'first'
    }).reset_index()
    print(f"\nSample-level data ({len(sample_means)} samples):")
    print(sample_means.to_string())

    # Separate by response group
    responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'Responsed']['CEACAM6'].values
    non_responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'No-response']['CEACAM6'].values

    print(f"\nResponders (n={len(responder_vals)}): mean={np.mean(responder_vals):.3f}")
    print(f"Non-Responders (n={len(non_responder_vals)}): mean={np.mean(non_responder_vals):.3f}")

    # Mann-Whitney U test
    stat, pval = stats.mannwhitneyu(non_responder_vals, responder_vals, alternative='greater')
    print(f"\nOne-tailed Mann-Whitney U test (NR > R): U={stat:.2f}, p={pval:.4f}")

    # Create figure with exact dimensions
    fig_width = PANEL_WIDTH_CM * CM_TO_INCH
    fig_height = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Prepare data for boxplot
    data = [responder_vals, non_responder_vals]
    positions = [1, 2]

    # Create boxplot - 4× scaled line widths
    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showfliers=True,
        boxprops=dict(linewidth=0.5),
        whiskerprops=dict(color='black', linewidth=0.5),
        capprops=dict(color='black', linewidth=0.5),
        flierprops=dict(marker='o', markerfacecolor='white', markersize=4,
                       linestyle='none', markeredgecolor='black', markeredgewidth=0.5)
    )

    # Style boxes
    bp['boxes'][0].set_facecolor(COLORS['Responsed'])
    bp['boxes'][0].set_edgecolor('black')
    bp['boxes'][1].set_facecolor(COLORS['No-response'])
    bp['boxes'][1].set_edgecolor('black')

    # Color median lines - 4× scaled
    bp['medians'][0].set_color(MEDIAN_COLORS['Responsed'])
    bp['medians'][0].set_linewidth(0.8)
    bp['medians'][1].set_color(MEDIAN_COLORS['No-response'])
    bp['medians'][1].set_linewidth(0.8)

    # Add significance bracket - 4× scaled
    y_max = max(np.max(responder_vals), np.max(non_responder_vals))
    y_bracket = y_max * 1.15
    ax.plot([1, 1, 2, 2], [y_bracket, y_bracket*1.05, y_bracket*1.05, y_bracket],
            'k-', linewidth=0.5)
    pval_text = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
    ax.text(1.5, y_bracket*1.08, pval_text, ha='center', va='bottom', fontsize=6 * SCALE)

    # Labels - 4× scaled fonts
    ax.set_title('CEACAM6', fontsize=7 * SCALE, fontweight='normal')
    ax.set_ylabel('Expression', fontsize=6 * SCALE)
    ax.set_xticks(positions)
    ax.set_xticklabels(['R', 'NR'], fontsize=6 * SCALE)
    ax.set_ylim(0, y_max * 1.40)

    # Styling - 4× scaled (no grid per user preference)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.tick_params(axis='both', labelsize=6 * SCALE, width=0.5, length=4)

    # Adjust layout - top=0.88 to leave room for title
    fig.subplots_adjust(left=0.28, right=0.95, top=0.88, bottom=0.15)

    # Save
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output_path = os.path.join(OUTPUT_DIR, "ceacam6_pre_boxplot.png")
    plt.savefig(output_path, dpi=DPI, facecolor='white')
    plt.savefig(output_path.replace('.png', '.svg'), dpi=DPI, facecolor='white')
    print(f"\nSaved: {output_path}")

    output_pdf = os.path.join(OUTPUT_DIR, "ceacam6_pre_boxplot.pdf")
    plt.savefig(output_pdf, dpi=DPI, facecolor='white')
    print(f"Saved: {output_pdf}")

    plt.close()

if __name__ == '__main__':
    main()
