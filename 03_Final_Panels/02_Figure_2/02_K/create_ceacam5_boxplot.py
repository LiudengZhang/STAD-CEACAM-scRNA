#!/usr/bin/env python3
"""
Create CEACAM5 boxplot comparing Pre-treatment Responders vs Non-Responders
Sample-level aggregation with Mann-Whitney U test (raw counts)
4× scaling method for crisp text rendering
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

# Nature Cancer specifications - 4× scaling method
DPI = 300
PANEL_WIDTH_CM = 3.2 * 4   # 12.8 cm electronic → 3.2 cm print
PANEL_HEIGHT_CM = 3.0 * 4  # 12 cm electronic → 3.0 cm print (increased for title)
CM_TO_INCH = 1 / 2.54
SCALE = 4

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# DATA_PATH - now using EPITHELIAL_H5AD from central config
OUTPUT_DIR = BASE_DIR

# Colors: Blue for R, Red for NR
COLORS = {
    'Responsed': '#0072B2',
    'No-response': '#D55E00',
}
MEDIAN_COLORS = {
    'Responsed': '#005689',
    'No-response': '#A34700',
}

def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("Loading data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)

    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
    adata_filtered = adata[pre_mask & valid_mask].copy()
    print(f"Pre-treatment cells: {adata_filtered.n_obs}")

    idx = adata_filtered.raw.var_names.get_loc('CEACAM5')
    expr = adata_filtered.raw.X[:, idx].toarray().flatten()
    adata_filtered.obs['CEACAM5'] = expr

    sample_means = adata_filtered.obs.groupby('sample', observed=True).agg({
        'CEACAM5': 'mean',
        'stomach_pre_grouping': 'first'
    }).reset_index()

    responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'Responsed']['CEACAM5'].values
    non_responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'No-response']['CEACAM5'].values

    stat, pval = stats.mannwhitneyu(non_responder_vals, responder_vals, alternative='greater')
    print(f"Responders (n={len(responder_vals)}): mean={np.mean(responder_vals):.3f}")
    print(f"Non-Responders (n={len(non_responder_vals)}): mean={np.mean(non_responder_vals):.3f}")
    print(f"Mann-Whitney p={pval:.4f}")

    fig, ax = plt.subplots(figsize=(PANEL_WIDTH_CM * CM_TO_INCH, PANEL_HEIGHT_CM * CM_TO_INCH))

    data = [responder_vals, non_responder_vals]
    bp = ax.boxplot(data, positions=[1, 2], widths=0.6, patch_artist=True,
                    boxprops=dict(linewidth=0.5),
                    whiskerprops=dict(color='black', linewidth=0.5),
                    capprops=dict(color='black', linewidth=0.5),
                    flierprops=dict(marker='o', markerfacecolor='white', markersize=4,
                                   markeredgecolor='black', markeredgewidth=0.5))

    bp['boxes'][0].set_facecolor(COLORS['Responsed'])
    bp['boxes'][1].set_facecolor(COLORS['No-response'])
    bp['medians'][0].set_color(MEDIAN_COLORS['Responsed'])
    bp['medians'][0].set_linewidth(0.8)
    bp['medians'][1].set_color(MEDIAN_COLORS['No-response'])
    bp['medians'][1].set_linewidth(0.8)

    y_max = max(np.max(responder_vals), np.max(non_responder_vals))
    y_bracket = y_max * 1.15
    ax.plot([1, 1, 2, 2], [y_bracket, y_bracket*1.05, y_bracket*1.05, y_bracket], 'k-', linewidth=0.5)
    pval_text = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
    ax.text(1.5, y_bracket*1.08, pval_text, ha='center', va='bottom', fontsize=6 * SCALE)

    ax.set_title('CEACAM5', fontsize=7 * SCALE, fontweight='normal')
    ax.set_ylabel('Expression', fontsize=6 * SCALE)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['R', 'NR'], fontsize=6 * SCALE)
    ax.set_ylim(0, y_max * 1.40)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.tick_params(axis='both', labelsize=6 * SCALE, width=0.5, length=4)
    fig.subplots_adjust(left=0.28, right=0.95, top=0.88, bottom=0.15)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    plt.savefig(os.path.join(OUTPUT_DIR, "ceacam5_pre_boxplot.png"), dpi=DPI, facecolor='white')
    plt.savefig(os.path.join(OUTPUT_DIR, "ceacam5_pre_boxplot.svg"), dpi=DPI, facecolor='white')
    plt.savefig(os.path.join(OUTPUT_DIR, "ceacam5_pre_boxplot.pdf"), dpi=DPI, facecolor='white')
    print(f"Saved to {OUTPUT_DIR}")
    plt.close()

if __name__ == '__main__':
    main()
