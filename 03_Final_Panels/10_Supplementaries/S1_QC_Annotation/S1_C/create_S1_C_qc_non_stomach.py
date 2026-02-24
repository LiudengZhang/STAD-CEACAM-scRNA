#!/usr/bin/env python3
"""
S1_C: QC Metrics — Non-Stomach Samples (38)
Boxplots of QC metrics per sample, colored by tissue site.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import FULL_DATASET_H5AD

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# =============================================================================
# Configuration — 4× scaling
# =============================================================================
SCALE = 4
FONT_SIZE = 6 * SCALE
LINEWIDTH = 1 * SCALE
FIG_WIDTH_CM = 24 * SCALE
FIG_HEIGHT_CM = 5 * SCALE
DPI = 300

FIG_WIDTH = FIG_WIDTH_CM / 2.54
FIG_HEIGHT = FIG_HEIGHT_CM / 2.54

# Tissue colors
TISSUE_COLORS = {
    'Liver': '#e6ab02',
    'Peripheral blood': '#66a61e',
    'Lymph node': '#7570b3',
    'Ovary': '#e7298a',
    'Metastatic lymph node': '#a6761d',
}

OUTPUT_DIR = Path(__file__).parent
OUTPUT_FILE = OUTPUT_DIR / "panel_S1_C.png"

# =============================================================================
# Load data — non-stomach only
# =============================================================================
print("Loading data...")
adata = sc.read_h5ad(FULL_DATASET_H5AD, backed='r')

mask = adata.obs['Sample site'] != 'Stomach'
plot_df = adata.obs.loc[mask, ['sample', 'Sample site',
                                'n_genes_by_counts', 'total_counts', 'pct_counts_mt']].copy()

print(f"Non-stomach cells: {len(plot_df)}")

# =============================================================================
# Sample ordering — group by tissue
# =============================================================================
sample_info = plot_df[['sample', 'Sample site']].drop_duplicates('sample')

tissue_order = ['Liver', 'Peripheral blood', 'Lymph node', 'Ovary', 'Metastatic lymph node']
sample_info['sort_key'] = sample_info['Sample site'].map({t: i for i, t in enumerate(tissue_order)})
sample_info = sample_info.sort_values(['sort_key', 'sample'])
sample_order = sample_info['sample'].tolist()

sample_to_tissue = dict(zip(sample_info['sample'], sample_info['Sample site']))

for t in tissue_order:
    n = sum(sample_info['Sample site'] == t)
    if n > 0:
        print(f"  {t}: {n} samples")

# Find tissue boundaries for vertical separators
boundaries = []
prev_tissue = None
for i, s in enumerate(sample_order):
    t = sample_to_tissue[s]
    if prev_tissue is not None and t != prev_tissue:
        boundaries.append(i - 0.5)
    prev_tissue = t

# =============================================================================
# Create figure
# =============================================================================
plt.rcParams.update({
    'font.size': FONT_SIZE,
    'axes.linewidth': LINEWIDTH,
    'axes.labelsize': FONT_SIZE,
    'axes.titlesize': FONT_SIZE,
    'xtick.labelsize': FONT_SIZE * 0.5,
    'ytick.labelsize': FONT_SIZE * 0.7,
    'xtick.major.width': LINEWIDTH * 0.5,
    'ytick.major.width': LINEWIDTH * 0.5,
    'xtick.major.size': 2 * SCALE,
    'ytick.major.size': 2 * SCALE,
    'legend.fontsize': FONT_SIZE * 0.6,
    'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

fig, axes = plt.subplots(1, 3, figsize=(FIG_WIDTH, FIG_HEIGHT))
plt.subplots_adjust(top=0.85, bottom=0.30, left=0.08, right=0.82, wspace=0.55)

titles = ['Genes per Cell', 'Total Counts', 'Mitochondrial %']
metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
ylabels = ['Number of genes', 'Total UMI counts', '% mitochondrial']

for ax, metric, title, ylabel in zip(axes, metrics, titles, ylabels):
    colors = [TISSUE_COLORS[sample_to_tissue[s]] for s in sample_order]

    box = ax.boxplot(
        [plot_df.loc[plot_df['sample'] == s, metric].values for s in sample_order],
        positions=list(range(len(sample_order))),
        widths=0.7,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color='black', linewidth=LINEWIDTH * 0.5),
        whiskerprops=dict(color='black', linewidth=LINEWIDTH * 0.3),
        capprops=dict(color='black', linewidth=LINEWIDTH * 0.3),
        boxprops=dict(linewidth=LINEWIDTH * 0.3),
    )

    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        patch.set_edgecolor('black')

    for b in boundaries:
        ax.axvline(x=b, color='gray', linestyle=':', linewidth=LINEWIDTH * 0.5, alpha=0.7)

    if metric == 'pct_counts_mt':
        ax.axhline(y=25, color='#666666', linestyle='--', linewidth=LINEWIDTH * 0.7, zorder=10)

    ax.set_title(title, fontsize=FONT_SIZE * 0.9, pad=5)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Sample')
    ax.set_xticks(list(range(len(sample_order))))
    ax.set_xticklabels(sample_order, rotation=90, ha='center')
    ax.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=LINEWIDTH * 0.3)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-0.5, len(sample_order) - 0.5)

# Legend
legend_handles = [
    mpatches.Patch(facecolor=TISSUE_COLORS[t], edgecolor='black', alpha=0.8, label=t)
    for t in tissue_order if t in sample_to_tissue.values()
]
fig.legend(handles=legend_handles, loc='upper right', ncol=1,
           bbox_to_anchor=(1.0, 0.85), frameon=True, fancybox=False, edgecolor='gray')

plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_FILE.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

print(f"Saved: {OUTPUT_FILE}")
print("Done!")
