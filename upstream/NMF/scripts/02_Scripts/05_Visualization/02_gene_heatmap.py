#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""02_gene_heatmap.py - Gene-MetaProgram heatmap"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
MP_DIR = BASE_DIR / "03_Output" / "03_MetaPrograms"
OUTPUT_DIR = BASE_DIR / "03_Output" / "05_Figures"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

sig_matrix = pd.read_csv(MP_DIR / "signature_gene_matrix.csv", index_col=0)
sig_dir = MP_DIR / "metaprogram_signatures"

gene_scores = {}
for mp_col in sig_matrix.columns:
    scores_file = sig_dir / f"{mp_col}_gene_scores.csv"
    if scores_file.exists():
        scores_df = pd.read_csv(scores_file)
        gene_scores[mp_col] = dict(zip(scores_df['gene'], scores_df['frequency']))

weighted_matrix = sig_matrix.copy().astype(float)
for mp in gene_scores:
    if mp in weighted_matrix.columns:
        for gene in weighted_matrix.index:
            if weighted_matrix.loc[gene, mp] > 0:
                weighted_matrix.loc[gene, mp] = gene_scores[mp].get(gene, 1)

present_genes = sig_matrix.sum(axis=1) > 0
filtered_matrix = weighted_matrix.loc[present_genes]

top_n = 20
selected_genes = set()
for mp in gene_scores:
    mp_genes = sorted(gene_scores[mp].items(), key=lambda x: -x[1])[:top_n]
    selected_genes.update([g[0] for g in mp_genes])
selected_genes = [g for g in selected_genes if g in filtered_matrix.index]
plot_matrix = filtered_matrix.loc[selected_genes]

if len(plot_matrix) > 2:
    gene_linkage = linkage(pdist(plot_matrix.values, metric='euclidean'), method='average')
    gene_order = dendrogram(gene_linkage, no_plot=True)['leaves']
    plot_matrix = plot_matrix.iloc[gene_order]

fig, ax = plt.subplots(figsize=(10, max(8, len(plot_matrix) * 0.25)))
sns.heatmap(plot_matrix, ax=ax, cmap=sns.color_palette("YlOrRd", as_cmap=True), linewidths=0.5, linecolor='white', cbar_kws={'label': 'Gene Frequency'})
ax.set_xlabel('Meta-Program', fontsize=12)
ax.set_ylabel('Genes', fontsize=12)
ax.set_title('Gene-MetaProgram Signature Heatmap', fontsize=14, fontweight='bold')
ax.set_yticklabels(ax.get_yticklabels(), fontsize=8)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "02_gene_metaprogram_heatmap.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print("Figure 2 complete")
