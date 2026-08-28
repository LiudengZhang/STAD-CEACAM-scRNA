#!/usr/bin/env python3
"""
Regenerate SVG output for panels 01_B and 01_C.
Loads the full dataset ONCE to avoid double-loading 12.5 GB.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from paths import FULL_DATASET_H5AD

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR.parent.parent / '04_Final_Panels/99_Code/00_Shared'))
try:
    from figure_config import CELL_TYPE_COLORS
except ImportError:
    CELL_TYPE_COLORS = {
        'B cells': '#E41A1C', 'Plasma cells': '#377EB8', 'T/NK cells': '#4DAF4A',
        'Dendritic cells': '#984EA3', 'MoMac': '#FF7F00', 'Mast cells': '#FFFF33',
        'Neutrophils': '#A65628', 'Fibroblasts': '#F781BF', 'Pericytes': '#999999',
        'Endothelial cells': '#66C2A5', 'Epithelial': '#FC8D62', 'Hepatocytes': '#8DA0CB'
    }

sc.set_figure_params(dpi=150, frameon=False, figsize=(12, 10), facecolor='white')
plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42})

print("=" * 60)
print("Loading full dataset ONCE for panels 01_B + 01_C...")
print("=" * 60)
adata = sc.read_h5ad(FULL_DATASET_H5AD)
print(f"Loaded: {adata.shape[0]:,} cells")

# Filter to stomach
adata_stomach = adata[adata.obs['Sample site'] == 'Stomach'].copy()
del adata  # Free memory
print(f"Stomach: {adata_stomach.shape[0]:,} cells, {adata_stomach.obs['sample'].nunique()} samples")

# Merge cell types (shared between B and C)
adata_stomach.obs['major_cell_type_merged'] = adata_stomach.obs['major_cell_type'].replace({
    'CD4+ T cells': 'T/NK cells',
    'CD8+ T cells': 'T/NK cells',
    'NK cells': 'T/NK cells',
    'DC cells': 'Dendritic cells',
    'Fibroblast': 'Fibroblasts',
    'Pericyte': 'Pericytes',
    'Hepatocyte': 'Hepatocytes'
})

# ============================================================
# Panel 01_B: UMAP by major cell type
# ============================================================
print("\n--- Panel 01_B: UMAP ---")
OUT_B = SCRIPT_DIR / '01_B'

cell_counts = adata_stomach.obs['major_cell_type_merged'].value_counts().sort_values(ascending=False)
cell_type_order = list(cell_counts.index)
colors = [CELL_TYPE_COLORS.get(ct, '#999999') for ct in cell_type_order]
adata_stomach.uns['major_cell_type_merged_colors'] = colors

fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.umap(adata_stomach, color='major_cell_type_merged', ax=ax, show=False,
           legend_loc='right margin', title='', frameon=False, size=2, alpha=0.6)
plt.tight_layout()

plt.savefig(OUT_B / 'umap_major_cell_type_stomach.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(OUT_B / 'umap_major_cell_type_stomach.svg', format='svg', bbox_inches='tight', facecolor='white')
plt.savefig(OUT_B / 'umap_major_cell_type_stomach.pdf', format='pdf', bbox_inches='tight')
plt.close()
print(f"  Saved 01_B: PNG + SVG + PDF")

# ============================================================
# Panel 01_C: Stacked bar chart
# ============================================================
print("\n--- Panel 01_C: Stacked bar ---")
OUT_C = SCRIPT_DIR / '01_C'

proportions = adata_stomach.obs.groupby(
    ['sample', 'major_cell_type_merged'], observed=True
).size().unstack(fill_value=0)
proportions = proportions.div(proportions.sum(axis=1), axis=0)

# Sort samples
sample_meta = adata_stomach.obs.groupby('sample', observed=True).agg({
    'Treatment phase': 'first',
    'stomach_pre_grouping': 'first',
    'stomach_post_grouping': 'first',
    'Patient ID': 'first'
}).reset_index()

def get_sort_key(row):
    phase = row['Treatment phase']
    if phase == 'Pre':
        response = row['stomach_pre_grouping']
    else:
        response = row['stomach_post_grouping']
    if phase == 'Pre' and response == 'Responsed':
        return (0, row['Patient ID'])
    elif phase == 'Pre' and response == 'No-response':
        return (1, row['Patient ID'])
    elif phase == 'Post' and response == 'Responsed':
        return (2, row['Patient ID'])
    else:
        return (3, row['Patient ID'])

sample_meta['sort_key'] = sample_meta.apply(get_sort_key, axis=1)
sample_meta = sample_meta.sort_values('sort_key')
sample_order = sample_meta['sample'].tolist()
proportions = proportions.reindex(sample_order)

ct_order = proportions.mean().sort_values(ascending=False).index.tolist()
proportions = proportions[ct_order]

fig, ax = plt.subplots(figsize=(10, 4))
bottom = np.zeros(len(proportions))
for cell_type in ct_order:
    color = CELL_TYPE_COLORS.get(cell_type, '#999999')
    ax.bar(range(len(proportions)), proportions[cell_type], bottom=bottom,
           color=color, label=cell_type, width=0.8)
    bottom += proportions[cell_type].values

ax.set_xticks(range(len(proportions)))
ax.set_xticklabels(proportions.index, fontsize=5, rotation=90, ha='center')
ax.set_ylabel('Proportion', fontsize=7)
ax.set_ylim(0, 1)

pre_r_end = len(sample_meta[sample_meta['sort_key'].apply(lambda x: x[0] == 0)])
pre_nr_end = pre_r_end + len(sample_meta[sample_meta['sort_key'].apply(lambda x: x[0] == 1)])
post_r_end = pre_nr_end + len(sample_meta[sample_meta['sort_key'].apply(lambda x: x[0] == 2)])

for x in [pre_r_end - 0.5, pre_nr_end - 0.5, post_r_end - 0.5]:
    if 0 < x < len(proportions):
        ax.axvline(x=x, color='black', linewidth=1.5, linestyle='-')

ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=5, frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

plt.savefig(OUT_C / '01_C_sample_balance_stomach.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(OUT_C / '01_C_sample_balance_stomach.svg', format='svg', bbox_inches='tight', facecolor='white')
plt.savefig(OUT_C / '01_C_sample_balance_stomach.pdf', format='pdf', bbox_inches='tight')
plt.close()
print(f"  Saved 01_C: PNG + SVG + PDF")

print("\n" + "=" * 60)
print("Done! Both panels regenerated with SVG output.")
print("=" * 60)
