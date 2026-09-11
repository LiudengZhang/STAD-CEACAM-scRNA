#!/usr/bin/env python3
"""
Generate Cell Type Composition Stacked Bar - Stomach Samples Only (32 samples)
==============================================================================
Creates a stacked bar chart showing cell type proportions for stomach samples only.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FULL_DATASET_H5AD, MANUSCRIPT

# Sample ID mapping from Supplementary Table 1
ST1_CSV = (MANUSCRIPT / '04_Tables'
           / 'ST1_patient_sample_characteristics.csv')

# Set2 + Set3 palette (12 major cell types, consistent across Fig 1B & 1C)
CELL_TYPE_COLORS = {
    'Epithelial cells': '#e5c494',
    'T/NK cells': '#66c2a5',
    'Monocytes/Macrophages': '#fc8d62',
    'Plasma cells': '#8da0cb',
    'B cells': '#ffd92f',
    'Endothelial cells': '#a6d854',
    'Fibroblasts': '#e78ac3',
    'Neutrophils': '#b3b3b3',
    'Mast cells': '#fb8072',
    'Pericytes': '#bebada',
    'Dendritic cells': '#80b1d3',
    'Hepatocytes': '#bc80bd',
}

plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']})

# Paths
# The upstream integration pipeline produced this object; the deposited
# copy is FULL_DATASET_H5AD, which is what the line below reads.
DATA_FILE = FULL_DATASET_H5AD
OUTPUT_DIR = Path(__file__).parent

print("=" * 80)
print("GENERATING CELL TYPE STACKED BAR - STOMACH SAMPLES ONLY")
print("=" * 80)

# Load data
print(f"\nLoading data from: {DATA_FILE}")
adata = sc.read_h5ad(DATA_FILE)
print(f"Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Filter to stomach samples only
adata_stomach = adata[adata.obs['Sample site'] == 'Stomach'].copy()
print(f"Stomach samples: {adata_stomach.shape[0]:,} cells")
print(f"Number of samples: {adata_stomach.obs['sample'].nunique()}")

# Merge T/NK cells and standardize names
adata_stomach.obs['major_cell_type_merged'] = adata_stomach.obs['major_cell_type'].replace({
    'CD4+ T cells': 'T/NK cells',
    'CD8+ T cells': 'T/NK cells',
    'NK cells': 'T/NK cells',
    'DC cells': 'Dendritic cells',
    'Fibroblast': 'Fibroblasts',
    'Pericyte': 'Pericytes',
    'Hepatocyte': 'Hepatocytes'
})

# Calculate proportions per sample
proportions = adata_stomach.obs.groupby(['sample', 'major_cell_type_merged'], observed=True).size().unstack(fill_value=0)
proportions = proportions.div(proportions.sum(axis=1), axis=0)

# Specimen identifier -> the study sample ID printed in the paper.
# Supplementary Table 1 has no 'Original Sample ID' column, so the
# crosswalk comes from the dataset object, which carries both labels;
# ST1 decides which candidate is the published one. An unmapped
# specimen is an error, not a raw identifier on the axis.
st1 = pd.read_csv(ST1_CSV)
st1_labels = set(st1[st1['Anatomical site'] == 'Stomach']['Sample'].astype(str))
pairs = (adata_stomach.obs[['sample', 'Sample ID']]
         .drop_duplicates().astype(str))
orig_to_sample = dict(pairs[pairs['Sample ID'].isin(st1_labels)].values)
unmapped = sorted(set(pairs['sample']) - set(orig_to_sample))
if unmapped:
    raise SystemExit(
        'no Supplementary Table 1 label for: ' + ', '.join(unmapped))

# Sort samples alphabetically (no group splitting)
sample_order = sorted(proportions.index.tolist())
proportions = proportions.reindex(sample_order)

# Map to ST1 Sample IDs for x-axis labels
x_labels = [orig_to_sample[s] for s in sample_order]

# Define cell type order (most abundant first)
cell_type_order = proportions.mean().sort_values(ascending=False).index.tolist()
proportions = proportions[cell_type_order]

# Create figure
fig, ax = plt.subplots(figsize=(10, 4))

# Plot stacked bars — continuous block, no gaps
positions = np.arange(len(proportions))
bottom = np.zeros(len(proportions))
for cell_type in cell_type_order:
    color = CELL_TYPE_COLORS.get(cell_type, '#999999')
    ax.bar(positions, proportions[cell_type], bottom=bottom,
           color=color, label=cell_type, width=0.8, edgecolor='white', linewidth=0.3)
    bottom += proportions[cell_type].values

# Styling — patient sample IDs from ST1
ax.set_xticks(positions)
ax.set_xticklabels(x_labels, rotation=90, fontsize=12, ha='center')
ax.set_ylabel('Proportion', fontsize=7)
ax.set_ylim(0, 1)

# Legend
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=5, frameon=False)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

# Save
output_file = OUTPUT_DIR / '01_C_sample_balance_stomach.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_DIR / '01_C_sample_balance_stomach.svg', format='svg', bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_DIR / '01_C_sample_balance_stomach.pdf', format='pdf', bbox_inches='tight')
plt.close()

print(f"\n✓ Saved: {output_file}")
print(f"  Samples: {len(proportions)}")
print("=" * 80)
