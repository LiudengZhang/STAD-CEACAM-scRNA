#!/usr/bin/env python3
"""
S1_H: Cell Type Composition Stacked Bar - Non-Stomach Samples (38 samples)
==========================================================================
Creates a stacked bar chart showing cell type proportions for all non-stomach samples.
Adapted from Figure 1C (stomach-only version).
"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
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

DATA_FILE = FULL_DATASET_H5AD
OUTPUT_DIR = Path(__file__).parent

print("=" * 80)
print("GENERATING CELL TYPE STACKED BAR - NON-STOMACH SAMPLES")
print("=" * 80)

# Load data
print(f"\nLoading data from: {DATA_FILE}")
adata = sc.read_h5ad(DATA_FILE)
print(f"Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Filter to non-stomach samples
adata_other = adata[adata.obs['Sample site'] != 'Stomach'].copy()
print(f"Non-stomach samples: {adata_other.shape[0]:,} cells")
print(f"Number of samples: {adata_other.obs['sample'].nunique()}")
print(f"Tissue sites: {adata_other.obs['Sample site'].value_counts().to_dict()}")

# Merge T/NK cells and standardize names
adata_other.obs['major_cell_type_merged'] = adata_other.obs['major_cell_type'].replace({
    'CD4+ T cells': 'T/NK cells',
    'CD8+ T cells': 'T/NK cells',
    'NK cells': 'T/NK cells',
    'DC cells': 'Dendritic cells',
    'Fibroblast': 'Fibroblasts',
    'Pericyte': 'Pericytes',
    'Hepatocyte': 'Hepatocytes'
})

# Calculate proportions per sample
proportions = adata_other.obs.groupby(['sample', 'major_cell_type_merged'], observed=True).size().unstack(fill_value=0)
proportions = proportions.div(proportions.sum(axis=1), axis=0)

# Specimen identifier -> the study sample ID printed in the paper.
# Supplementary Table 1 is the authority for these labels. The dataset
# carries more than one candidate label for some specimens, so the
# candidate that appears in ST1 is the one used; anything unmatched is
# an error rather than a raw identifier on the axis.
st1 = pd.read_csv(ST1_CSV)
st1_labels = set(st1[st1['Anatomical site'] != 'Stomach']['Sample'].astype(str))
pairs = (adata_other.obs[['sample', 'Sample ID']]
         .drop_duplicates().astype(str))
orig_to_sample = dict(pairs[pairs['Sample ID'].isin(st1_labels)].values)
unmapped = sorted(set(pairs['sample']) - set(orig_to_sample))
if unmapped:
    raise SystemExit(
        'no Supplementary Table 1 label for: ' + ', '.join(unmapped))

# Sort samples alphabetically
sample_order = sorted(proportions.index.tolist())
proportions = proportions.reindex(sample_order)

# Map to ST1 Sample IDs for x-axis labels
x_labels = [orig_to_sample.get(s, s) for s in sample_order]

# Define cell type order (most abundant first)
cell_type_order = proportions.mean().sort_values(ascending=False).index.tolist()
proportions = proportions[cell_type_order]

# Create figure (wider for 38 samples)
fig, ax = plt.subplots(figsize=(12, 4))

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
output_file = OUTPUT_DIR / 'panel_S1_H.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_DIR / 'panel_S1_H.svg', format='svg', bbox_inches='tight', facecolor='white')
plt.close()

print(f"\n✓ Saved: {output_file}")
print(f"  Samples: {len(proportions)}")
print("=" * 80)
