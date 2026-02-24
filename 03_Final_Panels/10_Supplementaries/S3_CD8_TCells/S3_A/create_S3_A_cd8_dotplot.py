#!/usr/bin/env python3
"""
S3_A: CD8+ T cell marker dotplot by sub-cluster.
Top marker genes per CD8 sub-cluster.
Simple sc.pl.dotplot style matching S2_A / S1_F.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import TCD8_H5AD

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib

# =============================================================================
# Configuration
# =============================================================================
SCALE = 4
DPI = 300
OUTPUT_DIR = Path(__file__).parent
OUTPUT_PNG = OUTPUT_DIR / "panel_S3_A.png"

markers = {
    'Cytotoxic CCL': ['CCL4', 'CCL5'],
    'Cytotoxic DUSP1': ['DUSP1', 'JUN'],
    'MAIT': ['KLRB1', 'SLC4A10'],
    'Tcm': ['CCR7', 'SELL'],
    'TEMRA': ['KLRG1', 'CX3CR1'],
    'Proliferating': ['MKI67', 'TOP2A'],
    'Tex': ['PDCD1', 'HAVCR2', 'LAG3'],
    'ISG': ['ISG15', 'MX1'],
}

label_map = {
    'C0_CD8_Cytotoxic_CCL': 'Cytotoxic CCL',
    'C1_CD8_Cytotoxic_DUSP1': 'Cytotoxic DUSP1',
    'C2_CD8_MAIT_KLRB1': 'MAIT',
    'C3_CD8_Tcm_CCR7': 'Tcm',
    'C4_CD8_Temra_KLRG1': 'TEMRA',
    'C5_CD8_Prolif_MKI67': 'Proliferating',
    'C6_CD8_Tex_PDCD1': 'Tex',
    'C7_CD8_ISG_ISG15': 'ISG',
}

# =============================================================================
# Load data
# =============================================================================
print("Loading CD8 data...")
adata = sc.read_h5ad(TCD8_H5AD)
adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
print(f"Stomach CD8: {adata.n_obs} cells")

# Clean labels
states = adata.obs['minor_cell_state'].astype(str)
adata.obs['cluster_label'] = states.map(label_map).fillna(states)
cluster_order = list(label_map.values())

# Verify markers
all_genes = [g for genes in markers.values() for g in genes]
available = [g for g in all_genes if g in adata.var_names]
missing = [g for g in all_genes if g not in adata.var_names]
if missing:
    print(f"Missing genes: {missing}")
    markers = {k: [g for g in v if g in adata.var_names] for k, v in markers.items()}
    markers = {k: v for k, v in markers.items() if v}
print(f"Using {len(available)}/{len(all_genes)} markers")

# =============================================================================
# Create dotplot — simple style matching S2_A / S1_F
# =============================================================================
plt.rcParams.update({
    'font.size': 5 * SCALE,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

sc.pl.dotplot(
    adata,
    var_names=available,
    groupby='cluster_label',
    categories_order=cluster_order,
    use_raw=(adata.raw is not None),
    cmap='Reds',
    show=False,
    save=None,
    standard_scale='var',
)

fig = plt.gcf()
fig.set_size_inches(14, 5)

print(f"Saving to {OUTPUT_PNG}")
fig.savefig(OUTPUT_PNG, dpi=DPI, bbox_inches='tight', facecolor='white')
fig.savefig(OUTPUT_PNG.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close('all')
print("Done!")
