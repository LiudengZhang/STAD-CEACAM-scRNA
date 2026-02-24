#!/usr/bin/env python3
"""
S2_A: Epithelial marker dotplot by sub-cluster.
Top marker genes per epithelial sub-cluster.
Simple sc.pl.dotplot style matching S1_F.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import EPITHELIAL_TUMOR_SCORED_H5AD

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib

# =============================================================================
# Configuration
# =============================================================================
SCALE = 4
DPI = 300
OUTPUT_DIR = Path(__file__).parent
OUTPUT_PNG = OUTPUT_DIR / "panel_S2_A.png"

# Markers per sub-cluster (from cluster naming + known biology)
markers = {
    'PTMA': ['PTMA', 'STMN1'],
    'KRT19': ['KRT19', 'KRT8'],
    'CEACAM5/6': ['CEACAM6', 'CEACAM5'],
    'Chief-like': ['PGC', 'LIPF'],
    'MUC5AC': ['MUC5AC', 'TFF1'],
    'Stem TPX2': ['TPX2', 'MKI67'],
    'Stem SPINK4': ['SPINK4', 'SPINK1'],
    'MT1E': ['MT1E', 'MT2A'],
    'CD74': ['CD74', 'HLA-DRA'],
}

label_map = {
    'C0_Epi_PTMA': 'PTMA',
    'C1_Epi_KRT19': 'KRT19',
    'C2_Epi_CEACAM6': 'CEACAM5/6',
    'C3_Epi_Chief_Like_PGC': 'Chief-like',
    'C4_Epi_MUC5AC': 'MUC5AC',
    'C5_Epi_Stem_Like_TPX2': 'Stem TPX2',
    'C6_Epi_Stem_Like_SPINK4': 'Stem SPINK4',
    'C7_Epi_MT1E': 'MT1E',
    'C8_Epi_CD74': 'CD74',
}

# =============================================================================
# Load data
# =============================================================================
print("Loading Epithelial data...")
adata = sc.read_h5ad(EPITHELIAL_TUMOR_SCORED_H5AD)
adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
print(f"Stomach epithelial: {adata.n_obs} cells")

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
# Create dotplot — simple style matching S1_F
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
