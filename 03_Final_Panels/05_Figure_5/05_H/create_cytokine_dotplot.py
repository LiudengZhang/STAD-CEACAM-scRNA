#!/usr/bin/env python3
"""
Panel H: Dotplot of TNF, IL6, IL1B, IL1A across 13 major cell types.
Style matched to Figure 1D (major_celltype_markers dotplot).
Cell types match those in Panel F radar (NF-κB NES).
"""

import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']})
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FULL_DATASET_H5AD

OUTPUT_DIR = Path(__file__).resolve().parent

# Genes of interest
GENES = ['TNF', 'IL6', 'IL1B', 'IL1A']

# Cell type order (matching 05_F radar, using original adata labels)
CELL_TYPE_ORDER = [
    'B cells',
    'DC cells',
    'Endothelial cells',
    'Epithelial cells',
    'Fibroblast',
    'Mast cells',
    'Monocytes/Macrophages',
    'Neutrophils',
    'NK cells',
    'Pericyte',
    'Plasma cells',
    'CD4+ T cells',
    'CD8+ T cells',
]


def main():
    print("=" * 60)
    print("Panel H: Cytokine Dotplot (TNF, IL6, IL1B, IL1A)")
    print("=" * 60)

    print("\n[1/3] Loading data...")
    adata = sc.read_h5ad(FULL_DATASET_H5AD)
    print(f"  {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # Filter to 13 cell types (exclude Hepatocyte)
    mask = adata.obs['major_cell_type'].isin(CELL_TYPE_ORDER)
    adata = adata[mask].copy()
    print(f"  After filtering: {adata.shape[0]:,} cells, {adata.obs['major_cell_type'].nunique()} types")

    # Validate genes
    var_names = adata.raw.var_names if adata.raw is not None else adata.var_names
    available = [g for g in GENES if g in var_names]
    missing = [g for g in GENES if g not in var_names]
    print(f"  Genes available: {available}")
    if missing:
        print(f"  WARNING — missing: {missing}")

    print("\n[2/3] Creating dotplot...")
    fig_width = max(6, len(available) * 1.2)
    fig_height = max(6, len(CELL_TYPE_ORDER) * 0.5)
    plt.figure(figsize=(fig_width, fig_height))

    sc.pl.dotplot(
        adata,
        var_names=available,
        groupby='major_cell_type',
        categories_order=CELL_TYPE_ORDER,
        use_raw=(adata.raw is not None),
        cmap='Reds',
        show=False,
        save=None,
        standard_scale='var',
    )

    plt.tight_layout()

    print("\n[3/3] Saving...")
    out = OUTPUT_DIR / 'cytokine_dotplot.png'
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(out.with_suffix('.svg'), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {out}")
    print("Done!")


if __name__ == "__main__":
    main()
