#!/usr/bin/env python3
"""
Panel E: MoMac Subset Marker Dotplot
=====================================

Generates a dotplot showing canonical markers for all 7 MoMac minor cell states.
3 markers per state, 21 genes total.

Adapted from 01_Figure_1/01_D/create_major_celltype_dotplot.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import OrderedDict
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD

# Output directory
OUTPUT_DIR = Path(__file__).resolve().parent

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=300, facecolor='white')
plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']})

# ============================================================================
# CANONICAL MARKERS BY MOMAC SUBSET
# ============================================================================

# Short display names for dotplot y-axis
SHORT_NAMES = OrderedDict([
    ('C0_Mac_Classic_TREM2',           'Mac_TREM2'),
    ('C1_Mono_Classic_CD14',           'Mono_CD14'),
    ('C2_MoMac_Intermediate_HLA-DRA',  'MoMac_Inter'),
    ('C3_Mac_Inflam_IL1B',             'Mac_IL1B'),
    ('C4_Mono_Alternative_CD16',       'Mono_CD16'),
    ('C5_Mac_Prolif_MKI67',            'Mac_Prolif'),
    ('C6_Mac_Metallothionein_MT1G',    'Mac_MT1G'),
])

CANONICAL_MARKERS = OrderedDict([
    ('Mac_TREM2',    ['TREM2', 'C1QA', 'APOE']),
    ('Mono_CD14',    ['CD14', 'S100A8', 'VCAN']),
    ('MoMac_Inter',  ['HLA-DRA', 'CST3', 'CLEC10A']),
    ('Mac_IL1B',     ['IL1B', 'TNF', 'CXCL8']),
    ('Mono_CD16',    ['CD16', 'CDKN1C', 'LST1']),
    ('Mac_Prolif',   ['MKI67', 'TOP2A', 'STMN1']),
    ('Mac_MT1G',     ['MT1G', 'MT2A', 'MT1X']),
])

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    print("=" * 80)
    print("MoMac Subset Marker Dotplot Generator")
    print("=" * 80)

    # 1. Load Data
    print("\n[1/4] Loading data...")
    print(f"  Path: {MOMAC_H5AD}")

    adata = sc.read_h5ad(MOMAC_H5AD)
    print(f"  Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    if 'minor_cell_state' not in adata.obs.columns:
        raise ValueError("minor_cell_state column not found in adata.obs")

    # Rename to short display names
    adata.obs['momac_short'] = adata.obs['minor_cell_state'].replace(SHORT_NAMES)

    # Rename genes for display (FCGR3A → CD16)
    GENE_DISPLAY = {'FCGR3A': 'CD16'}
    adata.var_names = pd.Index([GENE_DISPLAY.get(g, g) for g in adata.var_names])
    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        raw_adata.var_names = pd.Index([GENE_DISPLAY.get(g, g) for g in raw_adata.var_names])
        adata.raw = raw_adata

    # 2. Validate Markers
    print("\n[2/4] Validating marker genes...")
    all_markers = []
    for cell_type, markers in CANONICAL_MARKERS.items():
        all_markers.extend(markers)

    if adata.raw is not None:
        var_names = adata.raw.var_names
    else:
        var_names = adata.var_names

    available_markers = [m for m in all_markers if m in var_names]
    missing_markers = [m for m in all_markers if m not in var_names]

    print(f"  Available: {len(available_markers)}/{len(all_markers)} markers")
    if missing_markers:
        print(f"  Missing markers: {missing_markers}")

    filtered_markers = OrderedDict()
    for cell_type, markers in CANONICAL_MARKERS.items():
        avail = [m for m in markers if m in var_names]
        if avail:
            filtered_markers[cell_type] = avail

    # 3. Prepare
    print("\n[3/4] Preparing visualization...")
    cell_type_order = list(CANONICAL_MARKERS.keys())

    plot_genes = []
    for cell_type in cell_type_order:
        if cell_type in filtered_markers:
            plot_genes.extend(filtered_markers[cell_type])

    # 4. Create Dotplot
    print("\n[4/4] Creating dotplot...")
    fig_width = max(10, len(plot_genes) * 0.45)
    fig_height = max(5, len(cell_type_order) * 0.6)

    plt.figure(figsize=(fig_width, fig_height))

    sc.pl.dotplot(
        adata,
        var_names=plot_genes,
        groupby='momac_short',
        categories_order=cell_type_order,
        use_raw=(adata.raw is not None),
        cmap='Reds',
        show=False,
        save=None,
        standard_scale='var',
    )

    # De-rasterize all axes (scanpy dotplot may rasterize internally)
    for ax in plt.gcf().get_axes():
        for coll in ax.collections:
            coll.set_rasterized(False)
        for img in ax.images:
            img.set_rasterized(False)

    plt.tight_layout()

    output_path = OUTPUT_DIR / 'momac_marker_dotplot.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")
    print(f"  Saved: {output_path.with_suffix('.svg')}")

    # Summary table
    summary_data = []
    for cell_type in cell_type_order:
        if cell_type in filtered_markers:
            markers = filtered_markers[cell_type]
            summary_data.append({
                'Cell_State': cell_type,
                'Markers': ', '.join(markers),
                'N_Markers': len(markers)
            })

    summary_df = pd.DataFrame(summary_data)
    summary_path = OUTPUT_DIR / 'momac_marker_summary.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f"  Saved: {summary_path}")
    print("\nComplete!")


if __name__ == "__main__":
    main()
