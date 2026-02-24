#!/usr/bin/env python3
"""
Create Cell State Correlation Matrix Heatmap (Figure 04 Panel A)

Generates a correlation heatmap showing cell-cell relationships organized by modules.
"""

import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']})

# Parameters
FIGURE_SIZE = (12, 12)
DPI = 300
COLORMAP = 'RdBu_r'
VMIN, VMAX = -1, 1
CENTER = 0
MASK_DIAGONAL = True
BOUNDARY_COLOR = 'black'
BOUNDARY_WIDTH = 2
CELL_LABEL_SIZE = 8
MODULE_LABEL_SIZE = 12
COLORBAR_SHRINK = 0.6

# Short display names for y-axis labels
SHORT_NAMES = {
    # Module 1 — T/NK/DC
    'C2_CD8_MAIT_KLRB1':              'C2_CD8_MAIT',
    'C9_CD4_Th17_IL17A':              'C9_CD4_Th17',
    'C0_CD4_Treg_FOXP3':              'C0_CD4_Treg',
    'C8_CD4_Temra_KLRG1':             'C8_CD4_Temra',
    'C2_CD4_Tcm_CCR7':                'C2_CD4_Tcm',
    'C3_CD8_Tcm_CCR7':                'C3_CD8_Tcm',
    'C0_CD8_Cytotoxic_CCL':           'C0_CD8_CCL',
    'C2_NKT_CD3D':                    'C2_NKT',
    'C3_Prolif_MKI67':                'C3_NK_Prolif',
    'C6_CD8_Tex_PDCD1':               'C6_CD8_Tex',
    'C5_CD4_Prolif_MKI67':            'C5_CD4_Prolif',
    'C5_CD8_Prolif_MKI67':            'C5_CD8_Prolif',
    'C1_DC_Conventional_S100A4_high':  'C1_DC_Conv',
    'C3_DC_CCR7+':                    'C3_DC_CCR7',
    'C7_CD8_ISG_ISG15':               'C7_CD8_ISG',
    'C4_CD8_Temra_KLRG1':             'C4_CD8_Temra',
    # Module 2 — MoMac
    'C4_Mono_Alternative_CD16':       'C4_Mono_CD16',
    'C2_MoMac_Intermediate_HLA-DRA':  'C2_MoMac_Inter',
    'C1_Mono_Classic_CD14':           'C1_Mono_CD14',
    'C3_Mac_Inflam_IL1B':             'C3_Mac_IL1B',
    'C0_Mac_Classic_TREM2':           'C0_Mac_TREM2',
    'C6_Mac_Metallothionein_MT1G':    'C6_Mac_MT1G',
    # Module 3 — Mixed immune
    'C4_CD4_Tfh_PDCD1':              'C4_CD4_Tfh',
    'C1_CD4_Tem_IL7R':               'C1_CD4_Tem',
    'C1_CD8_Cytotoxic_DUSP1':        'C1_CD8_DUSP1',
    'C2_DC_Activated_IL1B_high':      'C2_DC_Act',
    'C5_DC_WDFY4+':                  'C5_DC_WDFY4',
    'C6_DC_C1Q+_STMN1+':            'C6_DC_C1Q_STMN1',
    'C5_Mac_Prolif_MKI67':           'C5_Mac_Prolif',
    'C4_DC_C1Q+_STMN1-':            'C4_DC_C1Q',
    'C7_DC_LTB+':                    'C7_DC_LTB',
    # Module 5 — B/Plasma
    'C6_B_Germinal_Center':          'C6_B_GC',
    'C2_B_Memory_IGHA1':             'C2_B_Mem_IGHA1',
    'C4_B_Naive_TCL1A':              'C4_B_Naive_TCL1A',
    'C3_B_Naive_NR4A1':              'C3_B_Naive_NR4A1',
    'C1_B_Memory_NR4A1':             'C1_B_Mem_NR4A1',
}

def main():
    """Generate correlation heatmap."""
    print("=" * 60)
    print("Cell State Correlation Matrix - Figure 04A")
    print("=" * 60)

    # Get script directory
    script_dir = Path(__file__).parent

    # Load data
    print("\nLoading data...")
    correlation_matrix = pd.read_csv(script_dir / 'correlation_matrix.csv', index_col=0)
    module_mappings = pd.read_csv(script_dir / 'module_mappings.csv')

    print(f"  Correlation matrix: {correlation_matrix.shape}")
    print(f"  Module mappings: {module_mappings.shape}")

    # Order by modules
    print("\nOrdering by modules...")
    module_mappings_sorted = module_mappings.sort_values('module')
    ordered_states = module_mappings_sorted['minor_cell_state'].tolist()
    correlation_ordered = correlation_matrix.loc[ordered_states, ordered_states]

    # Apply short display names and strip Cx_ prefix
    strip_cx = lambda s: re.sub(r'^C\d+_', '', s)
    short_index = [strip_cx(SHORT_NAMES.get(s, s)) for s in correlation_ordered.index]
    correlation_ordered.index = short_index
    correlation_ordered.columns = short_index

    # Calculate module boundaries
    module_counts = module_mappings_sorted.groupby('module').size()
    module_boundaries = np.cumsum(module_counts.values)[:-1]
    print(f"  Module counts: {dict(module_counts)}")

    # Create figure
    print("\nCreating heatmap...")
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Create mask for diagonal
    mask = np.eye(len(correlation_ordered), dtype=bool) if MASK_DIAGONAL else None

    # Create heatmap
    sns.heatmap(
        correlation_ordered,
        cmap=COLORMAP,
        vmin=VMIN,
        vmax=VMAX,
        center=CENTER,
        mask=mask,
        square=True,
        cbar_kws={'label': 'Spearman Correlation', 'shrink': COLORBAR_SHRINK},
        ax=ax,
        xticklabels=False,
        yticklabels=correlation_ordered.index,
        linewidths=0.5,
        linecolor='lightgray'
    )
    # De-rasterize heatmap (seaborn uses pcolormesh which defaults to rasterized in SVG)
    for coll in ax.collections:
        coll.set_rasterized(False)

    # Adjust cell labels
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=CELL_LABEL_SIZE)

    # Add module boundaries
    for boundary in module_boundaries:
        ax.axhline(y=boundary, color=BOUNDARY_COLOR, linewidth=BOUNDARY_WIDTH)
        ax.axvline(x=boundary, color=BOUNDARY_COLOR, linewidth=BOUNDARY_WIDTH)

    # Add module labels at bottom
    module_positions = []
    cumsum = 0
    for count in module_counts:
        module_positions.append(cumsum + count / 2)
        cumsum += count

    ax.set_xticks(module_positions)
    ax.set_xticklabels(
        ['IM-T/NK/DC', 'IM-MoMac', 'IM-Mixed', 'IM-Neutrophil', 'IM-B/Plasma'],
        fontsize=MODULE_LABEL_SIZE
    )
    ax.xaxis.tick_bottom()

    plt.tight_layout()

    # Save
    output_path = script_dir / 'correlation_heatmap_k5.png'
    fig.savefig(output_path, dpi=DPI, bbox_inches='tight')
    fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')

    file_size = output_path.stat().st_size / 1024
    print(f"\n✓ Saved: {output_path.name}")
    print(f"  Size: {file_size:.1f} KB")
    print("=" * 60)

    plt.close()

if __name__ == "__main__":
    main()
