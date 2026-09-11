#!/usr/bin/env python3
"""
Figure 4 panel A - Spearman correlation between every pair of minor cell
states, ordered by module, with the five modules named along the bottom.

  printed panel  Figure 4 A       (PROVENANCE.csv - the directory is "04_A";
                                   do NOT read the directory as the letter)

The correlation matrix, the module ordering, the short row names, the masked
diagonal, the colour map, the limits and the module boundaries are unchanged.
Only the canvas and the type change: the panel is drawn at the millimetre
rectangle it prints in and set in the figure's one type system.

ONE ROW NAME IS RESTORED
    The page names the C3 state MoMac_IL1B, the form the lineage analysis
    settled on; the script still carried the earlier Mac_IL1B. The row, its
    position and its correlations do not move - only the printed name.

THE ROW LABELS ARE AN EXEMPTION
    Fifty-six cell states are named down the side of a square matrix that has
    to fit a 94.3 mm slot, which leaves 1.57 mm - about 4.5 pt - of pitch per
    row. The names are kept at the size they print at, 3.32 pt: they are
    neither shrunk nor abbreviated, because every alternative loses a cell-state
    name a reader needs. Every other string on the panel is set to the figure's
    type spec. This is the one place in Figure 4 where 6 pt is not reached.

MARK
    The earlier drawing set no SCALE. It drew a 12 x 12 inch canvas and cropped
    it to the ink on save, and that crop reached the page at a fixed fraction of
    its natural size. The fraction is legible in the published panel itself:
    8 pt row labels print at 3.322 pt, 12 pt module labels at 4.984 pt, and
    matplotlib's own 10 pt colour-bar ticks at 4.153 pt - one ratio, 0.41525.
    So

        MARK = ROW_LABEL_PT / CELL_LABEL_SIZE

    and the two lengths the drawing sets in points, the cell grid and the module
    boundary, are written here as they were written there with `* MARK`
    appended.

Input : correlation_matrix.csv; module_mappings.csv (beside this script)
Output: this directory / correlation_heatmap_k5.{svg,pdf,png}
"""

import pandas as pd
import numpy as np
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
import panel_style_cns as style           # noqa: E402
import slots                              # noqa: E402

PANEL_LETTER = "A"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(4, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(4, PANEL_LETTER)

# Parameters
COLORMAP = 'RdBu_r'
VMIN, VMAX = -1, 1
CENTER = 0
MASK_DIAGONAL = True
BOUNDARY_COLOR = 'black'
BOUNDARY_WIDTH = 2
CELL_LABEL_SIZE = 8
COLORBAR_SHRINK = 0.6

#: The size the fifty-six row labels print at on the published page. Held
#: deliberately; see THE ROW LABELS ARE AN EXEMPTION above.
ROW_LABEL_PT = 3.322

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
    'C3_Mac_Inflam_IL1B':             'C3_MoMac_IL1B',
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

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = ROW_LABEL_PT / CELL_LABEL_SIZE
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")
    print(f"  row labels held at {ROW_LABEL_PT:g} pt by the author's exemption")

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
    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

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
        linewidths=0.5 * MARK,
        linecolor='lightgray'
    )
    # De-rasterize heatmap (seaborn uses pcolormesh which defaults to rasterized in SVG)
    for coll in ax.collections:
        coll.set_rasterized(False)

    # Adjust cell labels
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=ROW_LABEL_PT)

    # Add module boundaries
    for boundary in module_boundaries:
        ax.axhline(y=boundary, color=BOUNDARY_COLOR, linewidth=BOUNDARY_WIDTH * MARK)
        ax.axvline(x=boundary, color=BOUNDARY_COLOR, linewidth=BOUNDARY_WIDTH * MARK)

    # Add module labels at bottom
    module_positions = []
    cumsum = 0
    for count in module_counts:
        module_positions.append(cumsum + count / 2)
        cumsum += count

    ax.set_xticks(module_positions)
    ax.set_xticklabels(
        ['IM-T/NK/DC', 'IM-MoMac', 'IM-Mixed', 'IM-Neutrophil', 'IM-B/Plasma'])
    ax.xaxis.tick_bottom()

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    stem = 'correlation_heatmap_k5'
    style.save_panel(fig, script_dir / stem)
    print(f"\nSaved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")
    print("=" * 60)


if __name__ == "__main__":
    main()
