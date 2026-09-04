#!/usr/bin/env python3
"""
Figure 4 panel A, RESTYLED (Version B) - Spearman correlation matrix of the 56
minor cell states, ordered and boxed by the five interaction modules.

Version A is
`03_Revised_Panels/Main_Figures/04_Figure_4/04_A/create_correlation_heatmap.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every correlation, every module assignment, every ordering, every short name
and every string is Version A's. The drawing code is the same code.

  printed panel  Figure 4 A     (PROVENANCE.csv; NOT inferred from "04_A")
  printed rect   120.1 x 94.3 mm   (panel_rects.csv)
  Version B box  190.0 x 152.0 mm

MARK
    Version A drew a 12 x 12 inch canvas with no SCALE constant, and set two
    type sizes on it: `CELL_LABEL_SIZE = 8` on the 56 row labels and
    `MODULE_LABEL_SIZE = 12` on the five module names. The assembler fitted the
    saved 812.5 x 650.1 pt SVG into 120.1 x 94.3 mm, a fit of 0.4110, so the row
    labels printed at 3.29 pt. The smallest body type is therefore SMALL_PT = 8
    with SCALE = 1:

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 8 = 0.875

    The two non-type lengths are the module boundary rules (2 pt, printed at
    0.82 pt, now 1.75 pt) and the 0.5 pt grid between cells (now 0.4375 pt).
    Both are written as Version A wrote them with `* MARK` appended.

THE BOX - why this panel runs to the full page width
    This is the one panel in Figure 4 that the page constrains rather than the
    other way round. 56 row labels at 7 pt want about 7.5 pt of pitch each,
    which is 148 mm of height; the published panel gives them 94.3 mm, i.e.
    4.8 pt a row, which is why they were set at 8 pt and printed at 3.29 pt.
    `square=True` then makes the heatmap as wide as it is tall, and seaborn's
    colourbar takes a further 15% of the axes width, so the height that can be
    reached is fixed by the width available: at the full 190 mm page, with
    21 mm of gutter for 'B_Naive_TCL1A' at 7 pt and 4 mm for the colourbar
    label, the heatmap comes out ~139 mm tall, or 7.0 pt a row. That is the
    ceiling, and it is where this panel is drawn: 190 x 152 mm, 0.5 mm inside
    the page. The row labels clear each other at that pitch but not by much.
    The published 1.27 aspect becomes 1.25. PANEL_SPEC anticipates the growth:
    "Panels will grow."

    The margins are set *before* `sns.heatmap` runs. seaborn builds the
    colourbar with `plt.colorbar(..., ax=ax)`, which steals its space out of
    the axes' position at that moment and parks the colourbar axes outside the
    subplot grid; calling `subplots_adjust` afterwards would move the heatmap
    back over it. Same margins, applied in the only order that composes.

THE FRAME IS cnsplots'
    Version A set no frame widths of its own, so nothing is dropped here;
    cnsplots' `axes.linewidth` and tick settings govern, as everywhere else in
    Version B.

Input : this directory / correlation_matrix.csv, module_mappings.csv
        Version A reads these two files from beside its own script. The path is
        part of the code PANEL_SPEC forbids changing, so the two files were
        copied here byte for byte (md5 7471d50f... and 06b2...) rather than the
        read being repointed at the frozen tree. Nothing was moved or deleted.
Output: this directory / correlation_heatmap_k5.{svg,pdf,png}
"""

import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
import panel_style_cns as style              # noqa: E402

# Version A's canvas convention, kept only so MARK can be derived from it.
SCALE = 1
SMALL_PT = 8.0                               # Version A's `CELL_LABEL_SIZE`

PRINTED_MM = (120.1, 94.3)                   # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 190.0, 152.0
MARGIN = dict(left=21.0, right=4.0, top=1.5, bottom=5.5)

# Parameters
DPI = 300
COLORMAP = 'RdBu_r'
VMIN, VMAX = -1, 1
CENTER = 0
MASK_DIAGONAL = True
BOUNDARY_COLOR = 'black'
BOUNDARY_WIDTH = 2
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

    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

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

    # Margins first - see the module docstring. seaborn's colourbar is cut out
    # of the axes' position at draw time and does not move with a later
    # subplots_adjust.
    style.margins_mm(fig, **MARGIN)

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
    ax.set_yticklabels(ax.get_yticklabels())

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
        ['IM-T/NK/DC', 'IM-MoMac', 'IM-Mixed', 'IM-Neutrophil', 'IM-B/Plasma']
    )
    ax.xaxis.tick_bottom()

    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, script_dir / 'correlation_heatmap_k5')
    print(f"\nSaved: {script_dir / 'correlation_heatmap_k5'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")
    print("=" * 60)

if __name__ == "__main__":
    main()
