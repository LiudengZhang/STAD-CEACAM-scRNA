#!/usr/bin/env python3
"""Figure 4A: Spearman correlation between every pair of minor cell states,
ordered by module, with the five module names along the bottom.
The 56 row labels are 4 pt (ROW_LABEL_PT) so that each fits its row.

Input:  correlation_matrix.csv, module_mappings.csv (beside this script)
Output: correlation_heatmap_k5.{svg,pdf,png} in this directory
Run:    python create_correlation_heatmap.py
"""

import pandas as pd
import numpy as np
import re
import matplotlib
import matplotlib.transforms as mtransforms
matplotlib.use('Agg')
import seaborn as sns
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from cnsfig import Panel, style  # noqa: E402

ROW_LABEL_PT = 4.0                # row label size, below the 6 pt body size

PANEL = Panel(4, "A", __file__, small_pt=ROW_LABEL_PT)

COLORMAP = 'RdBu_r'
VMIN, VMAX = -1, 1
CENTER = 0
MASK_DIAGONAL = True
BOUNDARY_COLOR = 'black'
LABEL_SPLIT_MM = 0.75             # IM-Neutrophil left, IM-B/Plasma right (see main)
# Colour bar size in mm, centred on the heatmap.
CBAR_W_MM, CBAR_H_MM, CBAR_GAP_MM = 1.4, 36.0, 1.5

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
    PANEL.apply()

    correlation_matrix = pd.read_csv(PANEL.dir / 'correlation_matrix.csv', index_col=0)
    module_mappings = pd.read_csv(PANEL.dir / 'module_mappings.csv')

    module_mappings_sorted = module_mappings.sort_values('module')
    ordered_states = module_mappings_sorted['minor_cell_state'].tolist()
    correlation_ordered = correlation_matrix.loc[ordered_states, ordered_states]

    # Drop the Cx_ prefix unless that would make names collide.
    strip_cx = lambda s: re.sub(r'^C\d+_', '', s)
    short_index = [strip_cx(SHORT_NAMES.get(s, s)) for s in correlation_ordered.index]
    dup = {n for n in short_index if short_index.count(n) > 1}
    short_index = [SHORT_NAMES.get(s, s) if strip_cx(SHORT_NAMES.get(s, s)) in dup
                   else strip_cx(SHORT_NAMES.get(s, s))
                   for s in correlation_ordered.index]
    correlation_ordered.index = short_index
    correlation_ordered.columns = short_index

    module_counts = module_mappings_sorted.groupby('module').size()
    module_boundaries = np.cumsum(module_counts.values)[:-1]

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)
    cax = fig.add_axes([0.9, 0.5, 0.01, 0.2])     # placed below, after the fit

    mask = np.eye(len(correlation_ordered), dtype=bool) if MASK_DIAGONAL else None

    sns.heatmap(
        correlation_ordered,
        cmap=COLORMAP,
        vmin=VMIN,
        vmax=VMAX,
        center=CENTER,
        mask=mask,
        square=True,
        cbar_ax=cax,
        cbar_kws={'label': 'Spearman Correlation'},
        ax=ax,
        xticklabels=False,
        yticklabels=True,
        linewidths=style.EDGE_PT,
        linecolor='lightgray'
    )
    # De-rasterize heatmap (seaborn uses pcolormesh which defaults to rasterized in SVG)
    for coll in ax.collections:
        coll.set_rasterized(False)

    ax.set_yticklabels(short_index, fontsize=ROW_LABEL_PT, rotation=0)
    ax.tick_params(axis='y', length=1.0, width=style.RULE_PT, pad=1.0)

    for boundary in module_boundaries:
        ax.axhline(y=boundary, color=BOUNDARY_COLOR, linewidth=style.RULE_PT)
        ax.axvline(x=boundary, color=BOUNDARY_COLOR, linewidth=style.RULE_PT)

    module_positions = []
    cumsum = 0
    for count in module_counts:
        module_positions.append(cumsum + count / 2)
        cumsum += count

    ax.set_xticks(module_positions)
    ax.set_xticklabels(
        ['IM-T/NK/DC', 'IM-MoMac', 'IM-Mixed', 'IM-Neutrophil', 'IM-B/Plasma'])
    ax.xaxis.tick_bottom()
    # Shift the IM-Neutrophil and IM-B/Plasma labels apart so they do not touch.
    labels = ax.get_xticklabels()
    for lab, sign in ((labels[3], -1), (labels[4], +1)):
        lab.set_transform(lab.get_transform() + mtransforms.ScaledTranslation(
            sign * LABEL_SPLIT_MM / 25.4, 0, fig.dpi_scale_trans))

    cax.tick_params(labelsize=style.tick_pt(), width=style.RULE_PT, length=1.5)
    cax.yaxis.label.set_size(style.tick_pt())
    PANEL.fit(fig)
    # square=True shrinks the axes at draw time, so place the bar after a draw.
    fig.canvas.draw()
    W, H = PANEL.w_mm, PANEL.h_mm
    bb = ax.get_window_extent(fig.canvas.get_renderer()).transformed(
        fig.transFigure.inverted())
    cax.set_position([bb.x1 + CBAR_GAP_MM / W,
                      (bb.y0 + bb.y1) / 2 - CBAR_H_MM / H / 2,
                      CBAR_W_MM / W, CBAR_H_MM / H])
    PANEL.finish(fig, 'correlation_heatmap_k5')


if __name__ == "__main__":
    main()
