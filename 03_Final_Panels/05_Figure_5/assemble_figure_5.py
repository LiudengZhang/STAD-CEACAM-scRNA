#!/usr/bin/env python3
"""
Assemble Figure 5 — MoMac NF-κB Signaling & Downstream Consequences (4-row layout)

Row 1 (38mm): A (MoMac barplot) | B (pathway UMAP) | C (gene UMAP) | D (BACH1 TF)
Row 2 (42mm): E (NFKB1 TF) | F (cytokine dotplot) | G (violin) | H (radar) | I (GSEA)
Row 3 (36mm): J (CD274 MoMac) | K (CD274 Epi) | L (CD274 Fib) | M (CD274 DC) |
              N (Tex CD8) | O (IL-6 CD4) | P (Th17 CD4)
Row 4 (46mm): Q1 (GSEA MoMac) | Q2 (GSEA Epi) | Q3 (GSEA Fib) | Q4 (GSEA DC) | Legend

Inter-row gap: 2mm.  Total: 38+2+42+2+36+2+46 = 168mm (≤170mm limit).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from shared.svg_assembler import VectorAssembler

NC_DOUBLE_COL_MM = 183
NC_MAX_HEIGHT_MM = 170

FIG_WIDTH_MM = 180
ROW_GAP = 2       # inter-row gap (compressed from 3)
PANEL_GAP = 3     # intra-row gap between panels

ROW1_HEIGHT_MM = 38
ROW2_HEIGHT_MM = 42
ROW3_HEIGHT_MM = 36
ROW4_HEIGHT_MM = 46
TOTAL_HEIGHT_MM = (ROW1_HEIGHT_MM + ROW2_HEIGHT_MM + ROW3_HEIGHT_MM + ROW4_HEIGHT_MM
                   + 3 * ROW_GAP)  # 168mm

BASE_DIR = Path(__file__).parent

# Map: figure label → (folder, filename)
PANELS = {
    # Row 1 — Pathway Enrichment + Gene UMAPs + BACH1
    'A': ('05_A',  'momac_pathway_enrichment.svg'),
    'B': ('05_B',  'momac_4pathway_umap.svg'),
    'C': ('05_D',  'tnf_il1b_il6_il1a_cytokines.svg'),
    'D': ('05_TF', 'bach1_tf_4group.svg'),
    # Row 2 — NFKB1 + Cytokine Dotplot + Violin + NF-kB Overview
    'E': ('05_TF', 'nfkb1_tf_4group.svg'),
    'F': ('05_H',  'cytokine_dotplot.svg'),
    'G': ('05_C',  'macrophage_top_ligands_violin.svg'),
    'H': ('05_F',  'nfkb_radar_celltype_enrichment.svg'),
    'I': ('05_G',  'panel_g_gsea_4types.svg'),
    # Row 3 — CD274 boxplots + T-cell state scatters + IL-6 score
    'J': ('05_I',        'cd274_momac_boxplot.svg'),
    'K': ('05_J',        'cd274_epithelial_boxplot.svg'),
    'L': ('05_K',        'cd274_fibroblast_boxplot.svg'),
    'M': ('05_DC_CD274', 'cd274_dc_boxplot.svg'),
    'N': ('05_O',        'tex_nfkb_scatter_cd8.svg'),
    'O': ('05_IL6_CD4',  'il6_stat3_cd4t_boxplot.svg'),
    'P': ('05_Q',        'th17_nfkb_scatter_cd4.svg'),
}

# Row 4 GSEA panels (separate dict — no single-letter labels)
GSEA_PANELS = {
    'Q1': ('05_GSEA_Summary', 'gsea_dotplot_momac.svg'),
    'Q2': ('05_GSEA_Summary', 'gsea_dotplot_epithelial.svg'),
    'Q3': ('05_GSEA_Summary', 'gsea_dotplot_fibroblast.svg'),
    'Q4': ('05_GSEA_Summary', 'gsea_dotplot_dc.svg'),
    'QL': ('05_GSEA_Summary', 'gsea_legend.svg'),
}


def _panel_path(label):
    folder, fname = PANELS[label]
    return BASE_DIR / folder / fname


def _gsea_path(label):
    folder, fname = GSEA_PANELS[label]
    return BASE_DIR / folder / fname


def assemble_figure():
    """Assemble Figure 5 as vector SVG/PDF/PNG."""
    n_panels = len(PANELS) + len(GSEA_PANELS)
    print(f"Assembling Figure 5 — {n_panels} panels, 4 rows...")

    asm = VectorAssembler(FIG_WIDTH_MM, TOTAL_HEIGHT_MM)

    y_r1 = 0
    y_r2 = ROW1_HEIGHT_MM + ROW_GAP
    y_r3 = y_r2 + ROW2_HEIGHT_MM + ROW_GAP
    y_r4 = y_r3 + ROW3_HEIGHT_MM + ROW_GAP

    # === ROW 1: A | B | C | D  (4 panels) ===
    w_a, w_b, w_c, w_d = 42, 48, 46, 35
    g = PANEL_GAP
    x_a = 0
    x_b = w_a + g
    x_c = x_b + w_b + g
    x_d = x_c + w_c + g

    asm.place_panel('A', _panel_path('A'),
                    x_mm=x_a, y_mm=y_r1, w_mm=w_a, h_mm=ROW1_HEIGHT_MM)
    asm.place_panel('B', _panel_path('B'),
                    x_mm=x_b, y_mm=y_r1, w_mm=w_b, h_mm=ROW1_HEIGHT_MM,
                    force_png=True)
    asm.place_panel('C', _panel_path('C'),
                    x_mm=x_c, y_mm=y_r1, w_mm=w_c, h_mm=ROW1_HEIGHT_MM,
                    force_png=True)
    asm.place_panel('D', _panel_path('D'),
                    x_mm=x_d, y_mm=y_r1, w_mm=w_d, h_mm=ROW1_HEIGHT_MM)

    # === ROW 2: E | F | G | H | I  (5 panels) ===
    w_e, w_f, w_g, w_h, w_i = 30, 35, 25, 37, 41
    x_e = 0
    x_f = w_e + g
    x_g = x_f + w_f + g
    x_h = x_g + w_g + g
    x_i = x_h + w_h + g

    asm.place_panel('E', _panel_path('E'),
                    x_mm=x_e, y_mm=y_r2, w_mm=w_e, h_mm=ROW2_HEIGHT_MM)
    asm.place_panel('F', _panel_path('F'),
                    x_mm=x_f, y_mm=y_r2, w_mm=w_f, h_mm=ROW2_HEIGHT_MM)
    asm.place_panel('G', _panel_path('G'),
                    x_mm=x_g, y_mm=y_r2, w_mm=w_g, h_mm=ROW2_HEIGHT_MM)
    asm.place_panel('H', _panel_path('H'),
                    x_mm=x_h, y_mm=y_r2, w_mm=w_h, h_mm=ROW2_HEIGHT_MM)
    asm.place_panel('I', _panel_path('I'),
                    x_mm=x_i, y_mm=y_r2, w_mm=w_i, h_mm=ROW2_HEIGHT_MM)

    # === ROW 3: J|K|L|M (boxplots) | N (scatter) | O (boxplot) | P (scatter) ===
    # 4 boxplots @21mm + 1 scatter @27mm + 1 boxplot @21mm + 1 scatter @27mm
    # + 6 gaps @3mm = 21*4 + 27*2 + 21 + 6*3 = 84+54+21+18 = 177mm
    asm.place_panel('J', _panel_path('J'),
                    x_mm=0, y_mm=y_r3, w_mm=21, h_mm=ROW3_HEIGHT_MM)
    asm.place_panel('K', _panel_path('K'),
                    x_mm=24, y_mm=y_r3, w_mm=21, h_mm=ROW3_HEIGHT_MM)
    asm.place_panel('L', _panel_path('L'),
                    x_mm=48, y_mm=y_r3, w_mm=21, h_mm=ROW3_HEIGHT_MM)
    asm.place_panel('M', _panel_path('M'),
                    x_mm=72, y_mm=y_r3, w_mm=21, h_mm=ROW3_HEIGHT_MM)
    asm.place_panel('N', _panel_path('N'),
                    x_mm=96, y_mm=y_r3, w_mm=27, h_mm=ROW3_HEIGHT_MM)
    asm.place_panel('O', _panel_path('O'),
                    x_mm=126, y_mm=y_r3, w_mm=21, h_mm=ROW3_HEIGHT_MM)
    asm.place_panel('P', _panel_path('P'),
                    x_mm=150, y_mm=y_r3, w_mm=27, h_mm=ROW3_HEIGHT_MM)

    # === ROW 4: 4 GSEA dotplots + legend ===
    # 4 × 38mm + 3 × 2mm gaps + 18mm legend = 176mm
    gsea_w = 38
    gsea_gap = 2
    legend_w = 18
    x = 0
    for label in ['Q1', 'Q2', 'Q3', 'Q4']:
        asm.place_panel(None, _gsea_path(label),
                        x_mm=x, y_mm=y_r4, w_mm=gsea_w, h_mm=ROW4_HEIGHT_MM)
        x += gsea_w + gsea_gap
    asm.place_panel(None, _gsea_path('QL'),
                    x_mm=x, y_mm=y_r4, w_mm=legend_w, h_mm=ROW4_HEIGHT_MM)

    # Save
    asm.save('Figure_05_merged', output_dir=BASE_DIR)

    print(f"\nCancer Discovery max: {NC_DOUBLE_COL_MM} x {NC_MAX_HEIGHT_MM} mm")
    print(f"Figure size: {FIG_WIDTH_MM} x {TOTAL_HEIGHT_MM} mm")
    if FIG_WIDTH_MM <= NC_DOUBLE_COL_MM and TOTAL_HEIGHT_MM <= NC_MAX_HEIGHT_MM:
        print("Dimensions within limits")
    else:
        print("WARNING: Dimensions exceed limits!")


if __name__ == '__main__':
    assemble_figure()
