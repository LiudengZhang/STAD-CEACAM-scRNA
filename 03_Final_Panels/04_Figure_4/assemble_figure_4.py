#!/usr/bin/env python3
"""
Assemble Figure 4 for Nature Cancer submission (8 Panels A-H)
MoMac Module - Acquired Resistance

Layout (3 rows) — panel letters match Figure Legends:
  Row 1: A (correlation heatmap) + B (network graphs)
  Row 2: C (module2 boxplot) + E (marker dotplot)
  Row 3: D (MoMac UMAP) + F (density plots) + G (forest plot) + H (external validation)

Output: Vector PDF via svglib/reportlab + PNG preview.
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from shared.svg_assembler import VectorAssembler

# Nature Cancer specifications
NC_DOUBLE_COL_MM = 183
NC_MAX_HEIGHT_MM = 247

# Figure dimensions (mm)
FIG_WIDTH_MM = 175
GAP_MM = 2.5

ROW1_HEIGHT_MM = 90
ROW2_HEIGHT_MM = 34
ROW3_HEIGHT_MM = 58
TOTAL_HEIGHT_MM = ROW1_HEIGHT_MM + ROW2_HEIGHT_MM + ROW3_HEIGHT_MM + 2 * GAP_MM

BASE_DIR = Path(__file__).parent

# Panel files — .svg preferred, .png fallback handled by VectorAssembler
PANELS = {
    'A': ('04_A', 'correlation_heatmap_k5.svg'),
    'B': ('04_B', 'merged_networks_k5.svg'),
    'C': ('04_C', 'module2_post_response_boxplot_k5.svg'),
    'E': ('04_E', 'momac_marker_dotplot.svg'),
    'D': ('04_E', 'momac_umap_minor_states.svg'),
    'F': ('04_F', 'momac_density_panel.svg'),
    'G': ('04_G', 'momac_forest_plot_final.svg'),
}

# Panel H sub-panels (stacked vertically)
PANEL_H_TOP = ('04_H', 'mac3_boxplot_linghua.svg')
PANEL_H_BOTTOM = ('04_H', 'mac3_boxplot_kumar.svg')


def assemble_figure():
    """Assemble Figure 4 with 8 panels (A-H) as vector PDF."""
    print("Assembling Figure 4 (vector PDF)...")

    asm = VectorAssembler(FIG_WIDTH_MM, TOTAL_HEIGHT_MM)

    # Y positions (from top, in mm)
    y_row1 = 0
    y_row2 = ROW1_HEIGHT_MM + GAP_MM
    y_row3 = ROW1_HEIGHT_MM + GAP_MM + ROW2_HEIGHT_MM + GAP_MM

    # === ROW 1: A (152mm) + B (20.5mm) ===
    w_a, w_b = 152, 20.5
    asm.place_panel('A', BASE_DIR / PANELS['A'][0] / PANELS['A'][1],
                    x_mm=0, y_mm=y_row1, w_mm=w_a, h_mm=ROW1_HEIGHT_MM)
    asm.place_panel('B', BASE_DIR / PANELS['B'][0] / PANELS['B'][1],
                    x_mm=w_a + GAP_MM, y_mm=y_row1, w_mm=w_b, h_mm=ROW1_HEIGHT_MM)

    # === ROW 2: C (22mm) + E (150.5mm, marker dotplot) ===
    w_c, w_e = 22, 150.5
    x_c = 0
    x_e = w_c + GAP_MM
    asm.place_panel('C', BASE_DIR / PANELS['C'][0] / PANELS['C'][1],
                    x_mm=x_c, y_mm=y_row2, w_mm=w_c, h_mm=ROW2_HEIGHT_MM)
    asm.place_panel('E', BASE_DIR / PANELS['E'][0] / PANELS['E'][1],
                    x_mm=x_e, y_mm=y_row2, w_mm=w_e, h_mm=ROW2_HEIGHT_MM)

    # === ROW 3: D (48mm, UMAP) + F (52mm) + G (45mm) + H (22.5mm) ===
    w_d, w_f, w_g, w_h = 48, 52, 45, 22.5
    x_d = 0
    x_f = w_d + GAP_MM
    x_g = x_f + w_f + GAP_MM
    x_h = x_g + w_g + GAP_MM
    asm.place_panel('D', BASE_DIR / PANELS['D'][0] / PANELS['D'][1],
                    x_mm=x_d, y_mm=y_row3, w_mm=w_d, h_mm=ROW3_HEIGHT_MM)
    asm.place_panel('F', BASE_DIR / PANELS['F'][0] / PANELS['F'][1],
                    x_mm=x_f, y_mm=y_row3, w_mm=w_f, h_mm=ROW3_HEIGHT_MM)
    asm.place_panel('G', BASE_DIR / PANELS['G'][0] / PANELS['G'][1],
                    x_mm=x_g, y_mm=y_row3, w_mm=w_g, h_mm=ROW3_HEIGHT_MM)

    # Panel H: two sub-panels stacked vertically
    h_half = (ROW3_HEIGHT_MM - 1) / 2  # 1mm gap between sub-panels
    asm.place_panel('H', BASE_DIR / PANEL_H_TOP[0] / PANEL_H_TOP[1],
                    x_mm=x_h, y_mm=y_row3, w_mm=w_h, h_mm=h_half)
    asm.place_panel(None, BASE_DIR / PANEL_H_BOTTOM[0] / PANEL_H_BOTTOM[1],
                    x_mm=x_h, y_mm=y_row3 + h_half + 1, w_mm=w_h, h_mm=h_half)

    # Save
    asm.save('Figure_04_merged', output_dir=BASE_DIR)

    # Dimension check
    print(f"\nNature Cancer max: {NC_DOUBLE_COL_MM} x {NC_MAX_HEIGHT_MM} mm")
    if FIG_WIDTH_MM <= NC_DOUBLE_COL_MM and TOTAL_HEIGHT_MM <= NC_MAX_HEIGHT_MM:
        print("Dimensions within Nature Cancer limits")
    else:
        print("WARNING: Dimensions exceed limits!")


if __name__ == '__main__':
    assemble_figure()
