#!/usr/bin/env python3
"""
Assemble Figure 3 — CEACAM Correlations + CD8 UMAP + Spatial Analysis (14 panels A-N)

Layout (4 rows) — panel letters match Figure Legends:
  Row 1-2 left:  A (CD274 scRNA) | D (Tex scRNA)   }  C (CD8 UMAP, right,
                 B (CD274 TCGA)  | E (Tex TCGA)    }   spans rows 1-2)
  Row 3: F (CEACAM ratio) | G (Epi density) | H (box epi) | I (box stroma) | J (box immune)
  Row 4: K (dist stroma) | L (dist immune) | M (MoMac frac) | N (bar, 2× wide)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from shared.svg_assembler import VectorAssembler

NC_DOUBLE_COL_MM = 183
NC_MAX_HEIGHT_MM = 247

# Figure dimensions (mm)
FIG_WIDTH_MM = 180
GAP_MM = 3

# Row heights
ROW1_HEIGHT_MM = 40     # A + B scRNA correlation scatters
ROW2_HEIGHT_MM = 40     # C + D TCGA validation scatters
ROW3_HEIGHT_MM = 42     # 5 equal panels: 2 spatial + 3 boxplots
ROW4_HEIGHT_MM = 42     # 3 spatial + 1 double-width barplot
TOTAL_HEIGHT_MM = ROW1_HEIGHT_MM + ROW2_HEIGHT_MM + ROW3_HEIGHT_MM + ROW4_HEIGHT_MM + 3 * GAP_MM

# Rows 1-2: UMAP on right, scatters on left
UMAP_HEIGHT_MM = ROW1_HEIGHT_MM + ROW2_HEIGHT_MM + GAP_MM  # spans both rows
UMAP_WIDTH_MM = 60      # narrower UMAP, gives more room to scatters
SCATTER_ZONE_MM = FIG_WIDTH_MM - UMAP_WIDTH_MM - GAP_MM   # remaining for scatters
SCATTER_HALF_MM = (SCATTER_ZONE_MM - GAP_MM) / 2           # each scatter panel

# Rows 3-4: 5-column grid
UNIT_WIDTH_MM = (FIG_WIDTH_MM - 4 * GAP_MM) / 5  # 5 panels, 4 gaps
DOUBLE_UNIT_MM = 2 * UNIT_WIDTH_MM + GAP_MM       # N spans 2 units

BASE_DIR = Path(__file__).parent

PANELS = {
    # Rows 1-2 left: scRNA + TCGA correlations (legend order: A,B=CD274 pair; D,E=Tex pair)
    'A': ('03_A', 'ceacam_cd274_scatter.png'),
    'D': ('03_B', 'ceacam_tex_scatter.png'),
    'B': ('03_C', 'tcga_ceacam_cd274.png'),
    'E': ('03_D', 'tcga_ceacam_tex.png'),
    # Rows 1-2 right: CD8 UMAP
    'C': ('03_E', 'cd8_umap_minor_states.png'),
    # Row 3: spatial panels
    'F': ('03_F', 'spatial_ceacam_ratio.png'),
    'G': ('03_G', 'spatial_epi_density.png'),
    'H': ('03_H', 'ceacam_spatial_boxplot.png'),
    'I': ('03_I', 'ceacam_spatial_boxplot.png'),
    'J': ('03_J', 'ceacam_spatial_boxplot.png'),
    # Row 4: spatial + barplot
    'K': ('03_K', 'spatial_dist_stroma.png'),
    'L': ('03_L', 'spatial_dist_immune.png'),
    'M': ('03_M', 'spatial_momac_fraction.png'),
    'N': ('03_N', 'ceacam_immune_region_summary.png'),
}


def assemble_figure():
    """Assemble Figure 3 with 14 panels (A-N)."""
    print("Assembling Figure 3 (14 panels)...")
    print(f"  Canvas: {FIG_WIDTH_MM} x {TOTAL_HEIGHT_MM} mm")
    print(f"  Scatter zone: {SCATTER_ZONE_MM:.1f} mm, each half: {SCATTER_HALF_MM:.1f} mm")
    print(f"  UMAP: {UMAP_WIDTH_MM} x {UMAP_HEIGHT_MM} mm")
    print(f"  Unit width: {UNIT_WIDTH_MM:.1f} mm, Double: {DOUBLE_UNIT_MM:.1f} mm")

    asm = VectorAssembler(FIG_WIDTH_MM, TOTAL_HEIGHT_MM)

    # Y positions
    y_row1 = 0
    y_row2 = ROW1_HEIGHT_MM + GAP_MM
    y_row3 = ROW1_HEIGHT_MM + ROW2_HEIGHT_MM + 2 * GAP_MM
    y_row4 = ROW1_HEIGHT_MM + ROW2_HEIGHT_MM + ROW3_HEIGHT_MM + 3 * GAP_MM

    # X positions for scatter columns
    x_scatter_left = 0
    x_scatter_right = SCATTER_HALF_MM + GAP_MM
    x_umap = SCATTER_ZONE_MM + GAP_MM

    # ---- Row 1: A (CD274 scRNA) + D (Tex scRNA) ----
    asm.place_panel('A', BASE_DIR / PANELS['A'][0] / PANELS['A'][1],
                    x_mm=x_scatter_left, y_mm=y_row1,
                    w_mm=SCATTER_HALF_MM, h_mm=ROW1_HEIGHT_MM)

    asm.place_panel('D', BASE_DIR / PANELS['D'][0] / PANELS['D'][1],
                    x_mm=x_scatter_right, y_mm=y_row1,
                    w_mm=SCATTER_HALF_MM, h_mm=ROW1_HEIGHT_MM)

    # ---- Row 2: B (CD274 TCGA) + E (Tex TCGA) ----
    asm.place_panel('B', BASE_DIR / PANELS['B'][0] / PANELS['B'][1],
                    x_mm=x_scatter_left, y_mm=y_row2,
                    w_mm=SCATTER_HALF_MM, h_mm=ROW2_HEIGHT_MM)

    asm.place_panel('E', BASE_DIR / PANELS['E'][0] / PANELS['E'][1],
                    x_mm=x_scatter_right, y_mm=y_row2,
                    w_mm=SCATTER_HALF_MM, h_mm=ROW2_HEIGHT_MM)

    # ---- UMAP: C spans rows 1-2 on the right ----
    asm.place_panel('C', BASE_DIR / PANELS['C'][0] / PANELS['C'][1],
                    x_mm=x_umap, y_mm=y_row1,
                    w_mm=UMAP_WIDTH_MM, h_mm=UMAP_HEIGHT_MM)

    # ---- Row 3: 5 equal columns — F, G, H, I, J ----
    for i, label in enumerate(['F', 'G', 'H', 'I', 'J']):
        x = i * (UNIT_WIDTH_MM + GAP_MM)
        asm.place_panel(label, BASE_DIR / PANELS[label][0] / PANELS[label][1],
                        x_mm=x, y_mm=y_row3, w_mm=UNIT_WIDTH_MM, h_mm=ROW3_HEIGHT_MM)

    # ---- Row 4: K, L, M (1 unit each) + N (2 units) ----
    for i, label in enumerate(['K', 'L', 'M']):
        x = i * (UNIT_WIDTH_MM + GAP_MM)
        asm.place_panel(label, BASE_DIR / PANELS[label][0] / PANELS[label][1],
                        x_mm=x, y_mm=y_row4, w_mm=UNIT_WIDTH_MM, h_mm=ROW4_HEIGHT_MM)

    # N: double width, starting at column 4
    x_n = 3 * (UNIT_WIDTH_MM + GAP_MM)
    asm.place_panel('N', BASE_DIR / PANELS['N'][0] / PANELS['N'][1],
                    x_mm=x_n, y_mm=y_row4, w_mm=DOUBLE_UNIT_MM, h_mm=ROW4_HEIGHT_MM)

    # Save
    asm.save('Figure_03_merged', output_dir=BASE_DIR)

    print(f"\nNature Cancer max: {NC_DOUBLE_COL_MM} x {NC_MAX_HEIGHT_MM} mm")
    if FIG_WIDTH_MM <= NC_DOUBLE_COL_MM and TOTAL_HEIGHT_MM <= NC_MAX_HEIGHT_MM:
        print("Dimensions within Nature Cancer limits")
    else:
        print("WARNING: Dimensions exceed limits!")


if __name__ == '__main__':
    assemble_figure()
