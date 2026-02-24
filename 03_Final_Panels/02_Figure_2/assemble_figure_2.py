#!/usr/bin/env python3
"""
Assemble Figure 2 — CEACAM5/6 Intrinsic Resistance (15 panels A–Q)

Layout matches manuscript reference (02212026):

  Row 1:  A (large UMAP) | C/D stacked | G/H stacked | B (large Milo)
  Row 2:  I (Jaccard heatmap, triangular) — J/K overlap upper-right dead space
  Row 3:  M (lollipop, ~60%) | N1/N2 over O1/O2 (2×2 grid, ~40%)
  Row 4:  P (IHC images, ~78%) | Q (IHC boxplot)

Panels E (Tumor Score) and F (CNV Score) removed per manuscript layout.
J, K overlap I's empty upper-right triangle to save vertical space.

Panel mapping (panel letter → source directory):
  A   Epithelial UMAP (minor states)         02_A
  B   Milo differential abundance             02_F
  C   CEACAM5 expression UMAP                 02_B
  D   CEACAM6 expression UMAP                 02_C
  G   CEACAM5 vs CEACAM6 correlation          02_D
  H   CEACAM5/6 proportion boxplot            02_E
  I   NMF Jaccard heatmap (stomach)           02_P
  J   MP4/5 gene loading dotplot              02_H
  K   MP4/5 boxplots (combined S-MP4 + S-MP5) 02_G
  M   Checkpoint lollipop                     02_I
  N1  CEACAM6 R vs NR boxplot                 02_J
  N2  CEACAM5 R vs NR boxplot                 02_K
  O1  CEACAM6 external validation             02_N
  O2  CEACAM5 external validation             02_O
  P   IHC representative images               02_L
  Q   IHC combined boxplot                    02_M
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from shared.svg_assembler import VectorAssembler

NC_DOUBLE_COL_MM = 183
NC_MAX_HEIGHT_MM = 247

# ── Canvas ──────────────────────────────────────────────────────────
FIG_WIDTH_MM = 180
GAP = 3  # mm gap between rows

# Row heights
R1_H = 65   # Row 1: A | C/D | G/H | B
R2_H = 70   # Row 2: I (Jaccard) with J/K overlapping
R3_H = 42   # Row 3: M | N/O grid
R4_H = 32   # Row 4: P | Q

TOTAL_HEIGHT_MM = R1_H + R2_H + R3_H + R4_H + 3 * GAP  # 218 mm

BASE_DIR = Path(__file__).parent

# Panel label → (source_dir, filename)
PANELS = {
    'A':  ('02_A',  'epithelial_umap_minor_states.svg'),
    'B':  ('02_F',  'epithelial_milo_pre_rvsnr.png'),       # PNG (R-generated)
    'C':  ('02_B',  'ceacam5_expression_umap.svg'),
    'D':  ('02_C',  'ceacam6_expression_umap.svg'),
    'G':  ('02_D',  'ceacam_correlation_scatter.svg'),
    'H':  ('02_E',  'c2_proportion_pre_boxplot.svg'),
    'I':  ('02_P',  'stomach_jaccard_heatmap.svg'),
    'J':  ('02_H',  'mp45_gene_dotplot.svg'),
    'K':  ('02_G',  'mp45_horizontal_boxplot.svg'),
    'M':  ('02_I',  'checkpoint_dotplot_flipped.svg'),
    'N1': ('02_J',  'ceacam6_pre_boxplot.svg'),
    'N2': ('02_K',  'ceacam5_pre_boxplot.svg'),
    'O1': ('02_N',  'ceacam6_prjeb25780_boxplot.svg'),
    'O2': ('02_O',  'ceacam5_prjeb25780_boxplot.svg'),
    'P':  ('02_L',  'ihc_representative_1x6.svg'),
    'Q':  ('02_M',  'ihc_combined_boxplot.svg'),
}


def _path(label):
    folder, fname = PANELS[label]
    return BASE_DIR / folder / fname


def assemble_figure():
    """Assemble Figure 2 with reference manuscript layout."""
    print(f"Assembling Figure 2 (15 panels)...")
    print(f"  Canvas: {FIG_WIDTH_MM} x {TOTAL_HEIGHT_MM} mm")

    asm = VectorAssembler(FIG_WIDTH_MM, TOTAL_HEIGHT_MM)

    # ── Row 1 (y=0): A | C/D stacked | G/H stacked | B ────────────
    y1 = 0
    # Column widths: A(57) gap(2) C/D(22) gap(2) G/H(22) gap(2) B(73) = 180
    w_a = 57
    w_cd = 22
    w_gh = 22
    w_b = 73
    g = 2  # inner gap for row 1

    x_a = 0
    x_cd = w_a + g
    x_gh = x_cd + w_cd + g
    x_b = x_gh + w_gh + g

    # Stacked panel heights within row 1
    h_half = (R1_H - g) / 2  # 31mm each

    # A — large epithelial UMAP
    asm.place_panel('A', _path('A'),
                    x_mm=x_a, y_mm=y1, w_mm=w_a, h_mm=R1_H)

    # C above D (CEACAM5 / CEACAM6 UMAPs)
    asm.place_panel('C', _path('C'),
                    x_mm=x_cd, y_mm=y1, w_mm=w_cd, h_mm=h_half)
    asm.place_panel('D', _path('D'),
                    x_mm=x_cd, y_mm=y1 + h_half + g, w_mm=w_cd, h_mm=h_half)

    # G above H (correlation scatter / proportion boxplot)
    asm.place_panel('G', _path('G'),
                    x_mm=x_gh, y_mm=y1, w_mm=w_gh, h_mm=h_half)
    asm.place_panel('H', _path('H'),
                    x_mm=x_gh, y_mm=y1 + h_half + g, w_mm=w_gh, h_mm=h_half)

    # B — Milo differential abundance (large)
    asm.place_panel('B', _path('B'),
                    x_mm=x_b, y_mm=y1, w_mm=w_b, h_mm=R1_H)

    # ── Row 2 (y=y2): I (Jaccard) with J/K overlapping upper-right ─
    y2 = R1_H + GAP

    # I — NMF Jaccard heatmap (triangular, lower-left has data)
    # Place as a large square; upper-right triangle is empty white space
    w_i = 95
    asm.place_panel('I', _path('I'),
                    x_mm=0, y_mm=y2, w_mm=w_i, h_mm=R2_H)

    # J — NMF gene loading dotplot (overlaps I's empty upper-right)
    x_j = 46
    w_j = FIG_WIDTH_MM - x_j
    h_j = 30
    asm.place_panel('J', _path('J'),
                    x_mm=x_j, y_mm=y2, w_mm=w_j, h_mm=h_j)

    # K — combined MP4/5 boxplots (below J, overlapping I's right side)
    x_k = 98
    y_k = y2 + h_j + 2
    w_k = FIG_WIDTH_MM - x_k
    h_k = R2_H - h_j - 2
    asm.place_panel('K', _path('K'),
                    x_mm=x_k, y_mm=y_k, w_mm=w_k, h_mm=h_k)

    # ── Row 3 (y=y3): M (lollipop) | N1/N2 over O1/O2 ────────────
    y3 = y2 + R2_H + GAP

    # M — checkpoint lollipop (wide, ~60%)
    w_m = 108
    asm.place_panel('M', _path('M'),
                    x_mm=0, y_mm=y3, w_mm=w_m, h_mm=R3_H)

    # 2×2 grid: N1/N2 on top, O1/O2 on bottom
    x_grid = w_m + GAP
    w_cell = (FIG_WIDTH_MM - x_grid - g) / 2  # ~34mm each
    h_top = (R3_H - g) / 2  # ~19.5mm each

    asm.place_panel('N1', _path('N1'),
                    x_mm=x_grid, y_mm=y3, w_mm=w_cell, h_mm=h_top)
    asm.place_panel('N2', _path('N2'),
                    x_mm=x_grid + w_cell + g, y_mm=y3, w_mm=w_cell, h_mm=h_top)
    asm.place_panel('O1', _path('O1'),
                    x_mm=x_grid, y_mm=y3 + h_top + g, w_mm=w_cell, h_mm=h_top)
    asm.place_panel('O2', _path('O2'),
                    x_mm=x_grid + w_cell + g, y_mm=y3 + h_top + g,
                    w_mm=w_cell, h_mm=h_top)

    # ── Row 4 (y=y4): P (IHC images) | Q (IHC boxplot) ───────────
    y4 = y3 + R3_H + GAP

    w_p = 140
    w_q = FIG_WIDTH_MM - w_p - GAP
    asm.place_panel('P', _path('P'),
                    x_mm=0, y_mm=y4, w_mm=w_p, h_mm=R4_H)
    asm.place_panel('Q', _path('Q'),
                    x_mm=w_p + GAP, y_mm=y4, w_mm=w_q, h_mm=R4_H)

    # Save
    asm.save('Figure_02_merged', output_dir=BASE_DIR)

    print(f"\nNature Cancer max: {NC_DOUBLE_COL_MM} x {NC_MAX_HEIGHT_MM} mm")
    print(f"Figure size: {FIG_WIDTH_MM} x {TOTAL_HEIGHT_MM} mm")
    if FIG_WIDTH_MM <= NC_DOUBLE_COL_MM and TOTAL_HEIGHT_MM <= NC_MAX_HEIGHT_MM:
        print("Dimensions within Nature Cancer limits")
    else:
        print("WARNING: Dimensions exceed limits!")


if __name__ == '__main__':
    assemble_figure()
