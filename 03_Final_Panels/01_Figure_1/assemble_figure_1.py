#!/usr/bin/env python3
"""
Assemble Figure 1 for Nature Cancer submission (3 panels, A-C)

Layout (2 rows):
  Row 1: A (Schematic, full width)
  Row 2: B (UMAP, ~45%) + C (Sample balance stacked bar, ~55%)

Output: Vector PDF via svglib/reportlab + PNG preview.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from shared.svg_assembler import VectorAssembler

# Nature Cancer specifications
NC_DOUBLE_COL_MM = 183
NC_MAX_HEIGHT_MM = 247

# Figure dimensions (mm)
FIG_WIDTH_MM = 180
GAP_MM = 3

# Row heights (mm)
ROW1_HEIGHT_MM = 125  # A schematic
ROW2_HEIGHT_MM = 48   # B + C
TOTAL_HEIGHT_MM = ROW1_HEIGHT_MM + ROW2_HEIGHT_MM + GAP_MM

BASE_DIR = Path(__file__).parent

# Panel files
PANELS = {
    'A': ('01_A', 'image.png'),
    'B': ('01_B', 'umap_major_cell_type_stomach.svg'),
    'C': ('01_C', '01_C_sample_balance_stomach.svg'),
}


def assemble_figure():
    """Assemble Figure 1 with 3 panels (A-C) as vector PDF."""
    print("Assembling Figure 1 (vector PDF)...")

    asm = VectorAssembler(FIG_WIDTH_MM, TOTAL_HEIGHT_MM)

    # Y positions (from top, mm)
    y_row1 = 0
    y_row2 = ROW1_HEIGHT_MM + GAP_MM

    # Row 1: A (schematic, full width)
    asm.place_panel('A', BASE_DIR / PANELS['A'][0] / PANELS['A'][1],
                    x_mm=0, y_mm=y_row1, w_mm=FIG_WIDTH_MM, h_mm=ROW1_HEIGHT_MM)

    # Row 2: B (45%) + C (55%)
    w_b = FIG_WIDTH_MM * 0.44
    w_c = FIG_WIDTH_MM * 0.54
    asm.place_panel('B', BASE_DIR / PANELS['B'][0] / PANELS['B'][1],
                    x_mm=0, y_mm=y_row2, w_mm=w_b, h_mm=ROW2_HEIGHT_MM)
    asm.place_panel('C', BASE_DIR / PANELS['C'][0] / PANELS['C'][1],
                    x_mm=w_b + GAP_MM, y_mm=y_row2, w_mm=w_c, h_mm=ROW2_HEIGHT_MM)

    # Save
    asm.save('Figure_01_merged', output_dir=BASE_DIR)

    # Dimension check
    print(f"\nNature Cancer max: {NC_DOUBLE_COL_MM} x {NC_MAX_HEIGHT_MM} mm")
    if FIG_WIDTH_MM <= NC_DOUBLE_COL_MM and TOTAL_HEIGHT_MM <= NC_MAX_HEIGHT_MM:
        print("Dimensions within Nature Cancer limits")
    else:
        print("WARNING: Dimensions exceed limits!")


if __name__ == '__main__':
    assemble_figure()
