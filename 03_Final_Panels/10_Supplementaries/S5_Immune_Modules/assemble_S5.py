#!/usr/bin/env python3
"""
Assemble Supplementary Figure S5: Immune Module Comparisons

Layout (2 rows x 3 columns, all panels same size):
  Row 1: A (M1 Adaptive) + B (M2 MoMac) + C (M3 Mixed)
  Row 2: D (M4 Neutrophils) + E (M5 B/Plasma) + [blank]

Output: Vector SVG + PDF + PNG preview.
"""

import sys
from pathlib import Path
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from shared.svg_assembler import VectorAssembler

FIG_WIDTH_MM = 183
GAP_MM = 3

SCRIPT_DIR = Path(__file__).parent
PANEL_A = SCRIPT_DIR / "S5_A" / "panel_S5_A.svg"
PANEL_B = SCRIPT_DIR / "S5_B" / "panel_S5_B.svg"
PANEL_C = SCRIPT_DIR / "S5_C" / "panel_S5_C.svg"
PANEL_D = SCRIPT_DIR / "S5_D" / "panel_S5_D.svg"
PANEL_E = SCRIPT_DIR / "S5_E" / "panel_S5_E.svg"


def get_height_mm(panel_path, width_mm):
    png = panel_path.with_suffix('.png')
    if png.exists():
        img = Image.open(png)
        return width_mm * (img.height / img.width)
    return width_mm * 0.5


def assemble_figure():
    print("Assembling Supplementary Figure S5: Immune Module Comparisons")

    all_panels = [
        ('A', PANEL_A), ('B', PANEL_B), ('C', PANEL_C),
        ('D', PANEL_D), ('E', PANEL_E),
    ]
    for name, path in all_panels:
        svg_ok = path.exists()
        png_ok = path.with_suffix('.png').exists()
        print(f"  Panel {name}: SVG={'OK' if svg_ok else 'MISSING'}, PNG={'OK' if png_ok else 'MISSING'}")

    # 2 rows x 3 columns — all panels same size, position (2,3) blank
    w_third = (FIG_WIDTH_MM - 2 * GAP_MM) / 3

    # Row 1: A + B + C
    h_a = get_height_mm(PANEL_A, w_third)
    h_b = get_height_mm(PANEL_B, w_third)
    h_c = get_height_mm(PANEL_C, w_third)
    h_r1 = max(h_a, h_b, h_c)

    # Row 2: D + E + [blank]
    h_d = get_height_mm(PANEL_D, w_third)
    h_e = get_height_mm(PANEL_E, w_third)
    h_r2 = max(h_d, h_e)

    total_h = h_r1 + GAP_MM + h_r2 + 2
    print(f"\nFigure: {FIG_WIDTH_MM} x {total_h:.0f} mm")

    asm = VectorAssembler(FIG_WIDTH_MM, total_h)

    y = 0
    asm.place_panel('A', PANEL_A, x_mm=0, y_mm=y, w_mm=w_third, h_mm=h_r1)
    asm.place_panel('B', PANEL_B, x_mm=w_third + GAP_MM, y_mm=y, w_mm=w_third, h_mm=h_r1)
    asm.place_panel('C', PANEL_C, x_mm=2 * (w_third + GAP_MM), y_mm=y, w_mm=w_third, h_mm=h_r1)

    y += h_r1 + GAP_MM
    asm.place_panel('D', PANEL_D, x_mm=0, y_mm=y, w_mm=w_third, h_mm=h_r2)
    asm.place_panel('E', PANEL_E, x_mm=w_third + GAP_MM, y_mm=y, w_mm=w_third, h_mm=h_r2)

    asm.save('S5_Immune_Modules', output_dir=SCRIPT_DIR)
    print("Done!")


if __name__ == '__main__':
    assemble_figure()
