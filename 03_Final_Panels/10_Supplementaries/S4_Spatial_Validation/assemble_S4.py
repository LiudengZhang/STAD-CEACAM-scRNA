#!/usr/bin/env python3
"""
Assemble Supplementary Figure S4: Spatial Validation

Layout (3 rows x 2 columns):
  Row 1: A (CEACAM ratio gallery) + B (MoMac fraction gallery)
  Row 2: C (Epithelial density gallery) + D (Distance to stroma)
  Row 3: E (Distance to immune) + F (Stacked bar)

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
PANEL_A = SCRIPT_DIR / "S4_A" / "panel_S4_A.svg"
PANEL_B = SCRIPT_DIR / "S4_B" / "panel_S4_B.svg"
PANEL_C = SCRIPT_DIR / "S4_C" / "panel_S4_C.svg"
PANEL_D = SCRIPT_DIR / "S4_D" / "panel_S4_D.svg"
PANEL_E = SCRIPT_DIR / "S4_E" / "panel_S4_E.svg"
PANEL_F = SCRIPT_DIR / "S4_F" / "panel_S4_F.svg"


def get_height_mm(panel_path, width_mm):
    png = panel_path.with_suffix('.png')
    if png.exists():
        img = Image.open(png)
        return width_mm * (img.height / img.width)
    return width_mm * 0.5


def assemble_figure():
    print("Assembling Supplementary Figure S4: Spatial Validation")

    all_panels = [
        ('A', PANEL_A), ('B', PANEL_B), ('C', PANEL_C),
        ('D', PANEL_D), ('E', PANEL_E), ('F', PANEL_F),
    ]
    for name, path in all_panels:
        svg_ok = path.exists()
        png_ok = path.with_suffix('.png').exists()
        print(f"  Panel {name}: SVG={'OK' if svg_ok else 'MISSING'}, PNG={'OK' if png_ok else 'MISSING'}")

    # 3 rows x 2 columns — equal halves
    w_half = (FIG_WIDTH_MM - GAP_MM) / 2

    # Row 1: A + B
    h_a = get_height_mm(PANEL_A, w_half)
    h_b = get_height_mm(PANEL_B, w_half)
    h_r1 = max(h_a, h_b)

    # Row 2: C + D
    h_c = get_height_mm(PANEL_C, w_half)
    h_d = get_height_mm(PANEL_D, w_half)
    h_r2 = max(h_c, h_d)

    # Row 3: E + F
    h_e = get_height_mm(PANEL_E, w_half)
    h_f = get_height_mm(PANEL_F, w_half)
    h_r3 = max(h_e, h_f)

    total_h = h_r1 + GAP_MM + h_r2 + GAP_MM + h_r3 + 2
    print(f"\nFigure: {FIG_WIDTH_MM} x {total_h:.0f} mm")

    asm = VectorAssembler(FIG_WIDTH_MM, total_h)

    # Row 1: A + B
    y = 0
    asm.place_panel('A', PANEL_A, x_mm=0, y_mm=y, w_mm=w_half, h_mm=h_r1)
    asm.place_panel('B', PANEL_B, x_mm=w_half + GAP_MM, y_mm=y, w_mm=w_half, h_mm=h_r1)

    # Row 2: C + D
    y += h_r1 + GAP_MM
    asm.place_panel('C', PANEL_C, x_mm=0, y_mm=y, w_mm=w_half, h_mm=h_r2)
    asm.place_panel('D', PANEL_D, x_mm=w_half + GAP_MM, y_mm=y, w_mm=w_half, h_mm=h_r2)

    # Row 3: E + F
    y += h_r2 + GAP_MM
    asm.place_panel('E', PANEL_E, x_mm=0, y_mm=y, w_mm=w_half, h_mm=h_r3)
    asm.place_panel('F', PANEL_F, x_mm=w_half + GAP_MM, y_mm=y, w_mm=w_half, h_mm=h_r3)

    asm.save('S4_Spatial_Validation', output_dir=SCRIPT_DIR)
    print("Done!")


if __name__ == '__main__':
    assemble_figure()
