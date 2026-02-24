#!/usr/bin/env python3
"""
Assemble Supplementary Figure S6: CD274 Remaining Cell Types

Layout: 2 rows — 5 + 4 panels (DC moved to main figure)
  Row 1: A (CD8+T) + B (CD4+T) + C (NK) + D (B cells) + E (Plasma)
  Row 2: F (Mast) + G (Neutrophils) + H (Endothelial) + I (Pericytes)

Output: Vector SVG + PDF + PNG preview.
"""

import sys
from pathlib import Path
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from shared.svg_assembler import VectorAssembler

FIG_WIDTH_MM = 183
GAP_MM = 2

SCRIPT_DIR = Path(__file__).parent
PANELS = {}
for letter in 'ABCDEFGHI':
    PANELS[letter] = SCRIPT_DIR / f"S6_{letter}" / f"panel_S6_{letter}.svg"


def get_height_mm(panel_path, width_mm):
    png = panel_path.with_suffix('.png')
    if png.exists():
        img = Image.open(png)
        return width_mm * (img.height / img.width)
    return width_mm * 0.8


def assemble_figure():
    print("Assembling Supplementary Figure S6: CD274 Remaining Cell Types")

    for name, path in PANELS.items():
        svg_ok = path.exists()
        png_ok = path.with_suffix('.png').exists()
        print(f"  Panel {name}: SVG={'OK' if svg_ok else 'MISSING'}, PNG={'OK' if png_ok else 'MISSING'}")

    # 5 panels per row
    w_panel = (FIG_WIDTH_MM - 4 * GAP_MM) / 5

    # Row 1: A-E
    row1_letters = list('ABCDE')
    row1_heights = [get_height_mm(PANELS[l], w_panel) for l in row1_letters]
    h_r1 = max(row1_heights)

    # Row 2: F-I (4 panels)
    row2_letters = list('FGHI')
    row2_heights = [get_height_mm(PANELS[l], w_panel) for l in row2_letters]
    h_r2 = max(row2_heights)

    total_h = h_r1 + GAP_MM + h_r2 + 2
    print(f"\nFigure: {FIG_WIDTH_MM} x {total_h:.0f} mm")

    asm = VectorAssembler(FIG_WIDTH_MM, total_h)

    # Row 1
    y = 0
    for i, letter in enumerate(row1_letters):
        x = i * (w_panel + GAP_MM)
        asm.place_panel(letter, PANELS[letter], x_mm=x, y_mm=y, w_mm=w_panel, h_mm=h_r1)

    # Row 2
    y += h_r1 + GAP_MM
    for i, letter in enumerate(row2_letters):
        x = i * (w_panel + GAP_MM)
        asm.place_panel(letter, PANELS[letter], x_mm=x, y_mm=y, w_mm=w_panel, h_mm=h_r2)

    asm.save('S6_CD274_Remaining', output_dir=SCRIPT_DIR)
    print("Done!")


if __name__ == '__main__':
    assemble_figure()
