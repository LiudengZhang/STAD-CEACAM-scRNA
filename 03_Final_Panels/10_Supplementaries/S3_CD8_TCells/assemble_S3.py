#!/usr/bin/env python3
"""
Assemble Supplementary Figure S3: CD8+ T Cell States

Layout:
  Row 1: A (CD8 marker dotplot — full width)
  Row 2: B (PDCD1 UMAP) + C (HAVCR2 UMAP) + D (Tex score UMAP)

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
PANEL_A = SCRIPT_DIR / "S3_A" / "panel_S3_A.svg"
PANEL_B = SCRIPT_DIR / "S3_B" / "panel_S3_B.svg"
PANEL_C = SCRIPT_DIR / "S3_C" / "panel_S3_C.svg"
PANEL_D = SCRIPT_DIR / "S3_D" / "panel_S3_D.svg"


def get_height_mm(panel_path, width_mm):
    png = panel_path.with_suffix('.png')
    if png.exists():
        img = Image.open(png)
        return width_mm * (img.height / img.width)
    return width_mm * 0.5


def assemble_figure():
    print("Assembling Supplementary Figure S3: CD8+ T Cell States")

    all_panels = [('A', PANEL_A), ('B', PANEL_B), ('C', PANEL_C), ('D', PANEL_D)]
    for name, path in all_panels:
        svg_ok = path.exists()
        png_ok = path.with_suffix('.png').exists()
        print(f"  Panel {name}: SVG={'OK' if svg_ok else 'MISSING'}, PNG={'OK' if png_ok else 'MISSING'}")

    # Row 1: A full width
    h_a = get_height_mm(PANEL_A, FIG_WIDTH_MM)

    # Row 2: B, C, D equal thirds
    w_panel = (FIG_WIDTH_MM - 2 * GAP_MM) / 3
    h_b = get_height_mm(PANEL_B, w_panel)
    h_c = get_height_mm(PANEL_C, w_panel)
    h_d = get_height_mm(PANEL_D, w_panel)
    h_r2 = max(h_b, h_c, h_d)

    total_h = h_a + GAP_MM + h_r2 + 2
    print(f"\nFigure: {FIG_WIDTH_MM} x {total_h:.0f} mm")

    asm = VectorAssembler(FIG_WIDTH_MM, total_h)

    y = 0
    asm.place_panel('A', PANEL_A, x_mm=0, y_mm=y, w_mm=FIG_WIDTH_MM, h_mm=h_a)

    y += h_a + GAP_MM
    asm.place_panel('B', PANEL_B, x_mm=0, y_mm=y, w_mm=w_panel, h_mm=h_r2)
    asm.place_panel('C', PANEL_C, x_mm=w_panel + GAP_MM, y_mm=y, w_mm=w_panel, h_mm=h_r2)
    asm.place_panel('D', PANEL_D, x_mm=2 * (w_panel + GAP_MM), y_mm=y, w_mm=w_panel, h_mm=h_r2)

    asm.save('S3_CD8_TCells', output_dir=SCRIPT_DIR)
    print("Done!")


if __name__ == '__main__':
    assemble_figure()
