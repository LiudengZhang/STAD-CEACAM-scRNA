#!/usr/bin/env python3
"""
Assemble Supplementary Figure S1: QC & Annotation

Layout (7 rows):
  Row 1: A (4 QC UMAPs — full width)
  Row 2: B (Stomach QC boxplots, half height — full width)
  Row 3: C (Non-stomach QC boxplots, half height — full width)
  Row 4: D (Harmony UMAPs — 50%) + E (Doublet score — 50%)
  Row 5: F (Major cell type markers — full width)
  Row 6: G (T cell subtype markers — full width)
  Row 7: H (Non-stomach stacked bar — full width)

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
PANEL_A = SCRIPT_DIR / "S1_A" / "panel_S1_A.svg"
PANEL_B = SCRIPT_DIR / "S1_B" / "panel_S1_B.svg"
PANEL_C = SCRIPT_DIR / "S1_C" / "panel_S1_C.svg"
PANEL_D = SCRIPT_DIR / "S1_D" / "panel_S1_D.svg"
PANEL_E = SCRIPT_DIR / "S1_E" / "panel_S1_E.svg"
PANEL_F = SCRIPT_DIR / "S1_F" / "panel_S1_F.svg"
PANEL_G = SCRIPT_DIR / "S1_G" / "panel_S1_G.svg"
PANEL_H = SCRIPT_DIR / "S1_H" / "panel_S1_H.svg"


def get_height_mm(panel_path, width_mm):
    png = panel_path.with_suffix('.png')
    if png.exists():
        img = Image.open(png)
        return width_mm * (img.height / img.width)
    return width_mm * 0.5


def assemble_figure():
    print("Assembling Supplementary Figure S1: QC & Annotation")

    all_panels = [
        ('A', PANEL_A), ('B', PANEL_B), ('C', PANEL_C),
        ('D', PANEL_D), ('E', PANEL_E), ('F', PANEL_F),
        ('G', PANEL_G), ('H', PANEL_H),
    ]
    for name, path in all_panels:
        svg_ok = path.exists()
        png_ok = path.with_suffix('.png').exists()
        print(f"  Panel {name}: SVG={'OK' if svg_ok else 'MISSING'}, PNG={'OK' if png_ok else 'MISSING'}")

    w_half = (FIG_WIDTH_MM - GAP_MM) / 2

    # Row 1: A — 4 QC UMAPs (full width)
    h_a = get_height_mm(PANEL_A, FIG_WIDTH_MM)

    # Row 2: B — Stomach QC boxplots (full width, half height)
    h_b = get_height_mm(PANEL_B, FIG_WIDTH_MM)

    # Row 3: C — Non-stomach QC boxplots (full width, half height)
    h_c = get_height_mm(PANEL_C, FIG_WIDTH_MM)

    # Row 4: D (Harmony UMAPs) + E (Doublet score) — 50/50
    h_d = get_height_mm(PANEL_D, w_half)
    h_e = get_height_mm(PANEL_E, w_half)
    h_r4 = max(h_d, h_e)

    # Row 5: F — Major cell type markers (full width)
    h_f = get_height_mm(PANEL_F, FIG_WIDTH_MM)

    # Row 6: G — T cell subtype markers (full width)
    h_g = get_height_mm(PANEL_G, FIG_WIDTH_MM)

    # Row 7: H — Non-stomach stacked bar (full width)
    h_h = get_height_mm(PANEL_H, FIG_WIDTH_MM)

    total_h = h_a + GAP_MM + h_b + GAP_MM + h_c + GAP_MM + h_r4 + GAP_MM + h_f + GAP_MM + h_g + GAP_MM + h_h + 2
    print(f"\nFigure: {FIG_WIDTH_MM} x {total_h:.0f} mm")

    asm = VectorAssembler(FIG_WIDTH_MM, total_h)

    # Row 1: A (full width)
    y = 0
    asm.place_panel('A', PANEL_A, x_mm=0, y_mm=y, w_mm=FIG_WIDTH_MM, h_mm=h_a)

    # Row 2: B (full width)
    y += h_a + GAP_MM
    asm.place_panel('B', PANEL_B, x_mm=0, y_mm=y, w_mm=FIG_WIDTH_MM, h_mm=h_b)

    # Row 3: C (full width)
    y += h_b + GAP_MM
    asm.place_panel('C', PANEL_C, x_mm=0, y_mm=y, w_mm=FIG_WIDTH_MM, h_mm=h_c)

    # Row 4: D + E (50/50)
    y += h_c + GAP_MM
    asm.place_panel('D', PANEL_D, x_mm=0, y_mm=y, w_mm=w_half, h_mm=h_r4)
    asm.place_panel('E', PANEL_E, x_mm=w_half + GAP_MM, y_mm=y, w_mm=w_half, h_mm=h_r4)

    # Row 5: F (full width)
    y += h_r4 + GAP_MM
    asm.place_panel('F', PANEL_F, x_mm=0, y_mm=y, w_mm=FIG_WIDTH_MM, h_mm=h_f)

    # Row 6: G (full width)
    y += h_f + GAP_MM
    asm.place_panel('G', PANEL_G, x_mm=0, y_mm=y, w_mm=FIG_WIDTH_MM, h_mm=h_g)

    # Row 7: H (full width)
    y += h_g + GAP_MM
    asm.place_panel('H', PANEL_H, x_mm=0, y_mm=y, w_mm=FIG_WIDTH_MM, h_mm=h_h)

    asm.save('S1_QC_Annotation', output_dir=SCRIPT_DIR)
    print("Done!")


if __name__ == '__main__':
    assemble_figure()
