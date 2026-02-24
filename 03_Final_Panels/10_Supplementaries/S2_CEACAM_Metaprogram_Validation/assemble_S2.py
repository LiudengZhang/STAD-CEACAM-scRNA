#!/usr/bin/env python3
"""
Assemble Supplementary Figure S2: CEACAM & Metaprogram Validation

Layout:
  Row 1: A (Epithelial dotplot — full width)
  Row 2: B (Tumor Score UMAP, 50%) + C (CNV Score UMAP, 50%)
  Row 3: D (TCGA survival, 40%) + E (PRJEB coexpression, 30%) + F (TCGA coexpression, 30%)
  Row 4: G (MP gene signature heatmap — full width)
  Row 5: H (Liver CEACAM5/6 boxplot, 80mm) + I (Liver C2_CEACAM proportion, 45mm)

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
PANEL_A = SCRIPT_DIR / "S2_A" / "panel_S2_A.svg"
PANEL_B = SCRIPT_DIR / "S2_B" / "panel_S2_B.svg"
PANEL_C = SCRIPT_DIR / "S2_C" / "panel_S2_C.svg"
PANEL_D = SCRIPT_DIR / "S2_D" / "panel_S2_D.svg"
PANEL_E = SCRIPT_DIR / "S2_E" / "panel_S2_E.svg"
PANEL_F = SCRIPT_DIR / "S2_F" / "panel_S2_F.svg"
PANEL_G = SCRIPT_DIR / "S2_G" / "panel_S2_G.svg"
PANEL_H = SCRIPT_DIR / "S2_H" / "panel_S2_H.svg"
PANEL_I = SCRIPT_DIR / "S2_I" / "panel_S2_I.svg"


def get_height_mm(panel_path, width_mm):
    png = panel_path.with_suffix('.png')
    if png.exists():
        img = Image.open(png)
        return width_mm * (img.height / img.width)
    return width_mm * 0.5


def assemble_figure():
    print("Assembling Supplementary Figure S2: CEACAM & Metaprogram Validation")

    all_panels = [
        ('A', PANEL_A), ('B', PANEL_B), ('C', PANEL_C),
        ('D', PANEL_D), ('E', PANEL_E), ('F', PANEL_F),
        ('G', PANEL_G), ('H', PANEL_H), ('I', PANEL_I),
    ]
    for name, path in all_panels:
        svg_ok = path.exists()
        png_ok = path.with_suffix('.png').exists()
        print(f"  Panel {name}: SVG={'OK' if svg_ok else 'MISSING'}, PNG={'OK' if png_ok else 'MISSING'}")

    # Row 1: A — full width (dotplot)
    h_a = get_height_mm(PANEL_A, FIG_WIDTH_MM)

    # Row 2: B + C — two small UMAPs (45mm each), left-aligned
    w_bc = 45
    h_b = get_height_mm(PANEL_B, w_bc)
    h_c = get_height_mm(PANEL_C, w_bc)
    h_r2 = max(h_b, h_c)

    # Row 3: D (40%) + E (30%) + F (30%)
    usable_r3 = FIG_WIDTH_MM - 2 * GAP_MM
    w_d = usable_r3 * 0.40
    w_ef = usable_r3 * 0.30
    h_d = get_height_mm(PANEL_D, w_d)
    h_e = get_height_mm(PANEL_E, w_ef)
    h_f = get_height_mm(PANEL_F, w_ef)
    h_r3 = max(h_d, h_e, h_f)

    # Row 4: G — full width (heatmap)
    h_g = get_height_mm(PANEL_G, FIG_WIDTH_MM)

    # Row 5: H (80mm) + I (45mm)
    w_h = 80
    w_i = 45
    h_h = get_height_mm(PANEL_H, w_h)
    h_i = get_height_mm(PANEL_I, w_i)
    h_r5 = max(h_h, h_i)

    total_h = h_a + GAP_MM + h_r2 + GAP_MM + h_r3 + GAP_MM + h_g + GAP_MM + h_r5 + 2
    print(f"\nFigure: {FIG_WIDTH_MM} x {total_h:.0f} mm")

    asm = VectorAssembler(FIG_WIDTH_MM, total_h)

    y = 0
    # Row 1: A full width
    asm.place_panel('A', PANEL_A, x_mm=0, y_mm=y, w_mm=FIG_WIDTH_MM, h_mm=h_a)

    # Row 2: B + C (two UMAPs)
    y += h_a + GAP_MM
    asm.place_panel('B', PANEL_B, x_mm=0, y_mm=y, w_mm=w_bc, h_mm=h_r2)
    asm.place_panel('C', PANEL_C, x_mm=w_bc + GAP_MM, y_mm=y, w_mm=w_bc, h_mm=h_r2)

    # Row 3: D (40%) + E (30%) + F (30%)
    y += h_r2 + GAP_MM
    asm.place_panel('D', PANEL_D, x_mm=0, y_mm=y, w_mm=w_d, h_mm=h_r3)
    asm.place_panel('E', PANEL_E, x_mm=w_d + GAP_MM, y_mm=y, w_mm=w_ef, h_mm=h_r3)
    asm.place_panel('F', PANEL_F, x_mm=w_d + GAP_MM + w_ef + GAP_MM, y_mm=y, w_mm=w_ef, h_mm=h_r3)

    # Row 4: G full width
    y += h_r3 + GAP_MM
    asm.place_panel('G', PANEL_G, x_mm=0, y_mm=y, w_mm=FIG_WIDTH_MM, h_mm=h_g)

    # Row 5: H + I
    y += h_g + GAP_MM
    asm.place_panel('H', PANEL_H, x_mm=0, y_mm=y, w_mm=w_h, h_mm=h_r5)
    asm.place_panel('I', PANEL_I, x_mm=w_h + GAP_MM, y_mm=y, w_mm=w_i, h_mm=h_r5)

    asm.save('S2_CEACAM_Metaprogram_Validation', output_dir=SCRIPT_DIR)
    print("Done!")


if __name__ == '__main__':
    assemble_figure()
