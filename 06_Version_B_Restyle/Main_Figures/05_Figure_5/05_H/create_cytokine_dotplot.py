#!/usr/bin/env python3
"""
Figure 5 panel F, RESTYLED (Version B) - dotplot of TNF, IL6, IL1B and IL1A
across the thirteen major cell types.

Version A is
`03_Revised_Panels/Main_Figures/05_Figure_5/05_H/create_cytokine_dotplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows.

The gene list, the cell-type order, the `major_cell_type` filter (Hepatocyte
excluded), `use_raw`, `cmap='Reds'` and `standard_scale='var'` are Version A's,
unchanged. `sc.pl.dotplot` still computes and draws everything.

  printed panel  Figure 5 F        (PROVENANCE.csv - the directory is "05_H";
                 do NOT read the directory as the letter)
  printed rect   38.5 x 32.8 mm    (panel_rects.csv)
  Version B box  95.0 x 80.0 mm

    Thirteen cell-type names are set as y tick labels, and the longest,
    "Monocytes/Macrophages", is about 30 mm at 7 pt. The label column alone is
    therefore wider than the whole published panel; thirteen rows need about
    5 mm each on top of that.

MARK
    This is the one Figure 5 panel with no `SCALE` and no `N * SCALE` anywhere:
    Version A set only the font family and the font types, so every string it
    drew came out at matplotlib's default 10 pt on a canvas scanpy sized
    itself (75.7 x 141.0 mm here). Read against PANEL_SPEC that is SCALE = 1
    and SMALL_PT = 10, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 10 = 0.7

    and it is unused: the script specifies no marker area, line width or cap
    size of its own - scanpy sizes the dots from the axes it is given. The
    factor is stated anyway so the panel is auditable like the rest.

Three deviations from the usual recipe, forced by scanpy drawing the figure:

  1. `sc.pl.dotplot` builds its own figure, so the panel box is passed to it as
     `figsize=` rather than through `style.subplots_mm`. scanpy's `figsize` is
     the whole figure, measured here: asking for 95 x 80 mm gives a 95 x 80 mm
     canvas.
  2. Version A's `plt.figure(figsize=...)` created an empty figure that
     scanpy then ignored - `sc.pl.dotplot` always makes a new one, so the two
     computed sizes were never used for anything and nothing was ever drawn on
     that canvas. It is dropped rather than left to make a second, blank
     figure; nothing drawn changes.
  3. scanpy lays the dotplot out flush against the left edge and lets the y
     tick labels hang outside the canvas - measured at 15.4 mm of overflow
     before any adjustment. `style.margins_mm` moves scanpy's gridspec inwards
     to make room; `style.overflow_mm` then returns all zeros. Version A's
     `plt.tight_layout()` is what did this job when the save cropped to the
     ink, and a 1:1 save does not crop.
"""

import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import FULL_DATASET_H5AD
import panel_style_cns as style

OUTPUT_DIR = Path(__file__).resolve().parent

SCALE = 1                       # Version A had no SCALE: it drew 1:1 at 10 pt
SMALL_PT = 10.0                 # matplotlib's default font.size, unchanged by A

PRINTED_MM = (38.5, 32.8)       # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 95.0, 80.0
MARGIN = dict(left=34.0, right=2.0, top=2.0, bottom=8.0)

# Genes of interest
GENES = ['TNF', 'IL6', 'IL1B', 'IL1A']

# Cell type order (matching 05_F radar, using original adata labels)
CELL_TYPE_ORDER = [
    'B cells',
    'DC cells',
    'Endothelial cells',
    'Epithelial cells',
    'Fibroblast',
    'Mast cells',
    'Monocytes/Macrophages',
    'Neutrophils',
    'NK cells',
    'Pericyte',
    'Plasma cells',
    'CD4+ T cells',
    'CD8+ T cells',
]


def main():
    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f} (unused here)")

    print("=" * 60)
    print("Panel H: Cytokine Dotplot (TNF, IL6, IL1B, IL1A)")
    print("=" * 60)

    print("\n[1/3] Loading data...")
    adata = sc.read_h5ad(FULL_DATASET_H5AD)
    print(f"  {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # Filter to 13 cell types (exclude Hepatocyte)
    mask = adata.obs['major_cell_type'].isin(CELL_TYPE_ORDER)
    adata = adata[mask].copy()
    print(f"  After filtering: {adata.shape[0]:,} cells, {adata.obs['major_cell_type'].nunique()} types")

    # Validate genes
    var_names = adata.raw.var_names if adata.raw is not None else adata.var_names
    available = [g for g in GENES if g in var_names]
    missing = [g for g in GENES if g not in var_names]
    print(f"  Genes available: {available}")
    if missing:
        print(f"  WARNING — missing: {missing}")

    print("\n[2/3] Creating dotplot...")
    sc.pl.dotplot(
        adata,
        var_names=available,
        groupby='major_cell_type',
        categories_order=CELL_TYPE_ORDER,
        use_raw=(adata.raw is not None),
        cmap='Reds',
        show=False,
        save=None,
        standard_scale='var',
        figsize=style.figsize_mm(PANEL_W_MM, PANEL_H_MM),
    )

    fig = plt.gcf()
    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    print("\n[3/3] Saving...")
    style.save_panel(fig, OUTPUT_DIR / 'cytokine_dotplot')
    print(f"  Saved: cytokine_dotplot.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")
    print("Done!")


if __name__ == "__main__":
    main()
