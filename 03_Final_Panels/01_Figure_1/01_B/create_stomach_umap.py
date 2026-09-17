#!/usr/bin/env python3
"""
Figure 1 panel B - UMAP of the stomach samples coloured by major cell type,
with the cell-type names on the embedding.

  printed panel  Figure 1 B       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

DRAWN AT THE SIZE IT PRINTS, since 2026-09-15. Until then this panel was made
on a 12 x 10 inch canvas, cropped to its ink and scaled into a 254 mm
landscape page, which the journal prints at column width: its 8 pt labels
came out near 4 pt. Figure 1 is now assembled like Figures 2-5, into
measured slots on the 171.10 mm page (03_Final_Panels/panel_rects_v2.csv),
and this drawing is made at that slot's millimetres with the type set by
00_Config/panel_style_cns.py.

BESIDE C, AT C'S HEIGHT, since 2026-09-16 (the author's fifth reading: "1B
and 1C should be the same height"): the slot is 63 x 62 mm at the left of the
second row (build_grid_v2.REPAGED), the cloud a square that fills the slot -
its side the smaller of the width less the side margins and the height less
the letter band and a 1 mm foot - centred under the letter band at equal
aspect (adjustable='datalim' expands the y range about its centre). The names
print without plates. The history: for a day (2026-09-15, afternoon) B had a
full-width row and a 96 mm cloud; that evening B and C went on one row with
C's bars turned across the panel and the row 96 mm tall, so B's 56 mm cloud
hung from the top of a 96 mm slot over 35 mm of paper; now C's key stands
under its upright bars and both panels are 62 mm.

DRAWING READS A TABLE (cnsfig.cache). The embedding lives in
data/umap_cells.csv - one row per stomach cell: x, y, cell_type - written
from the full-dataset h5ad only when the table is absent or --recompute is
passed. The merge of T/NK, the renamings and the palette are the earlier
script's; not one cell moves. The specimen number is not written.

The points are rasterised inside the vector file, as Figure 2A's are: the
type, the labels and the frame stay vectors.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FULL_DATASET_H5AD                       # noqa: E402
import panel_style_cns as style                           # noqa: E402
import slots                                              # noqa: E402
from cnsfig import cache                                  # noqa: E402

PANEL_LETTER = "B"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(1, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(1, PANEL_LETTER)
OUTPUT_DIR = Path(__file__).resolve().parent
#: The cloud is square and fills the slot: its side is the smaller of the
#: width less the side margins and the height less the letter band (4.5 mm)
#: and a 1 mm foot; it is centred across the slot under the letter band.
CLOUD_MM = min(PANEL_W_MM - 2.0, PANEL_H_MM - 4.5 - 1.0)
_SIDE_MM = (PANEL_W_MM - CLOUD_MM) / 2.0
MARGIN = dict(left=_SIDE_MM, right=_SIDE_MM, top=4.5,
              bottom=PANEL_H_MM - 4.5 - CLOUD_MM)

#: A point's diameter on the page. The published panel's cells are a fine
#: grain; 0.18 mm keeps the twelve populations as fields of colour without
#: the overplotting a larger dot gives on 100,000 cells.
DOT_MM = 0.18

# Set2 + Set3 palette (12 major cell types, consistent across Fig 1B & 1C)
CELL_TYPE_COLORS = {
    'Epithelial cells': '#e5c494',
    'T/NK cells': '#66c2a5',
    'Monocytes/Macrophages': '#fc8d62',
    'Plasma cells': '#8da0cb',
    'B cells': '#ffd92f',
    'Endothelial cells': '#a6d854',
    'Fibroblasts': '#e78ac3',
    'Neutrophils': '#b3b3b3',
    'Mast cells': '#fb8072',
    'Pericytes': '#bebada',
    'Dendritic cells': '#80b1d3',
    'Hepatocytes': '#bc80bd',
}

MERGE = {
    'CD4+ T cells': 'T/NK cells',
    'CD8+ T cells': 'T/NK cells',
    'NK cells': 'T/NK cells',
    'DC cells': 'Dendritic cells',
    'Fibroblast': 'Fibroblasts',
    'Pericyte': 'Pericytes',
    'Hepatocyte': 'Hepatocytes',
}

# The on-plot names are the published page's (00_GROUND_TRUTH/figures/
# Figure 1.pdf, twelve Arial-Bold strings), declared in labels.py
# RENAMES_FIGURE_1 against the shorter forms the earlier script set.
DISPLAY_NAMES = {
    'Epithelial cells': 'Epithelial Cells',
    'T/NK cells': 'T/NK Cells',
    'Monocytes/Macrophages': 'Monocytes/Macrophages',
    'Plasma cells': 'Plasma Cells',
    'B cells': 'B Cells',
    'Endothelial cells': 'Endothelial Cells',
    'Fibroblasts': 'Fibroblasts',
    'Neutrophils': 'Neutrophils',
    'Mast cells': 'Mast Cells',
    'Pericytes': 'Pericytes',
    'Dendritic cells': 'Dendritic Cells',
    'Hepatocytes': 'Hepatocytes',
}

#: Where each name is set, in UMAP units, where the population's median is
#: not the place for it - as the published page places them by hand:
#: Hepatocytes (334 cells scattered through the epithelium) above the
#: epithelial mass, Pericytes above and Fibroblasts below-right of their
#: two touching clusters, Dendritic Cells above Monocytes/Macrophages,
#: Endothelial Cells above their cluster at the right edge. A name with no
#: entry is set at its median.
PINNED = {
    'Epithelial cells': (7.6, 5.0),
    'Hepatocytes': (7.3, 10.6),
    'Neutrophils': (1.9, 8.6),
    'Pericytes': (9.0, 20.3),
    'Fibroblasts': (12.2, 16.5),
    'Dendritic cells': (-0.5, 16.3),
    'Monocytes/Macrophages': (-2.2, 12.4),
    'Plasma cells': (8.8, 13.4),
    'Endothelial cells': (13.9, 14.7),
}
#: The x range, in UMAP units: the embedding spans -11.4 to 15.4 and the
#: names at its right edge need paper beyond it.
XLIM = (-12.5, 19.0)


def compute_cells():
    """The computing half: every stomach cell's UMAP position and type."""
    import scanpy as sc
    print(f"Loading data from: {FULL_DATASET_H5AD}")
    adata = sc.read_h5ad(FULL_DATASET_H5AD)
    print(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    st = adata[adata.obs['Sample site'] == 'Stomach']
    print(f"Stomach samples: {st.shape[0]:,} cells, "
          f"{st.obs['sample'].nunique()} samples")
    if 'X_umap' not in st.obsm:
        raise ValueError("UMAP coordinates not found in adata.obsm['X_umap']")
    cell_type = st.obs['major_cell_type'].astype(str).replace(MERGE)
    # Five decimals: 1e-5 UMAP units is 2e-5 mm on this panel, and the
    # float32 coordinates written in full made the table 13 MB.
    xy = np.round(np.asarray(st.obsm['X_umap'], dtype=float), 5)
    return pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1],
                         "cell_type": cell_type.to_numpy()})


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt")

    cells = cache.table(OUTPUT_DIR, "umap_cells", compute_cells)
    counts = cells['cell_type'].value_counts()
    print("Cell type distribution (stomach only):")
    for ct, n in counts.items():
        print(f"  {ct:25s}: {n:7,d} cells ({100 * n / len(cells):5.2f}%)")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    style.margins_mm(fig, **MARGIN)

    colors = cells['cell_type'].map(CELL_TYPE_COLORS).fillna('#999999').to_numpy()
    sc_ = ax.scatter(cells['x'], cells['y'], c=colors,
                     s=(DOT_MM * style.PT_PER_MM) ** 2, alpha=0.6,
                     edgecolors='none', linewidths=0, rasterized=True)
    ax.set_xlim(*XLIM)
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    # The names at each population's median, bold, straight on the points:
    # the white plates the earlier drawing set behind them went on
    # 2026-09-15 (the author's fourth reading). No halo either - a stroked
    # halo turns the glyphs into paths and the words stop being text.
    for ct, sub in cells.groupby('cell_type'):
        cx, cy = PINNED.get(ct, (float(sub['x'].median()), float(sub['y'].median())))
        ax.text(cx, cy,
                DISPLAY_NAMES.get(ct, ct), fontsize=style.tick_pt(),
                fontweight='bold', ha='center', va='center', zorder=10)

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, OUTPUT_DIR / 'umap_major_cell_type_stomach')
    print(f"Saved: umap_major_cell_type_stomach.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
