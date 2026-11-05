#!/usr/bin/env python3
"""
Figure 1B: UMAP of the stomach samples coloured by major cell type, with the
cell-type names on the embedding.

Reads data/umap_cells.csv (x, y, cell_type per stomach cell); the table is
computed from the full h5ad when absent or with --recompute.
Writes umap_major_cell_type_stomach.[svg|pdf|png] beside this script.

Run: python create_stomach_umap.py [--recompute]
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FULL_DATASET_H5AD                       # noqa: E402
from cnsfig import cache, Panel, style                    # noqa: E402
from cnsfig.layout import umap_corner_axis                 # noqa: E402

PANEL = Panel(1, "B", __file__, small_pt=6.0)
# Square plotting box: below a 4.5 mm letter band, with a 1 mm foot.
CLOUD_MM = min(PANEL.w_mm - 2.0, PANEL.h_mm - 4.5 - 1.0)
_SIDE_MM = (PANEL.w_mm - CLOUD_MM) / 2.0
MARGIN = dict(left=_SIDE_MM, right=_SIDE_MM, top=4.5,
              bottom=PANEL.h_mm - 4.5 - CLOUD_MM)

# Point diameter on the page; small enough to avoid overplotting ~100k cells.
DOT_MM = 0.18

# Major-cell-type palette shared with Figure 1C.
from shared.figure_config import CELL_TYPE_COLORS  # noqa: E402,F401

MERGE = {
    'CD4+ T cells': 'T/NK cells',
    'CD8+ T cells': 'T/NK cells',
    'NK cells': 'T/NK cells',
    'DC cells': 'Dendritic cells',
    'Fibroblast': 'Fibroblasts',
    'Pericyte': 'Pericytes',
    'Hepatocyte': 'Hepatocytes',
}

# The data's own names are drawn, as in 1C's key.
DISPLAY_NAMES = {}

# Hand-placed label positions (UMAP units) for scattered or touching
# populations; any other name is set at its population's median.
PINNED = {
    'Epithelial cells': (7.6, 5.0),
    'Hepatocytes': (7.3, 10.6),
    'Neutrophils': (1.9, 8.6),
    'Pericytes': (9.0, 20.3),
    'Fibroblasts': (12.2, 16.5),
    'Dendritic cells': (-0.5, 16.3),
    'Monocytes/Macrophages': (-2.2, 12.4),
    'Plasma cells': (8.8, 13.55),
    'Endothelial cells': (13.9, 14.7),
}
# Axis ranges: the data's extent plus room for the names at the right edge.
XLIM = (-12.0, 18.0)
YLIM = (-0.7, 22.4)


def compute_cells():
    """Every stomach cell's UMAP position and type."""
    import scanpy as sc
    adata = sc.read_h5ad(FULL_DATASET_H5AD)
    st = adata[adata.obs['Sample site'] == 'Stomach']
    if 'X_umap' not in st.obsm:
        raise ValueError("UMAP coordinates not found in adata.obsm['X_umap']")
    cell_type = st.obs['major_cell_type'].astype(str).replace(MERGE)
    # Five decimals is far below print resolution and keeps the table small.
    xy = np.round(np.asarray(st.obsm['X_umap'], dtype=float), 5)
    return pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1],
                         "cell_type": cell_type.to_numpy()})


def main():
    cells = cache.table(PANEL.dir, "umap_cells", compute_cells)

    fig, ax = PANEL.subplots()
    style.margins_mm(fig, **MARGIN)

    colors = cells['cell_type'].map(CELL_TYPE_COLORS).fillna('#999999').to_numpy()
    sc_ = ax.scatter(cells['x'], cells['y'], c=colors,
                     s=(DOT_MM * style.PT_PER_MM) ** 2, alpha=0.6,
                     edgecolors='none', linewidths=0, rasterized=True)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    for ct, sub in cells.groupby('cell_type'):
        cx, cy = PINNED.get(ct, (float(sub['x'].median()), float(sub['y'].median())))
        ax.text(cx, cy,
                DISPLAY_NAMES.get(ct, ct), **style.umap_label_kw())
    # White under-layer so the names read on top of the cloud.
    style.umap_label_underlays(ax)

    # UMAP1/UMAP2 arrows at the lower left.
    umap_corner_axis(ax)

    PANEL.finish(fig, 'umap_major_cell_type_stomach')


if __name__ == '__main__':
    main()
