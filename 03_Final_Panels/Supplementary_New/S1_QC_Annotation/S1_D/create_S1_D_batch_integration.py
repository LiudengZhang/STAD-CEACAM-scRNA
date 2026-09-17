#!/usr/bin/env python3
"""
S1 panel D - the 30,000-cell UMAP before and after Harmony, coloured by
sample, drawn at the size it prints at.

  printed panel  Supplementary Figure S1 D   (PROVENANCE.csv)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/S1_QC_Annotation/
                 S1_D/create_S1_D_batch_integration.py  (16 x 7 cm at 4x)

DRAWING READS A TABLE (cnsfig.cache). The predecessor subsamples 30,000 cells
(numpy seed 42), recomputes a UMAP on the uncorrected PCA of those cells
(scanpy neighbors k = 15, umap at its default random_state = 0) and reads
the corrected UMAP off the h5ad; data/umap_cells.csv holds both embeddings
and each cell's sample index in the Set2 cycle, and that recomputation runs
only when the table is absent or with --recompute. Not one cell moves.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import FULL_DATASET_H5AD                       # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S1_QC_Annotation", "S1_D"
N_CELLS = 30000
SEED = 42

# The printed box, millimetres: two maps beside E in one row. 40 mm tall
# since 2026-09-16 (S1 back on one page): the maps went from 37 to 31 mm.
W, H = 92.0, 40.0
MAP_MM = 31.0
LEFT_MM, TOP_MM, GAP_MM = 6.0, 3.6, 12.0
DOT_MM = 0.18


def compute_cells():
    """The computing half: the predecessor's subsample and its two UMAPs."""
    import scanpy as sc
    import anndata
    adata = sc.read_h5ad(FULL_DATASET_H5AD, backed="r")
    np.random.seed(SEED)
    n_plot = min(N_CELLS, adata.n_obs)
    idx = np.sort(np.random.choice(adata.n_obs, n_plot, replace=False))
    sample_ids = adata.obs["Sample ID"].iloc[idx].astype(str).to_numpy()
    pca = np.asarray(adata.obsm["X_pca_uncorrected"][idx])
    after = np.asarray(adata.obsm["X_umap"][idx], dtype=float)
    sub = anndata.AnnData(X=np.zeros((n_plot, 10)), obsm={"X_pca": pca})
    sub.obs["Sample ID"] = sample_ids
    sc.pp.neighbors(sub, use_rep="X_pca", n_neighbors=15)
    sc.tl.umap(sub)
    before = np.asarray(sub.obsm["X_umap"], dtype=float)
    samples = np.unique(sample_ids)
    to_idx = {s: i for i, s in enumerate(samples)}
    return pd.DataFrame({
        "before_x": before[:, 0], "before_y": before[:, 1],
        "after_x": after[:, 0], "after_y": after[:, 1],
        "sample_idx": [to_idx[s] for s in sample_ids],
    })


def draw(cells):
    import matplotlib.pyplot as plt
    base.apply_style()
    fig = style.figure_mm(W, H)
    cmap = plt.get_cmap("Set2", 8)
    colors = [cmap(int(i) % 8) for i in cells["sample_idx"]]
    for k, (prefix, title) in enumerate([("before", "Before Harmonization"),
                                         ("after", "After Harmonization")]):
        x = LEFT_MM + k * (MAP_MM + GAP_MM)
        ax = fig.add_axes([x / W, 1 - (TOP_MM + MAP_MM) / H, MAP_MM / W, MAP_MM / H])
        ax.scatter(cells[f"{prefix}_x"], cells[f"{prefix}_y"], c=colors,
                   s=(DOT_MM * style.PT_PER_MM) ** 2, alpha=0.6, rasterized=True,
                   edgecolors="none", linewidths=0)
        ax.set_title(title, pad=2)
        ax.set_xlabel("UMAP1", fontsize=style.tick_pt())
        ax.set_ylabel("UMAP2", fontsize=style.tick_pt())
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"ink outside the {W} x {H} mm canvas: {over}")
    return fig


def main():
    cells = cache.table(HERE, "umap_cells", compute_cells)
    base.save(draw(cells), FIG, PANEL, "panel_S1_D")
    return 0


if __name__ == "__main__":
    sys.exit(main())
