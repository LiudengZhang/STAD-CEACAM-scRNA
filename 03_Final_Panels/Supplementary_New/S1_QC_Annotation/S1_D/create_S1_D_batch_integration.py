#!/usr/bin/env python3
"""Supplementary Figure S1 D - the 30,000-cell UMAP before and after Harmony,
coloured by sample.

Reads data/umap_cells.csv (both embeddings and each cell's sample index). With
--recompute, or when that table is absent, it is rebuilt from the full dataset:
30,000 cells at numpy seed 42, a UMAP on the uncorrected PCA (k = 15, scanpy's
default random_state), and the stored Harmony UMAP.

Run: python create_S1_D_batch_integration.py [--recompute]"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import FULL_DATASET_H5AD                       # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S1_QC_Annotation", "S1_D"
N_CELLS = 30000
SEED = 42

# Panel layout in mm.
W, H = 92.0, 40.0
MAP_MM = 31.0
LEFT_MM, TOP_MM, GAP_MM = 6.0, 3.6, 12.0
DOT_MM = 0.18

def compute_cells():
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
    return fig

if __name__ == "__main__":
    base.save(draw(cache.table(HERE, "umap_cells", compute_cells)),
              FIG, PANEL, "panel_S1_D")
