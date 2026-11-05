#!/usr/bin/env python3
"""Supplementary Figure S1 A - the four QC-metric UMAPs (genes per cell, total
counts, mitochondrial %, doublet score).

Reads data/qc_cells.csv: a 30,000-cell subsample (numpy seed 42) of the full
dataset with UMAP x, y and the four metrics. The h5ad is opened only when that
table is absent or with --recompute.

Run: python create_S1_A_qc_umaps.py [--recompute]"""

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

FIG, PANEL = "S1_QC_Annotation", "S1_A"
N_CELLS = 30000
SEED = 42

# Panel layout in mm.
W, H = 171.0, 38.0
MAP_MM = 28.0                 # each map is square
LEFT_MM, TOP_MM = 5.0, 3.6    # room for "UMAP2" at the left, titles above
CBAR_W_MM, CBAR_GAP_MM = 1.2, 0.8
TICK_MM = 4.5                 # colour-bar tick labels ("20000")
COL_GAP_MM = 6.0              # the four cells spread over the row
DOT_MM = 0.18                 # a point's diameter on the page

PANELS = [
    ("n_genes", "Genes per Cell", None, None),
    ("total_counts", "Total Counts", None, None),
    ("pct_mt", "Mitochondrial %", 0, 25),
    ("doublet_score", "Doublet Score", 0, None),
]

def compute_cells():
    import scanpy as sc
    adata = sc.read_h5ad(FULL_DATASET_H5AD, backed="r")
    np.random.seed(SEED)
    n_plot = min(N_CELLS, adata.n_obs)
    idx = np.sort(np.random.choice(adata.n_obs, n_plot, replace=False))
    xy = np.asarray(adata.obsm["X_umap"][idx], dtype=float)
    obs = adata.obs.iloc[idx]
    return pd.DataFrame({
        "x": xy[:, 0], "y": xy[:, 1],
        "n_genes": obs["n_genes_by_counts"].to_numpy(dtype=float),
        "total_counts": obs["total_counts"].to_numpy(dtype=float),
        "pct_mt": obs["pct_counts_mt"].to_numpy(dtype=float),
        "doublet_score": obs["doublet_score"].to_numpy(dtype=float),
    })

def draw(cells):
    base.apply_style()
    fig = style.figure_mm(W, H)
    cell_w = MAP_MM + CBAR_GAP_MM + CBAR_W_MM + TICK_MM
    axes, caxes = [], []          # maps first, colour bars after
    for i in range(4):
        x = LEFT_MM + i * (cell_w + COL_GAP_MM)
        axes.append(fig.add_axes([x / W, 1 - (TOP_MM + MAP_MM) / H,
                                  MAP_MM / W, MAP_MM / H]))
    for i in range(4):
        x = LEFT_MM + i * (cell_w + COL_GAP_MM) + MAP_MM + CBAR_GAP_MM
        caxes.append(fig.add_axes([x / W, 1 - (TOP_MM + 0.85 * MAP_MM) / H,
                                   CBAR_W_MM / W, 0.7 * MAP_MM / H]))
    xy = cells[["x", "y"]].to_numpy()
    for ax, cax, (col, title, vmin, vmax) in zip(axes, caxes, PANELS):
        values = cells[col].to_numpy()
        order = np.argsort(values)           # high values drawn on top
        sc_ = ax.scatter(xy[order, 0], xy[order, 1], c=values[order],
                         cmap="viridis", s=(DOT_MM * style.PT_PER_MM) ** 2,
                         alpha=0.6, rasterized=True, vmin=vmin, vmax=vmax,
                         edgecolors="none", linewidths=0)
        ax.set_title(title, pad=2)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        cb = fig.colorbar(sc_, cax=cax)
        cb.ax.tick_params(labelsize=style.tick_pt(), width=style.RULE_PT,
                          length=1.5)
        cb.outline.set_linewidth(style.EDGE_PT)
    axes[0].set_xlabel("UMAP1", fontsize=style.tick_pt())
    axes[0].set_ylabel("UMAP2", fontsize=style.tick_pt())
    return fig

if __name__ == "__main__":
    base.save(draw(cache.table(HERE, "qc_cells", compute_cells)),
              FIG, PANEL, "panel_S1_A")
