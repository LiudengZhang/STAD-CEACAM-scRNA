#!/usr/bin/env python3
"""
Figure 2D: CEACAM5 versus CEACAM6 in pre-treatment stomach epithelium, as
kNN metacells of 10 cells built within each sample.

Reads the epithelial h5ad (EPITHELIAL_H5AD); expression from .raw, or .X when
the file has no .raw; metacells from the first 30 PCs of X_pca.
Writes ceacam_metacell_correlation.[svg|pdf|png], ceacam_metacell_values.csv
and ceacam_correlation_levels.csv (Spearman rho per cell, metacell and sample)
beside this script.

Run: python create_ceacam_metacell_correlation.py

Per-cell correlation is limited by dropout, so the panel shows metacells.
Nearest neighbours are exact, so the result is deterministic.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD  # noqa: E402
from shared.expression import assert_x_stands_in_for_raw  # noqa: E402
from cnsfig import Panel, style  # noqa: E402


SCALE = 1                           # sets the point size relative to the tick size
SMALL_PT = 4.6

PANEL = Panel(2, "D", __file__, small_pt=SMALL_PT)

K = 10
N_PCS = 30


def load():
    """Pre-treatment stomach epithelium, expression from .raw."""
    ad = sc.read_h5ad(EPITHELIAL_H5AD)
    # Raises if this file has no .raw and its .X is not known to equal it.
    assert_x_stands_in_for_raw(ad, EPITHELIAL_H5AD)
    ad = ad[(ad.obs["Sample site"] == "Stomach")
            & (ad.obs["Treatment phase"] == "Pre")].copy()
    # The deposited h5ad has no .raw; its .X holds the same log1p matrix.
    src = ad.raw if ad.raw is not None else ad
    out = {}
    for g in ("CEACAM5", "CEACAM6"):
        i = list(src.var_names).index(g)
        x = src.X[:, i]
        out[g] = x.toarray().flatten() if hasattr(x, "toarray") \
            else np.asarray(x).flatten()
    df = pd.DataFrame(out)
    df["sample"] = ad.obs["sample"].astype(str).values
    if "X_pca" not in ad.obsm:
        raise SystemExit("Epithelial.h5ad has no X_pca - cannot build metacells")
    return df, np.asarray(ad.obsm["X_pca"])[:, :N_PCS]


def pools_knn(emb, k):
    """A disjoint cover of the sample by kNN pools of exactly k cells."""
    nn = NearestNeighbors(n_neighbors=k).fit(emb)
    _, nbr = nn.kneighbors(emb)
    used, pools = np.zeros(len(emb), bool), []
    for row in nbr:
        if used[row[0]]:
            continue
        take = row[~used[row]]
        if len(take) < k:
            continue
        used[take] = True
        pools.append(take)
    return pools


def main():
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2

    df, emb = load()

    rho_cell, p_cell = stats.spearmanr(df["CEACAM5"], df["CEACAM6"])

    rows = []
    samples = df["sample"].values
    c5, c6 = df["CEACAM5"].values, df["CEACAM6"].values
    for s in sorted(set(samples)):
        idx = np.where(samples == s)[0]
        if len(idx) < K:
            continue
        for p in pools_knn(emb[idx], K):
            rows.append(dict(sample=s, CEACAM5=float(c5[idx[p]].mean()),
                             CEACAM6=float(c6[idx[p]].mean())))
    mc = pd.DataFrame(rows)
    if mc.empty:
        raise SystemExit("no metacells were formed")
    rho_mc, p_mc = stats.spearmanr(mc["CEACAM5"], mc["CEACAM6"])

    agg = df.groupby("sample")[["CEACAM5", "CEACAM6"]].mean()
    rho_s, p_s = stats.spearmanr(agg["CEACAM5"], agg["CEACAM6"])

    mc.to_csv(PANEL.dir / "ceacam_metacell_values.csv", index=False)
    pd.DataFrame([
        dict(level="cell", n=len(df), rho=float(rho_cell), p=float(p_cell)),
        dict(level=f"metacell_k{K}", n=len(mc), rho=float(rho_mc), p=float(p_mc)),
        dict(level="sample", n=len(agg), rho=float(rho_s), p=float(p_s)),
    ]).to_csv(PANEL.dir / "ceacam_correlation_levels.csv", index=False)

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)
    ax.scatter(mc["CEACAM5"], mc["CEACAM6"], s=1.4 * AREA, alpha=0.45,
               c="#3498db", edgecolors="none")
    ax.set_xlabel('CEACAM5', fontstyle='italic')
    ax.set_ylabel('CEACAM6', fontstyle='italic')

    hi = float(np.ceil(max(mc["CEACAM5"].max(), mc["CEACAM6"].max())))
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)
    ax.set_xticks(np.arange(0, hi + 1, 2))
    ax.set_yticks(np.arange(0, hi + 1, 2))
    ax.set_box_aspect(1)

    ax.text(0.04, 0.97, f"ρ = {rho_mc:.2f}",
            transform=ax.transAxes, fontsize=style.tick_pt(), va="top",
            linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85,
                      edgecolor="none"))

    PANEL.fit(fig)
    PANEL.finish(fig, 'ceacam_metacell_correlation')


if __name__ == "__main__":
    main()
