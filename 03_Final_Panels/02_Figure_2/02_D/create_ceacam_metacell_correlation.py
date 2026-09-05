#!/usr/bin/env python3
"""
Panel D: CEACAM5 versus CEACAM6 in pre-treatment stomach epithelium.

Replaces the per-cell scatter of create_ceacam_correlation.py. At this
sequencing depth a per-cell scatter is dominated by dropout: two thirds of the
cells in the shallowest quintile are double-negative, and the per-cell rho of
0.44 is a floor set by detection rather than a measurement of co-expression.
The panel therefore shows kNN metacells of 10 cells built within each sample,
and reports the per-cell value beside the metacell value so that neither is
mistaken for the other.

The pooling and the data path are copied from
04_Revision_Analyses/04_R1.5_CEACAM5_vs_CEACAM6/scripts/dropout_and_coexpression.py,
which also carries the control this panel does not draw: pooling at random
gives rho = 0.71 at k = 10, so the rise from 0.44 is largely the arithmetic of
averaging. The co-expression claim rests on the within-stratum odds ratio in
that module (3.4 to 12.2 across depth quintiles), not on this panel.

Panel is drawn to the slot measured in the submitted Figure 2: x 80.5-118.2 mm,
y 6.7-31.3 mm.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD  # noqa: E402
from shared.figure_config import use_panel_style

OUT = Path(__file__).parent
PANEL_W_MM, PANEL_H_MM = 34.0, 24.0
MM_TO_INCH = 1 / 25.4
DPI = 300
K = 10
N_PCS = 30
SEED = 42


def load():
    """Pre-treatment stomach epithelium, expression from .raw."""
    ad = sc.read_h5ad(EPITHELIAL_H5AD)
    ad = ad[(ad.obs["Sample site"] == "Stomach")
            & (ad.obs["Treatment phase"] == "Pre")].copy()
    # .raw where the file has one; the clean deposit promotes .raw.X to .X and
    # keeps no .raw, so there .X is that same log1p matrix.
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
    use_panel_style(font_pt=6, scale=1, **{'axes.labelsize': 6, 'xtick.labelsize': 5, 'ytick.labelsize': 5, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2, 'axes.spines.top': False, 'axes.spines.right': False})

    df, emb = load()
    print(f"Pre-treatment stomach epithelial cells: {len(df):,} "
          f"from {df['sample'].nunique()} samples")

    rho_cell, p_cell = stats.spearmanr(df["CEACAM5"], df["CEACAM6"])
    print(f"per cell        rho = {rho_cell:.4f}  P = {p_cell:.3g}")

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
    print(f"metacells k={K}  n = {len(mc):,}  rho = {rho_mc:.4f}  P = {p_mc:.3g}")

    agg = df.groupby("sample")[["CEACAM5", "CEACAM6"]].mean()
    rho_s, p_s = stats.spearmanr(agg["CEACAM5"], agg["CEACAM6"])
    print(f"per sample n={len(agg)}  rho = {rho_s:.4f}  P = {p_s:.3g}")

    mc.to_csv(OUT / "ceacam_metacell_values.csv", index=False)
    pd.DataFrame([
        dict(level="cell", n=len(df), rho=float(rho_cell), p=float(p_cell)),
        dict(level=f"metacell_k{K}", n=len(mc), rho=float(rho_mc), p=float(p_mc)),
        dict(level="sample", n=len(agg), rho=float(rho_s), p=float(p_s)),
    ]).to_csv(OUT / "ceacam_correlation_levels.csv", index=False)

    fig, ax = plt.subplots(figsize=(PANEL_W_MM * MM_TO_INCH,
                                    PANEL_H_MM * MM_TO_INCH))
    ax.scatter(mc["CEACAM5"], mc["CEACAM6"], s=1.4, alpha=0.45, c="#3498db",
               edgecolors="none")
    ax.set_xlabel(r"$\it{CEACAM5}$")
    ax.set_ylabel(r"$\it{CEACAM6}$")

    hi = float(np.ceil(max(mc["CEACAM5"].max(), mc["CEACAM6"].max())))
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)
    ax.set_xticks(np.arange(0, hi + 1, 2))
    ax.set_yticks(np.arange(0, hi + 1, 2))

    ax.text(0.04, 0.97,
            f"ρ = {rho_mc:.2f}, metacells of {K} (n = {len(mc):,})\n"
            f"ρ = {rho_cell:.2f} per cell (n = {len(df):,})",
            transform=ax.transAxes, fontsize=4.6, va="top", linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85,
                      edgecolor="none"))

    plt.tight_layout(pad=0.3)
    for ext in ("svg", "pdf", "png"):
        f = OUT / f"ceacam_metacell_correlation.{ext}"
        plt.savefig(f, dpi=DPI, bbox_inches="tight", facecolor="white")
        print("Saved:", f)
    plt.close()


if __name__ == "__main__":
    main()
