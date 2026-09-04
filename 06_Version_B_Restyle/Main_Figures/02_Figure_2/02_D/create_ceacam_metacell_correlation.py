#!/usr/bin/env python3
"""
Figure 2, printed panel D, RESTYLED (Version B) - CEACAM5 versus CEACAM6 in
pre-treatment stomach epithelium.

Replaces the per-cell scatter of create_ceacam_correlation.py. At this
sequencing depth a per-cell scatter is dominated by dropout: two thirds of the
cells in the shallowest quintile are double-negative, and the per-cell rho of
0.44 is a floor set by detection rather than a measurement of co-expression.
The panel therefore shows kNN metacells of 10 cells built within each sample,
and reports the per-cell value beside the metacell value so that neither is
mistaken for the other.

The pooling and the data path are copied from
02_New_Analyses/04_R1.5_CEACAM5_vs_CEACAM6/scripts/dropout_and_coexpression.py,
which also carries the control this panel does not draw: pooling at random
gives rho = 0.71 at k = 10, so the rise from 0.44 is largely the arithmetic of
averaging. The co-expression claim rests on the within-stratum odds ratio in
that module (3.4 to 12.2 across depth quintiles), not on this panel.

Version A is
`03_Final_Panels/02_Figure_2/02_D/create_ceacam_metacell_correlation.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`. K, N_PCS, SEED,
the sample filter, the pooling and both correlations are Version A's.

  printed panel  Figure 2 D     (PROVENANCE.csv; the directory letter happens
                                 to agree here - it was still looked up)
  printed rect   30.9 x 23.6 mm    (panel_rects.csv)
  Version B box  62.0 x 48.0 mm

MARK
    Version A drew this one at 1:1, not at 4x: `use_panel_style(font_pt=6,
    scale=1, ...)` on a 34.0 x 24.0 mm canvas, so SCALE = 1. Its smallest body
    type is the two-line rho annotation at `fontsize=4.6`.

        SCALE = 1, SMALL_PT = 4.6
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 4.6 = 1.52174
        AREA  = MARK ** 2

    The one non-type size is the metacell dot area, `s = 1.4`, which takes AREA.

    Version A's rcParams block (spine and tick widths, tick sizes) is NOT
    carried over: that is axes furniture, cnsplots has its own settings for it,
    and following the library rather than rescaling the old numbers is the
    standard-methods rule.

    Why the panel grew so much: the annotation inside the axes is a 37-character
    line ("rho = 0.93, metacells of 10 (n = ...)"). At 4.6 pt it fitted a 31 mm
    panel; at 7 pt it needs about 46 mm of axes width on its own, and the y
    label plus its tick labels take another 11 mm. Nothing about what is drawn
    changed - only the width the same string occupies at legible type.
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

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import EPITHELIAL_H5AD  # noqa: E402
import panel_style_cns as style    # noqa: E402

OUT = Path(__file__).parent

PRINTED_MM = (30.9, 23.6)           # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 62.0, 48.0
MARGIN = dict(left=11.0, right=2.5, top=2.5, bottom=9.5)

SCALE = 1                           # Version A's canvas multiplier (not 4 here)
SMALL_PT = 4.6                      # Version A's smallest body type

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
    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

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

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    ax.scatter(mc["CEACAM5"], mc["CEACAM6"], s=1.4 * AREA, alpha=0.45,
               c="#3498db", edgecolors="none")
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
            transform=ax.transAxes, fontsize=style.tick_pt(), va="top",
            linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85,
                      edgecolor="none"))

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUT / "ceacam_metacell_correlation")
    print(f"Saved: {OUT / 'ceacam_metacell_correlation'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
