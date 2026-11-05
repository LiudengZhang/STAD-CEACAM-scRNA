#!/usr/bin/env python3
"""Draw Supplementary Figure S7 D (CEACAM5/6 vs ESTIMATE immune score in
TCGA-STAD) and S7 E (purity-adjusted coefficients).

Reads the model table outputs/tcga_immune_exclusion.csv and rebuilds the
per-tumour frame from the ESTIMATE scores, BayesPrism epithelial expression
and ABSOLUTE purity, as
04_Revision_Analyses/05_R1.6_Spatial_Confounders/scripts/tcga_immune_exclusion.py
does. No model is refit.


    python draw_tcga_immune_exclusion.py [--check]"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _driver_base as base
import panel_style_cns as style

base.apply_style()

A = base.analysis("05_R1.6_Spatial_Confounders/scripts/tcga_immune_exclusion.py")
OUT = base.outputs_of(A)
FIG = "S7_Spatial_Confounders"

# Panel size, mm.
W, H = 80.0, 54.0

# Mark sizes, rescaled to keep their size relative to the type.
A_FIT = 0.2625
A_TYPE = 5.5 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2

PURITY_COL = "ABSOLUTE_purity"

def frames():
    """Per-tumour frame `d` and model table `t`, assembled as in the analysis."""
    t = pd.read_csv(base.require(
        OUT / "tcga_immune_exclusion.csv", "S9D purity-adjusted model table"))

    est = pd.read_csv(base.require(
        OUT / "tcga_estimate_scores.csv", "S9D stored ESTIMATE scores"),
        index_col=0)
    epi = pd.read_csv(base.require(
        A.TCGA_BAYESPRISM_EPI,
        "S9D BayesPrism deconvolved epithelial expression"),
        sep="\t", index_col=0)

    epi.index = [str(i)[:12] for i in epi.index]
    est.index = [str(i)[:12] for i in est.index]
    common = epi.index.intersection(est.index)
    d = pd.DataFrame(index=common)
    d["CEACAM"] = np.log2(epi.loc[common, ["CEACAM5", "CEACAM6"]].mean(axis=1) + 1)
    d["CEACAM5"] = np.log2(epi.loc[common, "CEACAM5"] + 1)
    d["CEACAM6"] = np.log2(epi.loc[common, "CEACAM6"] + 1)
    for c in est.columns:
        d[c] = est.loc[common, c]

    absolute = pd.read_csv(base.require(
        A.ABSOLUTE_PURITY, "S9D PanCanAtlas ABSOLUTE purity"), sep="\t")
    absolute["patient"] = absolute["array"].astype(str).str[:12]
    ab = (absolute.dropna(subset=["purity"])
          .groupby("patient")["purity"].mean())
    d[PURITY_COL] = ab.reindex(d.index)
    return dict(d=d, t=t)

def draw_D(fr, save=True):
    """S7 D: scatter with a least-squares line."""
    base.apply_style()
    d, t = fr["d"], fr["t"]
    fig, ax = style.subplots_mm(W, H)
    ax.scatter(d["CEACAM"], d["ImmuneScore"], s=4 * 4 * AREA, c="#4d4d4d",
               alpha=0.45, edgecolors="none")
    z = np.polyfit(d["CEACAM"].dropna(),
                   d.loc[d["CEACAM"].notna(), "ImmuneScore"], 1)
    xs = np.linspace(d["CEACAM"].min(), d["CEACAM"].max(), 50)
    ax.plot(xs, np.polyval(z, xs), color="#B2182B", linewidth=style.RULE_PT)
    ax.set_xlabel("Epithelial CEACAM5/6\n(log2, deconvolved)")
    ax.set_ylabel("ESTIMATE immune score")
    row = t[(t["exposure"] == "CEACAM") & (t["outcome"] == "ImmuneScore")]
    # All tumours with expression and ESTIMATE scores are plotted; the model
    # drops those without a purity call, and the title gives the model's n.
    ax.set_title(f"TCGA-STAD, purity-adjusted\n"
                 f"{style.p_exact(row['p_purity_adjusted'].iloc[0])} "
                 f"(n = {int(row['n'].iloc[0])})")
    # Unadjusted Spearman rho inside the box.
    from cnsfig.layout import corr_annotate
    corr_annotate(ax, row['spearman_rho'].iloc[0], loc="upper right")
    base.finish(ax)
    base.fit(fig)
    if save:
        base.save(fig, FIG, "S7_D", "S7_D_tcga_immune_exclusion")
    return fig

def draw_E(fr, save=True):
    """S7 E: purity-adjusted coefficients with 95% CIs."""
    base.apply_style()
    t = fr["t"]
    fig, ax = style.subplots_mm(W, H)
    sub = t[t["outcome"] == "ImmuneScore"]
    y = np.arange(len(sub))
    cols = ["#B2182B" if v < 0 else "#2166AC" for v in sub["beta_purity_adjusted"]]
    ax.barh(y, sub["beta_purity_adjusted"], color=cols, edgecolor="#333",
            linewidth=style.EDGE_PT, height=0.6)
    ax.errorbar(sub["beta_purity_adjusted"], y,
                xerr=[sub["beta_purity_adjusted"] - sub["ci_lo"],
                      sub["ci_hi"] - sub["beta_purity_adjusted"]],
                fmt="none", ecolor="#333", elinewidth=style.RULE_PT,
                capsize=2 * MARK)
    ax.set_yticks(y)
    # Composite exposure (mean of CEACAM5 and CEACAM6) is labelled CEACAM5/6.
    ax.set_yticklabels(["CEACAM5/6" if e == "CEACAM" else e
                        for e in sub["exposure"]])
    ax.axvline(0, color="#666", linewidth=style.RULE_PT)
    ax.set_xlabel("Immune score change\nper log2 CEACAM5/6\n(purity-adjusted)")
    base.finish(ax)
    base.fit(fig)
    if save:
        base.save(fig, FIG, "S7_E", "S7_E_tcga_immune_exclusion_coefficients")
    return fig

def main():
    fr = frames()
    # --check: D against axis 0 of the analysis panel, E against axis 1.
    original = lambda: A._panel(fr["d"], fr["t"], PURITY_COL)   # noqa: E731
    return base.run({
        "S7_D": (lambda: draw_D(fr), original, None, [0]),
        "S7_E": (lambda: draw_E(fr), original, None, [1]),
    })

if __name__ == "__main__":
    sys.exit(main())

