#!/usr/bin/env python3
"""
Restyle S11B - TCGA-STAD immune exclusion (Reviewer 1, point R1.6).

Analysis: 02_New_Analyses/05_R1.6_Spatial_Confounders/scripts/tcga_immune_exclusion.py

NOTHING IS RECOMPUTED HERE, and in particular ESTIMATE IS NOT RE-RUN. `main()`
(line 103) builds three things and passes them to `_panel(d, t, purity_col)`
(line 207); each is reconstructed from what `main()` already put on disk or
from the stored tables it read, never from the raw counts:

  `t`   is the model table, written verbatim to outputs/tcga_immune_exclusion.csv.
        Six Spearman tests and six purity-adjusted OLS fits; all read, none refit.

  `d`   is the per-tumour frame. `main()` assembles it from
          - `est`, the ESTIMATE score matrix, which it writes to
            outputs/tcga_estimate_scores.csv *before* truncating its barcodes.
            That CSV is read here instead of calling `run_estimate()`, which
            shells out to R and needs the 147 MB tumour-only count matrix and
            the Ensembl-to-symbol map from a raw GDC file. Neither is opened.
          - `epi`, the stored BayesPrism epithelial expression table, read with
            the same two lines `main()` uses. This is a deconvolution result
            that already exists on disk, not a raw count file and not a
            recomputation: only `main()`'s own log2(mean + 1) is applied.
          - `ab`, TCGA PanCanAtlas ABSOLUTE purity, likewise read as `main()`
            reads it.
        The barcode truncation, the intersection and the column order all
        follow `main()` line for line, so `d` has the same 407 rows in the same
        order. Verified: no duplicate 12-character barcode in either table.

  `purity_col` is hardcoded by `main()` - line 148, `purity_col =
        "ABSOLUTE_purity"` - not derived from the columns present, so it is
        hardcoded here too, next to the column that carries it. `_panel` takes
        the argument and does not use it; it is passed for fidelity.

    conda run -n Liudeng_Python_310 python restyle_tcga_immune_exclusion.py --check
    conda run -n Liudeng_Python_310 python restyle_tcga_immune_exclusion.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _driver_base as base
import panel_style_cns as style

A = base.analysis("05_R1.6_Spatial_Confounders/scripts/tcga_immune_exclusion.py")
OUT = base.outputs_of(A)
FIG = "S11_Affirmative_Analyses"

# Printed box, millimetres. Version A was 118 x 42. The left axes carries a
# two-line title and a two-line x label and the right one a three-line x label;
# at 8 pt those need 8 mm more depth than they did at 6.3 pt.
W, H = 118.0, 50.0

# Version A: 8.0 x 4.0 cm at SCALE 4, fitted into 118 x 42 mm at 0.2625.
# Smallest type is the exposure tick labels and the x labels at 5.5 * SCALE.
A_FIT = 0.2625
A_TYPE = 5.5 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2

# `main()` hardcodes this at line 148; it is not chosen from the columns found.
PURITY_COL = "ABSOLUTE_purity"


def frames():
    """`d` and `t` as `main()` built them. No model is refitted, ESTIMATE is not re-run."""
    t = pd.read_csv(base.require(
        OUT / "tcga_immune_exclusion.csv", "S11B purity-adjusted model table"))

    est = pd.read_csv(base.require(
        OUT / "tcga_estimate_scores.csv", "S11B stored ESTIMATE scores"),
        index_col=0)
    epi = pd.read_csv(base.require(
        A.TCGA / "tcga_bayesprism_epithelial_expression.tsv",
        "S11B BayesPrism deconvolved epithelial expression"),
        sep="\t", index_col=0)

    # main() lines 128-140, verbatim.
    epi.index = [str(i)[:12] for i in epi.index]
    est.index = [str(i)[:12] for i in est.index]
    common = epi.index.intersection(est.index)
    d = pd.DataFrame(index=common)
    d["CEACAM"] = np.log2(epi.loc[common, ["CEACAM5", "CEACAM6"]].mean(axis=1) + 1)
    d["CEACAM5"] = np.log2(epi.loc[common, "CEACAM5"] + 1)
    d["CEACAM6"] = np.log2(epi.loc[common, "CEACAM6"] + 1)
    for c in est.columns:
        d[c] = est.loc[common, c]

    # main() lines 145-148. The panel never reads this column, but `d` is the
    # frame main() built and the driver does not quietly ship a different one.
    absolute = pd.read_csv(base.require(
        A.ABSOLUTE_PURITY, "S11B TCGA PanCanAtlas ABSOLUTE purity"), sep="\t")
    absolute["patient"] = absolute["array"].astype(str).str[:12]
    ab = (absolute.dropna(subset=["purity"])
          .groupby("patient")["purity"].mean())
    d[PURITY_COL] = ab.reindex(d.index)
    return dict(d=d, t=t)


def draw_B(fr, save=True):
    style.apply()
    d, t = fr["d"], fr["t"]
    fig, axes = style.subplots_mm(W, H, 1, 2)

    ax = axes[0]
    ax.scatter(d["CEACAM"], d["ImmuneScore"], s=4 * 4 * AREA, c="#4d4d4d",
               alpha=0.45, edgecolors="none")
    # The fit is the analysis's own line, copied verbatim so --check can prove
    # the same regression line is drawn.
    z = np.polyfit(d["CEACAM"].dropna(),
                   d.loc[d["CEACAM"].notna(), "ImmuneScore"], 1)
    xs = np.linspace(d["CEACAM"].min(), d["CEACAM"].max(), 50)
    ax.plot(xs, np.polyval(z, xs), color="#B2182B", linewidth=1.2 * MARK)
    ax.set_xlabel("Epithelial CEACAM5/6\n(log2, deconvolved)")
    ax.set_ylabel("ESTIMATE immune score")
    row = t[(t["exposure"] == "CEACAM") & (t["outcome"] == "ImmuneScore")]
    ax.set_title(f"TCGA-STAD, n = {int(row['n'].iloc[0])}\n"
                 f"purity-adjusted P = {row['p_purity_adjusted'].iloc[0]:.3g}")

    ax = axes[1]
    sub = t[t["outcome"] == "ImmuneScore"]
    y = np.arange(len(sub))
    cols = ["#B2182B" if v < 0 else "#2166AC" for v in sub["beta_purity_adjusted"]]
    ax.barh(y, sub["beta_purity_adjusted"], color=cols, edgecolor="#333",
            linewidth=0.4, height=0.6)
    ax.errorbar(sub["beta_purity_adjusted"], y,
                xerr=[sub["beta_purity_adjusted"] - sub["ci_lo"],
                      sub["ci_hi"] - sub["beta_purity_adjusted"]],
                fmt="none", ecolor="#333", elinewidth=0.8 * MARK,
                capsize=2 * MARK)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["exposure"])
    ax.axvline(0, color="#666", linewidth=0.5)
    ax.set_xlabel("Immune score change\nper log2 CEACAM\n(purity-adjusted)")

    for ax in axes:
        ax.tick_params(axis="both", width=0.6, length=2)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    style.margins_mm(fig, left=15, right=3, top=8, bottom=13, wspace=0.60)
    if save:
        base.save(fig, FIG, "S11_B", "S11_B_tcga_immune_exclusion")
    return fig


def main():
    fr = frames()
    return base.run({
        "S11_B": (lambda: draw_B(fr),
                  lambda: A._panel(fr["d"], fr["t"], PURITY_COL)),
    })


if __name__ == "__main__":
    sys.exit(main())
