#!/usr/bin/env python3
"""
Restyle S11A - affirmative spatial evidence (Reviewer 1, point R1.6).

Analysis: 02_New_Analyses/05_R1.6_Spatial_Confounders/scripts/spatial_positive_evidence.py

NOTHING IS RECOMPUTED HERE. `main()` (line 466) calls three functions and hands
their return values straight to `_panel(med, t, g)` (line 504). Each of the
three writes the frame it returns to a CSV before returning it, so all three
reconstruct from disk:

    med = mediation()    -> returns `med`, written to outputs/spatial_mediation.csv
    t   = tiger()        -> returns `t`,   written to outputs/tiger_immune_exclusion.csv
    g   = gse246011()    -> returns `g`,   written to outputs/gse246011_replication.csv

That matters more here than anywhere else in the set, because re-running any of
them would be expensive and, in one case, unsafe:

  * `mediation()` fits four mixed models and then bootstraps the indirect effect
    200 times per outcome by refitting OLS on a resampled cohort. It is also
    seeded from a module-level `RNG`, so a second call in the same process would
    not even reproduce the first;
  * `tiger()` reads the BayesPrism fraction and epithelial-expression tables and
    fits fifteen models;
  * `gse246011()` opens four spatial h5ads, normalises them and builds a kd-tree.
    DRIVER_SPEC forbids a driver opening an h5ad, and this one must not: its own
    comment records that `.X` in those sections is a z-scored matrix and only
    `.raw` carries counts.

None of that happens. The three CSVs carry every column `_panel` draws -
`med`: indirect_effect / direct_effect / indirect_ci_lo / indirect_ci_hi;
`t`: cell_type / n / beta_purity_adjusted; `g`: model / ceacam_coef / ci_low /
ci_high - and nothing else is read.

The four other tables `main()` writes - spatial_compositionality.csv,
spatial_architecture_adjusted.csv, tiger_estimate_harmonised.csv and
tiger_ceacam_purity_adjusted.csv - are not touched, because S11A does not draw
them.

    conda run -n Liudeng_Python_310 python restyle_spatial_positive_evidence.py --check
    conda run -n Liudeng_Python_310 python restyle_spatial_positive_evidence.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _driver_base as base
import panel_style_cns as style

A = base.analysis("05_R1.6_Spatial_Confounders/scripts/spatial_positive_evidence.py")
OUT = base.outputs_of(A)
FIG = "S11_Affirmative_Analyses"

# Printed box, millimetres. Version A was 171 x 42. The eight TIGER cell-type
# labels in the middle axes were printed at 4.72 pt and are 7 pt here, so the
# stack of eight bars needs the extra depth, and the mediation legend no longer
# fits inside its axes (see draw_A).
W, H = 171.0, 52.0

# Version A: 12.0 x 4.0 cm at SCALE 4, fitted into 171 x 42 mm at 0.2625.
# Smallest type is the middle axes' cell-type labels at 4.5 * SCALE.
A_FIT = 0.2625
A_TYPE = 4.5 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT


def frames():
    """The three tables `main()` already wrote, in the order `_panel` takes them."""
    return dict(
        med=pd.read_csv(base.require(
            OUT / "spatial_mediation.csv", "S11A mediation (med)")),
        t=pd.read_csv(base.require(
            OUT / "tiger_immune_exclusion.csv", "S11A TIGER bulk cohort (t)")),
        g=pd.read_csv(base.require(
            OUT / "gse246011_replication.csv", "S11A GSE246011 replication (g)")),
    )


def draw_A(fr, save=True):
    style.apply()
    med, t, g = fr["med"], fr["t"], fr["g"]
    fig, axes = style.subplots_mm(W, H, 1, 3)

    ax = axes[0]
    x = np.arange(len(med))
    w = 0.38
    ax.bar(x - w / 2, med["indirect_effect"], width=w, color="#7B3294",
           edgecolor="#333", linewidth=0.5, label="Indirect (via density)")
    ax.bar(x + w / 2, med["direct_effect"], width=w, color="#c2a5cf",
           edgecolor="#333", linewidth=0.5, label="Direct")
    ax.errorbar(x - w / 2, med["indirect_effect"],
                yerr=[med["indirect_effect"] - med["indirect_ci_lo"],
                      med["indirect_ci_hi"] - med["indirect_effect"]],
                fmt="none", ecolor="#333", elinewidth=0.8 * MARK,
                capsize=2 * MARK)
    ax.axhline(0, color="#666", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(["to immune", "to stroma"])
    ax.set_ylabel("Effect on distance (SD)")
    ax.set_title("Mediation by epithelial density")
    # Version A let matplotlib place this legend ("best"), which at 5.09 pt
    # fitted inside the axes. At 7 pt it covers the "to stroma" error bar, and
    # the y limit cannot be raised to make room because a limit is content.
    # So the legend moves below the axes. Same two entries, same labels.
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
              handlelength=1.2, handletextpad=0.4, columnspacing=1.4)

    ax = axes[1]
    # `t` is read from disk and is never None here; base.require has already
    # raised if the table is missing, so the analysis's `if t is not None`
    # guard - which would silently draw an empty axes - has nothing to do.
    s = t.sort_values("beta_purity_adjusted")
    y = np.arange(len(s))
    cols = ["#B2182B" if v < 0 else "#2166AC" for v in s["beta_purity_adjusted"]]
    ax.barh(y, s["beta_purity_adjusted"], color=cols, edgecolor="#333",
            linewidth=0.4, height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace("Monocytes Macrophages", "Mono/Mac")
                        for c in s["cell_type"]])
    ax.axvline(0, color="#666", linewidth=0.5)
    ax.set_xlabel("Fraction change per log2 CEACAM\n(purity-adjusted)")
    ax.set_title(f"PRJEB25780, n = {int(s['n'].iloc[0])}")

    ax = axes[2]
    y = np.arange(len(g))
    ax.barh(y, g["ceacam_coef"], color="#7B3294", edgecolor="#333",
            linewidth=0.4, height=0.6)
    ax.errorbar(g["ceacam_coef"], y,
                xerr=[g["ceacam_coef"] - g["ci_low"],
                      g["ci_high"] - g["ceacam_coef"]],
                fmt="none", ecolor="#333", elinewidth=0.8 * MARK,
                capsize=2 * MARK)
    ax.set_yticks(y)
    ax.set_yticklabels(g["model"])
    ax.axvline(0, color="#666", linewidth=0.5)
    ax.set_xlabel("CEACAM coefficient")
    ax.set_title(f"GSE246011 replication")

    for ax in axes:
        ax.tick_params(axis="both", width=0.6, length=2)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    style.margins_mm(fig, left=15, right=3, top=6, bottom=13, wspace=0.58)
    if save:
        base.save(fig, FIG, "S11_A", "S11_A_spatial_positive_evidence")
    return fig


def main():
    fr = frames()
    return base.run({
        "S11_A": (lambda: draw_A(fr),
                  lambda: A._panel(fr["med"], fr["t"], fr["g"])),
    })


if __name__ == "__main__":
    sys.exit(main())
