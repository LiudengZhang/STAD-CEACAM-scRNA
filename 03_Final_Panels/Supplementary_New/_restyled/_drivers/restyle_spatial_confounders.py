#!/usr/bin/env python3
"""
Restyle S9A and S9B - spatial confounder adjustment (Reviewer 1, point R1.6).

Analysis: 04_Revision_Analyses/05_R1.6_Spatial_Confounders/scripts/spatial_confounders.py

NOTHING IS RECOMPUTED HERE. The analysis's `main()` fits six mixed models for
the sequential adjustment, six more for the density strata and one OLS per
sample, and writes each set to a CSV before handing it to a `_panel*` function.
This driver imports the module for its constants and reads back the two tables
those panels were given:

    _panel_adjustment(adjusted)   <- outputs/spatial_adjusted_models.csv
    _panel_stratified(strat)      <- outputs/spatial_stratified.csv

Confirmed against `main()`: `adjusted` is the frame written at line 94 and
`strat` the frame written at line 112, in that order and with those columns.
No mixed model is refitted, `SPATIAL_SPOT_DATA` is never opened, and the
per-sample table is not touched because neither panel draws it.

    conda run -n Liudeng_Python_310 python restyle_spatial_confounders.py --check
    conda run -n Liudeng_Python_310 python restyle_spatial_confounders.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _driver_base as base
import panel_style_cns as style

A = base.analysis("05_R1.6_Spatial_Confounders/scripts/spatial_confounders.py")
OUT = base.outputs_of(A)
FIG = "S9_Mechanism_Specificity"

# Printed boxes, millimetres. S9A grew from 171 x 44 because its y tick labels
# are set at 7 pt here against 6.72 pt as printed in Version A and its legend
# sits below the x label. S9B shrank from 118 x 42 - its type grew too, but
# Version A reserved 20% of the canvas below the tick labels for nothing.
A_W, A_H = 171.0, 46.0
B_W, B_H = 118.0, 38.0

# Mark sizes, rescaled so each mark keeps its Version A size *relative to the
# type beside it*. See DRIVER_SPEC.md step 4.
#   S9A: 9.0 x 3.6 cm at SCALE 4, fitted into 171 x 44 mm at 0.3056.
#        Smallest type is the y tick labels and legend at 5.5 * SCALE.
A_FIT = 0.3056
A_TYPE = 5.5 * 4 * A_FIT
A_MARK = (style.tick_pt() / A_TYPE) * A_FIT
A_AREA = A_MARK ** 2
#   S9B: 8.0 x 3.8 cm at SCALE 4, fitted into 118 x 42 mm at 0.2763.
#        Smallest type is the stratum tick labels at 5 * SCALE.
B_FIT = 0.2763
B_TYPE = 5.0 * 4 * B_FIT
B_MARK = (style.tick_pt() / B_TYPE) * B_FIT


def frames():
    """The two tables the analysis already wrote. No computation."""
    return dict(
        adjusted=pd.read_csv(base.require(
            OUT / "spatial_adjusted_models.csv", "S9A sequential adjustment")),
        strat=pd.read_csv(base.require(
            OUT / "spatial_stratified.csv", "S9B density-stratified models")),
    )


# --------------------------------------------------------------------------
# S9A - sequential adjustment, forest plot. Layout unchanged (single axes).
# --------------------------------------------------------------------------
def draw_A(fr, save=True):
    style.apply()
    adjusted = fr["adjusted"]
    fig, ax = style.subplots_mm(A_W, A_H)
    labels, y = [], []
    # Version A's colour map, keyed by the outcome strings in the table.
    colors = {"Distance to immune-rich regions": "#7B3294",
              "Distance to stroma": "#1B7837"}
    k = 0
    for outcome in adjusted["outcome"].unique():
        for _, r in adjusted[adjusted["outcome"] == outcome].iterrows():
            ax.plot([r["ci_low"], r["ci_high"]], [k, k],
                    color=colors[outcome], linewidth=1.2 * A_MARK)
            ax.scatter(r["ceacam_coef"], k, s=26 * A_AREA, color=colors[outcome],
                       zorder=3, edgecolors="white", linewidths=0.4 * A_MARK)
            # The three model labels repeat for each outcome, so the first row
            # of each block carries the outcome name to disambiguate them.
            short = ("immune" if "immune" in outcome else "stroma")
            labels.append(f"{short}:  {r['model']}")
            y.append(k)
            k += 1
        k += 0.8
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("CEACAM ratio coefficient (SD of distance per unit)")
    ax.tick_params(axis="x", width=0.6, length=2)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    handles = [plt.Line2D([], [], color=c, linewidth=1.2 * A_MARK, marker="o",
                          markersize=3 * A_MARK, label=o)
               for o, c in colors.items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.26),
              ncol=2, handlelength=1.6, handletextpad=0.5, columnspacing=2.0)
    style.margins_mm(fig, left=48, right=3, top=3, bottom=16)
    if save:
        base.save(fig, FIG, "S9_A", "S9_A_spatial_adjusted_models")
    return fig


# --------------------------------------------------------------------------
# S9B - CEACAM coefficient within tertiles of local epithelial density.
# Layout unchanged (1 x 2).
# --------------------------------------------------------------------------
def draw_B(fr, save=True):
    style.apply()
    strat = fr["strat"]
    outcomes = list(strat["outcome"].unique())
    fig, axes = style.subplots_mm(B_W, B_H, 1, len(outcomes))
    for ax, outcome in zip(np.atleast_1d(axes), outcomes):
        sub = strat[strat["outcome"] == outcome]
        xs = np.arange(len(sub))
        ax.bar(xs, sub["ceacam_coef"], color="#7B3294", alpha=0.75,
               edgecolor="#444444", linewidth=0.5, width=0.6)
        ax.errorbar(xs, sub["ceacam_coef"],
                    yerr=[sub["ceacam_coef"] - sub["ci_low"],
                          sub["ci_high"] - sub["ceacam_coef"]],
                    fmt="none", ecolor="#333333", elinewidth=0.8 * B_MARK,
                    capsize=2 * B_MARK)
        ax.axhline(0, color="#999999", linewidth=0.5)
        ax.set_xticks(xs)
        ax.set_xticklabels([s.replace(" epithelial density", "\nepith. density")
                            for s in sub["stratum"]])
        ax.set_title(outcome)
        ax.tick_params(axis="both", width=0.6, length=2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    np.atleast_1d(axes)[0].set_ylabel("CEACAM ratio coefficient")
    style.margins_mm(fig, left=13, right=3, top=5, bottom=6, wspace=0.34)
    if save:
        base.save(fig, FIG, "S9_B", "S9_B_spatial_density_stratified")
    return fig


def main():
    fr = frames()
    return base.run({
        "S9_A": (lambda: draw_A(fr), lambda: A._panel_adjustment(fr["adjusted"])),
        "S9_B": (lambda: draw_B(fr), lambda: A._panel_stratified(fr["strat"])),
    })


if __name__ == "__main__":
    sys.exit(main())
