#!/usr/bin/env python3
"""
Draw S8A - spatial confounder adjustment (Reviewer 1, point R1.6).

Analysis: 04_Revision_Analyses/05_R1.6_Spatial_Confounders/scripts/spatial_confounders.py

NOTHING IS RECOMPUTED HERE. The analysis's `main()` fits six mixed models for
the sequential adjustment and writes them to a CSV before handing that frame to
`_panel_adjustment`. This driver imports the module for its constants and reads
that table back:

    _panel_adjustment(adjusted)   <- outputs/spatial_adjusted_models.csv

No mixed model is refitted, the spatial spot table is never opened, and the
per-sample table is not touched because this panel does not draw it.

    python draw_spatial_confounders.py --check
    python draw_spatial_confounders.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import matplotlib.pyplot as plt

import _driver_base as base
import panel_style_cns as style

base.apply_style()

A = base.analysis("05_R1.6_Spatial_Confounders/scripts/spatial_confounders.py")
OUT = base.outputs_of(A)
FIG = "S8_Spatial_Confounders"

# The printed box, millimetres. Full width across the 183 mm page less its 6 mm
# margins. Six model rows, their labels set on the y axis, and a two-entry
# legend below the x axis label set the depth.
W, H = 171.0, 46.0

# Mark sizes, rescaled so each mark keeps its size relative to the type beside
# it. Previously 9.0 x 3.6 cm drawn at four times print size and fitted into a
# 171 x 44 mm box at 0.3056; the smallest type was the y tick labels and the
# legend at 5.5 x that scale.
A_FIT = 0.3056
A_TYPE = 5.5 * 4 * A_FIT
A_MARK = (style.tick_pt() / A_TYPE) * A_FIT
A_AREA = A_MARK ** 2


def frames():
    """The table the analysis already wrote. No computation."""
    return dict(
        adjusted=pd.read_csv(base.require(
            OUT / "spatial_adjusted_models.csv", "S8A sequential adjustment")),
    )


# --------------------------------------------------------------------------
# S8A - sequential adjustment, forest plot.
# --------------------------------------------------------------------------
def draw_A(fr, save=True):
    base.apply_style()
    adjusted = fr["adjusted"]
    fig, ax = style.subplots_mm(W, H)
    labels, y = [], []
    # The colour map, keyed by the outcome strings in the table.
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
    base.fit(fig)
    if save:
        base.save(fig, FIG, "S8_A", "S8_A_spatial_adjusted_models")
    return fig


def main():
    fr = frames()
    return base.run({
        "S8_A": (lambda: draw_A(fr), lambda: A._panel_adjustment(fr["adjusted"])),
    })


if __name__ == "__main__":
    sys.exit(main())
