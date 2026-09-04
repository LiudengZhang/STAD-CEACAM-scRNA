#!/usr/bin/env python3
"""
Restyle S8H - dropout and CEACAM5/CEACAM6 co-expression (Reviewer 1, R1.5).

Analysis: 02_New_Analyses/04_R1.5_CEACAM5_vs_CEACAM6/scripts/dropout_and_coexpression.py

NOTHING IS RECOMPUTED HERE. The driver imports that module for its constants
(KS, the four colours, SCALE) and reads the three tables its `main()` already
wrote:

    strata  outputs/dropout_depth_strata.csv        (main() line 325)
    ws      outputs/coexpression_within_sample.csv  (main() line 337)
    sweep   outputs/metacell_sweep.csv              (main() line 360)

`panel()` uses `ws` only for `len(ws)` in the suptitle, so the whole panel is
reconstructible from disk. It never opens Epithelial.h5ad, never builds a kNN
metacell cover and never draws a permutation, so the `default_rng(SEED)` stream
is not consumed and cannot move.

    conda run -n Liudeng_Python_310 python restyle_dropout_and_coexpression.py
    conda run -n Liudeng_Python_310 python restyle_dropout_and_coexpression.py --check
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _driver_base as base
import panel_style_cns as style

A = base.analysis("04_R1.5_CEACAM5_vs_CEACAM6/scripts/dropout_and_coexpression.py")
OUT = base.outputs_of(A)
FIG = "S8_CEACAM_Metaprogram"

# --------------------------------------------------------------------------
# LAYOUT CHANGED - 1 x 3 becomes 3 x 1. Content identical.
#
# S8H had the worst fit of any panel in the set: 680 x 200 mm of canvas fitted
# into a 171 x 32 mm box, a scale of 0.1600, so its 4.6 pt annotations printed
# at 2.94 pt and its axis titles at 3.84 pt. Set at 8 pt, the second axis title
# - "Co-detected far above independence, at every depth" - measures 66.3 mm
# (measured 2026-09-01, Nimbus Sans). Three cells across a 171 mm page are at
# most 57 mm wide before any margin, and about 44 mm after them, so that title
# overhangs its cell by roughly 11 mm on each side and collides with the axis
# titles either side of it. The first axis title needs 48.7 mm in the same
# 44 mm.
#
# The five "N.NN x expected" labels on the middle axis are the same story:
# 16.7 mm each at 7 pt, against a 44 mm axis that has to hold five of them.
#
# Stacking the three axes instead gives each one the full 151 mm of drawing
# width, which holds every title and spaces the five labels 34 mm apart. The
# three axes, their data, their limits, their scales, their tick labels, their
# annotations and their legends are untouched.
# --------------------------------------------------------------------------
# Height. This panel was first drawn 181 mm tall because the annotation offsets
# below were being read as panel content and so had to stay at Version A's
# point values - 8 pt and 18 pt on a canvas a quarter the size, which forced a
# tall axis to keep the labels apart. That reading was a fault in the harness,
# not a property of the panel: an Annotation's content is the point it points
# at, and `xytext` is layout. Since compare_panel_content was corrected the
# offsets rescale with the type like every other length here, and the three
# axes hold their labels at 30 mm each instead of 44 mm.
#
# That matters beyond this panel. S8 carries seven panels, and at 181 mm this
# one alone was 38% of the page; the figure came to 481 mm, which no printable
# area accepts without scaling the type back down to 3.6 pt. See
# RESTYLE_REPORT.md.
W, H = 171.0, 139.0

# Data-mark rescale, per DRIVER_SPEC step 4. The smallest type in Version A was
# 4.6 pt (legends and annotations), which at the 0.1600 fit printed at 2.94 pt.
FIT = 0.1600                                        # measured 2026-09-01
A_TYPE = 4.6 * A.SCALE * FIT                        # 2.944 pt as printed
MARK = (style.tick_pt() / A_TYPE) * FIT

# Reference rules take cnsplots' own line weight, as the spines do.
RULE_LW = 0.5


def frames():
    """The three tables the analysis already wrote. No computation."""
    return dict(
        strata=pd.read_csv(base.require(
            OUT / "dropout_depth_strata.csv", "S8H depth strata")),
        ws=pd.read_csv(base.require(
            OUT / "coexpression_within_sample.csv", "S8H within-sample rho")),
        sweep=pd.read_csv(base.require(
            OUT / "metacell_sweep.csv", "S8H metacell sweep")),
    )


def draw_H(fr, save=True):
    style.apply()
    strata, ws, sweep = fr["strata"], fr["ws"], fr["sweep"]
    fig, axes = style.subplots_mm(W, H, 3, 1)

    ax = axes[0]
    x = np.arange(len(strata))
    ax.plot(x, strata["rho"], "-o", color=A.COLOR_REAL, lw=1.2 * A.SCALE * 0.5 * MARK,
            ms=4 * A.SCALE * 0.5 * MARK, label="Spearman $\\rho$, per cell")
    ax.plot(x, strata["pct_double_of_expressing"] / 100, "-s", color=A.COLOR_5,
            lw=1.2 * A.SCALE * 0.5 * MARK, ms=4 * A.SCALE * 0.5 * MARK,
            label="double positive, of cells\nexpressing either gene")
    ax.plot(x, strata["pct_single_pos"] / 100, "-^", color=A.COLOR_6,
            lw=1.2 * A.SCALE * 0.5 * MARK, ms=4 * A.SCALE * 0.5 * MARK,
            label="single positive, of all cells")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:,.0f}" for v in strata["median_total_counts"]],
                       rotation=0)
    ax.set_xlabel("median UMI per cell, quintile")
    ax.set_ylabel("proportion, or $\\rho$")
    ax.set_title("Deeper cells look more double positive")
    ax.legend(frameon=False, loc="upper left", ncol=3, handlelength=1.4,
              handletextpad=0.5, columnspacing=1.6)
    ax.set_ylim(0, 1)

    # The odds ratio rather than observed/expected: the expected value itself
    # climbs with detection rate, so the ratio shrinks even as the association
    # strengthens. The odds ratio does not have that dependence.
    #
    # The annotation offsets are rescaled with the type, like every other
    # point-specified length in this driver: MARK, not a bare Version A value.
    # See the note on W, H above for why they were once left unscaled.
    ax = axes[1]
    ax.plot(x, strata["odds_ratio"], "-o", color=A.COLOR_REAL,
            lw=1.2 * A.SCALE * 0.5 * MARK, ms=4 * A.SCALE * 0.5 * MARK)
    ax.axhline(1.0, color="#999999", lw=RULE_LW, ls="--")
    ax.annotate("independence", (x[0], 1.0), textcoords="offset points",
                xytext=(2 * MARK, 4 * A.SCALE * MARK), fontsize=style.tick_pt(),
                color="#777777")
    # Labels sit below the line and alternate side, so they clear the y axis
    # on the left and the last point on the right.
    for i, r in strata.reset_index().iterrows():
        last = i == len(strata) - 1
        ax.annotate(f"{r['obs_over_expected']:.2f}× expected",
                    (i, r["odds_ratio"]), textcoords="offset points",
                    xytext=(-4 * A.SCALE * MARK if last
                            else 4 * A.SCALE * MARK,
                            -9 * A.SCALE * MARK),
                    ha="right" if last else "left", va="top",
                    fontsize=style.tick_pt())
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:,.0f}" for v in strata["median_total_counts"]])
    ax.set_yscale("log")
    ax.set_yticks([1, 2, 5, 10, 20])
    ax.set_yticklabels(["1", "2", "5", "10", "20"])
    ax.set_ylim(0.8, 30)
    ax.set_xlabel("median UMI per cell, quintile")
    ax.set_ylabel("odds ratio for co-detection")
    ax.set_title("Co-detected far above independence, at every depth")

    ax = axes[2]
    s = sweep.dropna(subset=["k"])
    ax.plot(s["k"], s["rho"], "-o", color=A.COLOR_REAL,
            lw=1.2 * A.SCALE * 0.5 * MARK, ms=4 * A.SCALE * 0.5 * MARK,
            label="kNN metacells")
    ax.plot(s["k"], s["rho_random_pools"], "-o", color=A.COLOR_SHUF,
            lw=1.2 * A.SCALE * 0.5 * MARK, ms=4 * A.SCALE * 0.5 * MARK,
            label="random pools, same size")
    lvl = float(sweep.loc[sweep["k"].isna(), "rho"].iloc[0])
    ax.axhline(lvl, color=A.COLOR_5, lw=RULE_LW, ls=":")
    # Anchored at k = KS[1] and 22 pt below the rule, both of which are fixed.
    # Set to the left of the anchor rather than the right: at 7 pt the label is
    # 34 mm long, and to the right of k = 5 that is exactly where both curves
    # climb through it. To the left they are still at rho 0.44-0.66, well below.
    ax.annotate(f"whole-sample means, $\\rho$ = {lvl:.2f}", (A.KS[1], lvl),
                textcoords="offset points", xytext=(0, -11 * A.SCALE * MARK),
                ha="right", fontsize=style.tick_pt(), color=A.COLOR_5)
    ax.set_xscale("log")
    ax.set_xticks(A.KS)
    ax.set_xticklabels([str(k) for k in A.KS])
    ax.set_xlabel("cells pooled per metacell")
    ax.set_ylabel("Spearman $\\rho$")
    ax.set_title("Pooling raises $\\rho$, neighbours or not")
    ax.legend(frameon=False, loc="lower right", handlelength=1.4,
              handletextpad=0.5)
    ax.set_ylim(0, 1)

    for ax in axes:
        ax.tick_params(axis="both", width=0.6, length=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.suptitle(
        f"Pre-treatment stomach epithelial cells (n = {int(strata['n_cells'].sum()):,}"
        f", {len(ws)} samples); expression from the log1p CP10K matrix",
        fontsize=style.tick_pt())
    # top 10 mm holds the figure title over the first axis title; the 14 mm
    # gap between rows holds the axis label of the row above and the axis title
    # of the row below. 181 = 10 + 3 x 44 + 2 x 14 + 11.
    style.margins_mm(fig, left=19, right=3, top=10, bottom=11,
                     hspace=14.0 / 44.0)
    if save:
        base.save(fig, FIG, "S8_H", "S8_H_dropout_and_coexpression")
    return fig


def main():
    fr = frames()
    return base.run({
        "S8_H": (lambda: draw_H(fr),
                 lambda: A.panel(fr["strata"], fr["ws"], fr["sweep"])),
    })


if __name__ == "__main__":
    sys.exit(main())
