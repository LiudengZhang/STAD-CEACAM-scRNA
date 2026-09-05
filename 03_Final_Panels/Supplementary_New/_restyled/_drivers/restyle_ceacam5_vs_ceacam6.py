#!/usr/bin/env python3
"""
Restyle S8B and S8C - CEACAM5 versus CEACAM6 (Reviewer 1, point R1.5).

Analysis: 04_Revision_Analyses/04_R1.5_CEACAM5_vs_CEACAM6/scripts/ceacam5_vs_ceacam6.py

NOTHING IS RECOMPUTED HERE. The driver imports that module for its constants
(STATES, COLOR_R/COLOR_NR, SCALE) and for `mw`, and reads the three tables its
`main()` already wrote. It never opens Epithelial.h5ad, so it never reaches the
`.raw`/`.X` guard at line 78 that 00_Data_Audit/FINDINGS.md sections 1 and 7
are about, and it never touches the IHC deconvolution CSV.

Two things are re-executed rather than reproduced, because re-executing them is
calling the *same function object* the original panel called:

  S8B  the original `_panel_states` recomputes its printed P inside the drawing
       loop with the module's own `mw()`. So does this driver, via `A.mw`.
       Reimplementing a Mann-Whitney here would be exactly the move the spec
       forbids.
  S8C  the original reads its P out of the tests table. So does this driver.

`piv` for S8C is read back from `ihc_per_marker_values.csv`. That file is
written at line 169, *before* the tests loop and before `_panel_ihc(piv, ...)`
is called, and nothing between the write and the call mutates `piv`. It carries
"patient", "group" and all five measure columns, so the three columns the panel
draws - "CEACAM5", "CEACAM6", "Summed (published)" - come back in exactly the
shape the panel was given. No pivot is rebuilt and the IHC source CSV is not
read.

    conda run -n Liudeng_Python_310 python restyle_ceacam5_vs_ceacam6.py
    conda run -n Liudeng_Python_310 python restyle_ceacam5_vs_ceacam6.py --check
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _driver_base as base
import panel_style_cns as style

A = base.analysis("04_R1.5_CEACAM5_vs_CEACAM6/scripts/ceacam5_vs_ceacam6.py")
OUT = base.outputs_of(A)
FIG = "S8_CEACAM_Metaprogram"

# Printed boxes, millimetres. Version A: S8B 171 x 40, S8C 124 x 36. Both grow
# in height only - the type is now 7/8 pt instead of the 4.76/4.71 pt these two
# panels actually printed at, and the tick labels, the "P = ..." bridge label
# and the axis titles all need the room. Widths are unchanged.
B_W, B_H = 171.0, 50.0
C_W, C_H = 124.0, 46.0

# Data-mark rescale, per DRIVER_SPEC "Restyling the drawing code", step 4. The
# marks must keep their size *relative to the type*, and the type is changing.
#
#   A_TYPE = smallest fontsize in the original x SCALE x the assembler fit
#          = what that type actually printed at in Version A
#
# S8B's smallest is the 5 pt "P = ..." label, S8C's the 5.5 pt tick labels.
B_FIT, C_FIT = 0.2381, 0.2143                       # measured 2026-09-01
B_TYPE = 5.0 * A.SCALE * B_FIT                      # 4.762 pt as printed
C_TYPE = 5.5 * A.SCALE * C_FIT                      # 4.715 pt as printed
B_MARK = (style.tick_pt() / B_TYPE) * B_FIT
C_MARK = (style.tick_pt() / C_TYPE) * C_FIT
B_AREA, C_AREA = B_MARK ** 2, C_MARK ** 2

# Box, whisker and cap outlines are structure rather than a mark whose size
# encodes a value, so they take cnsplots' own line weight the way the spines
# do; the median keeps its emphasis over them.
BOX_LW, MEDIAN_LW = 0.5, 0.8


def frames():
    """The tables the analysis already wrote. No computation."""
    return dict(
        frac=pd.read_csv(base.require(
            OUT / "ceacam_state_fractions.csv", "S8B state fractions")),
        piv=pd.read_csv(base.require(
            OUT / "ihc_per_marker_values.csv", "S8C per-marker IHC values")),
        ihc_tests=pd.read_csv(base.require(
            OUT / "ihc_per_marker_tests.csv", "S8C per-marker IHC tests")),
    )


# --------------------------------------------------------------------------
# S8B - single- and double-positive fractions, 1 x 4 boxplots. Layout unchanged.
#
# The jitter draws from `np.random.default_rng(1)` in the original order -
# state by state, responders before non-responders - so every point lands where
# it landed in Version A.
# --------------------------------------------------------------------------
def draw_B(fr, save=True):
    style.apply()
    frac = fr["frac"]
    fig, axes = style.subplots_mm(B_W, B_H, 1, len(A.STATES))
    rng = np.random.default_rng(1)
    for ax, s in zip(axes, A.STATES):
        r = frac.loc[frac["group"] == "R", s].values * 100
        nr = frac.loc[frac["group"] == "NR", s].values * 100
        bp = ax.boxplot([r, nr], positions=[0, 1], widths=0.55,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(linewidth=BOX_LW),
                        whiskerprops=dict(linewidth=BOX_LW),
                        capprops=dict(linewidth=BOX_LW),
                        medianprops=dict(color="black", linewidth=MEDIAN_LW))
        bp["boxes"][0].set_facecolor(A.COLOR_R); bp["boxes"][0].set_alpha(0.55)
        bp["boxes"][1].set_facecolor(A.COLOR_NR); bp["boxes"][1].set_alpha(0.55)
        for i, (vals, c) in enumerate(((r, A.COLOR_R), (nr, A.COLOR_NR))):
            ax.scatter(i + rng.uniform(-0.1, 0.1, len(vals)), vals,
                       s=9 * A.SCALE * B_AREA, c=c, zorder=3,
                       edgecolors="white",
                       linewidths=0.3 * A.SCALE * B_MARK)
        # The module's own Mann-Whitney, called the way the original called it.
        p, _ = A.mw(nr, r)
        top = max(np.max(r), np.max(nr))
        ax.plot([0, 0, 1, 1], [top * 1.06, top * 1.11, top * 1.11, top * 1.06],
                color="#444444", linewidth=BOX_LW)
        ax.text(0.5, top * 1.13, f"P = {p:.3f}", ha="center", va="bottom",
                fontsize=style.tick_pt())
        ax.set_ylim(0, top * 1.30)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["R", "NR"])
        ax.set_title(s)
        ax.tick_params(axis="both", width=0.6, length=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("% of pre-treatment\nepithelial cells")
    style.margins_mm(fig, left=16, right=2, top=6, bottom=8, wspace=0.42)
    if save:
        base.save(fig, FIG, "S8_B", "S8_B_ceacam_single_double_positive")
    return fig


# --------------------------------------------------------------------------
# S8C - IHC, each marker separately against the summed composite.
#
# `cols` is a literal in the original drawing function rather than a module
# constant, so it is copied verbatim; the P for each column is read from the
# tests table exactly as the original reads it. Jitter from `default_rng(2)`,
# column by column, responders first.
# --------------------------------------------------------------------------
def draw_C(fr, save=True):
    style.apply()
    piv, tests = fr["piv"], fr["ihc_tests"]
    cols = ["CEACAM5", "CEACAM6", "Summed (published)"]
    fig, axes = style.subplots_mm(C_W, C_H, 1, 3)
    rng = np.random.default_rng(2)
    for ax, col in zip(axes, cols):
        r = piv.loc[piv["group"] == "R", col].values
        nr = piv.loc[piv["group"] == "NR", col].values
        bp = ax.boxplot([r, nr], positions=[0, 1], widths=0.55,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(linewidth=BOX_LW),
                        whiskerprops=dict(linewidth=BOX_LW),
                        capprops=dict(linewidth=BOX_LW),
                        medianprops=dict(color="black", linewidth=MEDIAN_LW))
        bp["boxes"][0].set_facecolor(A.COLOR_R); bp["boxes"][0].set_alpha(0.55)
        bp["boxes"][1].set_facecolor(A.COLOR_NR); bp["boxes"][1].set_alpha(0.55)
        for i, (vals, c) in enumerate(((r, A.COLOR_R), (nr, A.COLOR_NR))):
            ax.scatter(i + rng.uniform(-0.1, 0.1, len(vals)), vals,
                       s=11 * A.SCALE * C_AREA, c=c, zorder=3,
                       edgecolors="white",
                       linewidths=0.3 * A.SCALE * C_MARK)
        p = tests.loc[tests["measure"] == col, "p_two_tailed"].iloc[0]
        top = max(np.max(r), np.max(nr))
        ax.plot([0, 0, 1, 1], [top * 1.05, top * 1.09, top * 1.09, top * 1.05],
                color="#444444", linewidth=BOX_LW)
        ax.text(0.5, top * 1.11, f"P = {p:.3f}", ha="center", va="bottom",
                fontsize=style.tick_pt())
        ax.set_ylim(0, top * 1.28)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["R", "NR"])
        ax.set_title(col)
        ax.tick_params(axis="both", width=0.6, length=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("DAB+ area (% of tissue)")
    fig.suptitle("Immunohistochemistry, n = 4 responders vs 4 non-responders "
                 "(two-sided Mann-Whitney)", fontsize=style.tick_pt())
    style.margins_mm(fig, left=17, right=2, top=10, bottom=8, wspace=0.42)
    if save:
        base.save(fig, FIG, "S8_C", "S8_C_ihc_per_marker")
    return fig


def main():
    fr = frames()
    return base.run({
        "S8_B": (lambda: draw_B(fr), lambda: A._panel_states(fr["frac"])),
        "S8_C": (lambda: draw_C(fr),
                 lambda: A._panel_ihc(fr["piv"], fr["ihc_tests"])),
    })


if __name__ == "__main__":
    sys.exit(main())
