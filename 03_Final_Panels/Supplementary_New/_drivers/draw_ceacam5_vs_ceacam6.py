#!/usr/bin/env python3
"""
Draw S7B and S7C - CEACAM5 versus CEACAM6 (Reviewer 1, point R1.5).

Analysis: 04_Revision_Analyses/04_R1.5_CEACAM5_vs_CEACAM6/scripts/ceacam5_vs_ceacam6.py

NOTHING IS RECOMPUTED HERE. The driver imports that module for its constants
(STATES, COLOR_R/COLOR_NR, SCALE) and for `mw`, and reads the three tables its
`main()` already wrote. It never opens Epithelial.h5ad, so it never reaches the
`.raw`/`.X` guard the data audit is about, and it never touches the IHC
deconvolution CSV.

Two things are re-executed rather than reproduced, because re-executing them is
calling the *same function object* the original panel called:

  S7B  the original `_panel_states` recomputes its printed P inside the drawing
       loop with the module's own `mw()`. So does this driver, via `A.mw`.
       Reimplementing a Mann-Whitney here would move the statistic out of the
       file that owns it.
  S7C  the original reads its P out of the tests table. So does this driver.

`piv` for S7C is read back from `ihc_per_marker_values.csv`. That file is
written before the tests loop and before `_panel_ihc(piv, ...)` is called, and
nothing between the write and the call mutates `piv`. It carries "patient",
"group" and all five measure columns, so the three columns the panel draws -
"CEACAM5", "CEACAM6", "Summed (published)" - come back in exactly the shape the
panel was given. No pivot is rebuilt and the IHC source CSV is not read.

    python draw_ceacam5_vs_ceacam6.py --check
    python draw_ceacam5_vs_ceacam6.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _driver_base as base
import panel_style_cns as style
from cnsfig import boxes

base.apply_style()

A = base.analysis("04_R1.5_CEACAM5_vs_CEACAM6/scripts/ceacam5_vs_ceacam6.py")
OUT = base.outputs_of(A)
FIG = "S8_CEACAM5_vs_CEACAM6"

# Printed boxes, millimetres. S7B is full width - four boxplot pairs across the
# 183 mm page less its 6 mm margins - and its depth is set by a two-line y axis
# label beside a bridge annotation above the tallest box. S7C carries three
# boxplot pairs under a wrapped figure-level caption, so it is set narrower and
# a little deeper per column.
B_W, B_H = 171.0, 46.0
C_W, C_H = 124.0, 46.0

# Data-mark rescale: the marks keep their size *relative to the type*, and the
# type is set here rather than inherited from a four-times canvas.
#
#   A_TYPE = the smallest fontsize in the original x its draw scale x the fit
#          = what that type actually printed at
#
# S7B's smallest is the 5 pt "P = ..." label, S7C's the 5.5 pt tick labels.
B_FIT, C_FIT = 0.2381, 0.2143
B_TYPE = 5.0 * A.SCALE * B_FIT                      # 4.762 pt as printed
C_TYPE = 5.5 * A.SCALE * C_FIT                      # 4.715 pt as printed
B_MARK = (style.tick_pt() / B_TYPE) * B_FIT
C_MARK = (style.tick_pt() / C_TYPE) * C_FIT
B_AREA, C_AREA = B_MARK ** 2, C_MARK ** 2

# ONE BOX, NO POINTS (2026-09-16, the author's fifth reading). The boxes are
# cnsfig.boxes.draw_boxes (Figure 2 H/I's: RULE_PT black lines, 0.79 mm open
# fliers, solid faces in the analysis's R/NR colours, as Figure 5 J/L print
# them) and the bracket is cnsfig.boxes.bracket, black, its label's ink 0.4 mm
# above the line. The jittered individual samples are no longer drawn ("box
# plots should look alike throughout; don't show every point"); the fliers
# now stand in for the values outside the whiskers. B_MARK/C_MARK are kept for
# the record of what the marks were rescaled by until then.


def frames():
    """The tables the analysis already wrote. No computation."""
    return dict(
        frac=pd.read_csv(base.require(
            OUT / "ceacam_state_fractions.csv", "S7B state fractions")),
        piv=pd.read_csv(base.require(
            OUT / "ihc_per_marker_values.csv", "S7C per-marker IHC values")),
        ihc_tests=pd.read_csv(base.require(
            OUT / "ihc_per_marker_tests.csv", "S7C per-marker IHC tests")),
    )


# --------------------------------------------------------------------------
# S7B - single- and double-positive fractions, 1 x 4 boxplots. The individual
# samples were jittered on top until 2026-09-16; the boxes and the P are
# unchanged.
# --------------------------------------------------------------------------
def draw_B(fr, save=True):
    base.apply_style()
    frac = fr["frac"]
    fig, axes = style.subplots_mm(B_W, B_H, 1, len(A.STATES))
    for ax, s in zip(axes, A.STATES):
        r = frac.loc[frac["group"] == "R", s].values * 100
        nr = frac.loc[frac["group"] == "NR", s].values * 100
        bxp = boxes.draw_boxes(ax, [r, nr], [0, 1], [A.COLOR_R, A.COLOR_NR],
                               width=0.55)
        # The module's own Mann-Whitney, called the way the original called it.
        p, _ = A.mw(nr, r)
        top = max(np.max(r), np.max(nr))
        # The original's bracket: foot at 1.06 x top, arms to 1.11 x top.
        _, text, _ = boxes.bracket(fig, ax, 0, 1, top, top, p, kind="pair",
                                   lift=0.06, arm=0.05)
        ax.set_ylim(0, top * 1.30)
        boxes.ylim_above(ax, text)
        ax.set_xlim(*boxes.box_xlim([0, 1], 0.55))
        ax.set_xticks([0, 1]); ax.set_xticklabels(["R", "NR"])
        ax.set_title(s)
        ax.tick_params(axis="both", width=0.6, length=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        boxes.assert_no_points(ax, bxp)
    axes[0].set_ylabel("% of pre-treatment\nepithelial cells")
    base.fit(fig, wspace=0.42)
    if save:
        base.save(fig, FIG, "S8_B", "S8_B_ceacam_single_double_positive")
    return fig


# --------------------------------------------------------------------------
# S7C - IHC, each marker separately against the summed composite.
#
# `cols` is a literal in the original drawing function rather than a module
# constant, so it is copied verbatim; the P for each column is read from the
# tests table exactly as the original reads it. The jittered samples were
# withdrawn on 2026-09-16.
# --------------------------------------------------------------------------
def draw_C(fr, save=True):
    base.apply_style()
    piv, tests = fr["piv"], fr["ihc_tests"]
    # The stored table keeps the analysis-era column name; the panel uses the
    # concise display label requested for the supplement.
    cols = ["CEACAM5", "CEACAM6", "Summed (published)"]
    fig, axes = style.subplots_mm(C_W, C_H, 1, 3)
    for ax, col in zip(axes, cols):
        r = piv.loc[piv["group"] == "R", col].values
        nr = piv.loc[piv["group"] == "NR", col].values
        bxp = boxes.draw_boxes(ax, [r, nr], [0, 1], [A.COLOR_R, A.COLOR_NR],
                               width=0.55)
        p = tests.loc[tests["measure"] == col, "p_two_tailed"].iloc[0]
        top = max(np.max(r), np.max(nr))
        # The original's bracket: foot at 1.05 x top, arms to 1.09 x top.
        _, text, _ = boxes.bracket(fig, ax, 0, 1, top, top, p, kind="pair",
                                   lift=0.05, arm=0.04)
        ax.set_ylim(0, top * 1.28)
        boxes.ylim_above(ax, text)
        ax.set_xlim(*boxes.box_xlim([0, 1], 0.55))
        ax.set_xticks([0, 1]); ax.set_xticklabels(["R", "NR"])
        ax.set_title("Summed" if col == "Summed (published)" else col)
        ax.tick_params(axis="both", width=0.6, length=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        boxes.assert_no_points(ax, bxp)
    axes[0].set_ylabel("DAB+ area (% of tissue)")
    base.fit(fig, wspace=0.42)
    if save:
        base.save(fig, FIG, "S8_C", "S8_C_ihc_per_marker")
    return fig


def main():
    fr = frames()
    # The P strings are declared renames (labels.py RENAMES_SUPPLEMENTARY):
    # the analysis sets three decimals, the drawing p_label's two.
    return base.run({
        "S8_B": (lambda: draw_B(fr), lambda: A._panel_states(fr["frac"]),
                 base.aliases),
        "S8_C": (lambda: draw_C(fr),
                 lambda: A._panel_ihc(fr["piv"], fr["ihc_tests"]),
                 base.aliases),
    })


if __name__ == "__main__":
    sys.exit(main())
