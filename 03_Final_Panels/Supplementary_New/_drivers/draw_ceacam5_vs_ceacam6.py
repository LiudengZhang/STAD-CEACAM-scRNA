#!/usr/bin/env python3
"""Supplementary Figure S4 B and C: CEACAM5 versus CEACAM6, pre-treatment R
against NR (two-sided Mann-Whitney U). B: single- and double-positive
epithelial fractions; C: IHC DAB+ area per marker and summed.

Reads the three tables written by
04_Revision_Analyses/04_R1.5_CEACAM5_vs_CEACAM6/scripts/ceacam5_vs_ceacam6.py;
writes S4_CEACAM5_vs_CEACAM6/S4_B and S4_C.

Run: python draw_ceacam5_vs_ceacam6.py [--check]"""

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
FIG = "S4_CEACAM5_vs_CEACAM6"
from shared.figure_config import PRE_PAIR  # noqa: E402
COLOR_R, COLOR_NR = PRE_PAIR      # pre-treatment colours

# Panel sizes in mm.
B_W, B_H = 171.0, 46.0
C_W, C_H = 124.0, 46.0

def frames():
    """The three tables the analysis wrote."""
    return dict(
        frac=pd.read_csv(base.require(
            OUT / "ceacam_state_fractions.csv", "S7B state fractions")),
        piv=pd.read_csv(base.require(
            OUT / "ihc_per_marker_values.csv", "S7C per-marker IHC values")),
        ihc_tests=pd.read_csv(base.require(
            OUT / "ihc_per_marker_tests.csv", "S7C per-marker IHC tests")),
    )

def draw_B(fr, save=True):
    """S4 B: single- and double-positive fractions, one box pair per state."""
    base.apply_style()
    frac = fr["frac"]
    fig, axes = style.subplots_mm(B_W, B_H, 1, len(A.STATES))
    for ax, s in zip(axes, A.STATES):
        r = frac.loc[frac["group"] == "R", s].values * 100
        nr = frac.loc[frac["group"] == "NR", s].values * 100
        bxp = boxes.draw_boxes(ax, [r, nr], [0, 1], [COLOR_R, COLOR_NR],
                               width=0.55)
        p, _ = A.mw(nr, r)
        top = max(np.max(r), np.max(nr))
        _, text, _ = boxes.bracket(fig, ax, 0, 1, top, top, p, kind="pair",
                                   lift=0.06, arm=0.05)
        ax.set_ylim(0, top * 1.30)
        boxes.ylim_above(ax, text)
        ax.set_xlim(*boxes.box_xlim([0, 1], 0.55))
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Pre-R", "Pre-NR"])
        ax.set_title(s)
        base.finish(ax)
        boxes.assert_no_points(ax, bxp)
    axes[0].set_ylabel("% of pre-treatment\nepithelial cells")
    base.fit(fig, wspace=0.42)
    if save:
        base.save(fig, FIG, "S4_B", "S4_B_ceacam_single_double_positive")
    return fig

def draw_C(fr, save=True):
    """S4 C: IHC, each marker separately and the summed composite."""
    base.apply_style()
    piv, tests = fr["piv"], fr["ihc_tests"]
    cols = ["CEACAM5", "CEACAM6", "Summed (published)"]
    fig, axes = style.subplots_mm(C_W, C_H, 1, 3)
    for ax, col in zip(axes, cols):
        r = piv.loc[piv["group"] == "R", col].values
        nr = piv.loc[piv["group"] == "NR", col].values
        bxp = boxes.draw_boxes(ax, [r, nr], [0, 1], [COLOR_R, COLOR_NR],
                               width=0.55)
        p = tests.loc[tests["measure"] == col, "p_two_tailed"].iloc[0]
        top = max(np.max(r), np.max(nr))
        _, text, _ = boxes.bracket(fig, ax, 0, 1, top, top, p, kind="pair",
                                   lift=0.05, arm=0.04)
        ax.set_ylim(0, top * 1.28)
        boxes.ylim_above(ax, text)
        ax.set_xlim(*boxes.box_xlim([0, 1], 0.55))
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Pre-R", "Pre-NR"])
        ax.set_title("Summed" if col == "Summed (published)" else col)
        base.finish(ax)
        boxes.assert_no_points(ax, bxp)
    axes[0].set_ylabel("DAB+ area (% of tissue)")
    base.fit(fig, wspace=0.42)
    if save:
        base.save(fig, FIG, "S4_C", "S4_C_ihc_per_marker")
    return fig

def main():
    fr = frames()
    return base.run({
        "S4_B": (lambda: draw_B(fr), lambda: A._panel_states(fr["frac"]),
                 base.aliases),
        "S4_C": (lambda: draw_C(fr),
                 lambda: A._panel_ihc(fr["piv"], fr["ihc_tests"]),
                 base.aliases),
    })

if __name__ == "__main__":
    sys.exit(main())
