#!/usr/bin/env python3
"""
Draw S9E - NF-kB regulon activity in monocytes and macrophages
(Reviewer 1, point R1.8, affirmative evidence).

Analysis: 04_Revision_Analyses/07_R1.8_NFkB_Specificity/scripts/nfkb_regulon_activity.py

NOTHING IS RECOMPUTED HERE. The driver does not call `load_regulons()` and so
never opens the stored AUCell matrix. That matters beyond speed: the deposited
AUCell matrix is a subset of the cells, so a driver that re-derived sample-level
means from it could quietly draw a different panel from the one `main()` drew.
It also never opens a h5ad, so `score_genes` is not re-run and the feedback
target scores are not re-derived.

Instead it reads the two tables `main()` already wrote that this panel needs:

    outputs/nfkb_sample_values.csv      -> sv, the per-sample values that the
                                          boxes and the points actually plot
    outputs/nfkb_regulon_activity.csv   -> reg_df, from which the panel takes
                                          the post-treatment NFKB1(+)/MoMac
                                          P value it prints as its title

ONE BLOCK OF THE ORIGINAL PANEL, NOT BOTH. The analysis's `_panel` builds a
1 x 2 figure: the left axes is the NFKB1 regulon in monocytes and macrophages
at both timepoints, the right one the NF-kB negative-feedback target score by
compartment. Only the left block is drawn here; the right one is withdrawn from
the figure set.

The two blocks share no state. `_panel` filters `sv` twice and independently -
`source == "regulon"` for the left, `source == "panel"` for the right - and
rebinds `series` and `colors` to fresh lists between them, so nothing the right
block builds is read by the left. `_box` takes its axes as an argument, seeds
its own generator per call and returns nothing. `fb` is accepted by `_panel`
and never read at all. Dropping the right block therefore changes no value the
left block draws, and `--check` compares this panel against exactly the left
axes of the original to prove it.

`_box` is copied in below as a local helper. It is drawing code only - it
computes nothing except the jitter, from its own seeded generator - and its
behaviour is unchanged.

    python draw_nfkb_regulon_activity.py --check
    python draw_nfkb_regulon_activity.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _driver_base as base
import panel_style_cns as style

base.apply_style()

A = base.analysis("07_R1.8_NFkB_Specificity/scripts/nfkb_regulon_activity.py")
OUT = base.outputs_of(A)
FIG = "S9_MoMac_Identity_NFkB"

# The printed box, millimetres. Four sample-level boxes, each with a two-line
# tick label under it, a two-line y axis label and a P value as the title. Half
# the width of the two-block original, because half of it is what is drawn.
W, H = 64.0, 50.0

# Previously 9.0 x 4.6 cm drawn at four times print size and fitted into a
# 171 x 44 mm box at 0.2391; the smallest type, `_box`'s x tick labels, was set
# at 5.0 x that scale and printed at 4.782 pt. Every mark and stroke below is
# rescaled by the factor that keeps its size relative to the type unchanged.
A_FIT = 0.2391
A_TYPE = 5.0 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2


def frames():
    """The tables the analysis already wrote. No computation."""
    return dict(
        reg_df=pd.read_csv(base.require(
            OUT / "nfkb_regulon_activity.csv", "S9E regulon activity")),
        fb=pd.read_csv(base.require(
            OUT / "nfkb_feedback_targets.csv", "S9E feedback targets")),
        sv=pd.read_csv(base.require(
            OUT / "nfkb_sample_values.csv", "S9E per-sample values")),
    )


def _box(ax, series, colors, ylabel, title=None):
    """`nfkb_regulon_activity._box`, redrawn. Same marks, same numbers.

    One sample-level box per group with every sample drawn on top. The groups
    here have five or six samples, so the points are the honest display and the
    box is only there to carry the median and the spread.
    """
    rng = np.random.default_rng(0)
    labels = [lab for lab, _ in series]
    data = [vals for _, vals in series]
    bp = ax.boxplot(data, widths=0.55, showfliers=False, patch_artist=True,
                    medianprops=dict(color="#333333", linewidth=1.0 * MARK),
                    whiskerprops=dict(color="#666666", linewidth=0.8 * MARK),
                    capprops=dict(color="#666666", linewidth=0.8 * MARK),
                    boxprops=dict(linewidth=0.6 * MARK, edgecolor="#333333"))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.35)
    for i, (vals, c) in enumerate(zip(data, colors), start=1):
        if not len(vals):
            continue
        jitter = rng.uniform(-0.13, 0.13, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=9 * 4 * AREA, c=c,
                   edgecolors="white", linewidths=0.4 * MARK, zorder=3)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels([f"{lab}\nn = {len(v)}" for lab, v in zip(labels, data)])
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="both", width=0.6, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    if title:
        ax.set_title(title)


def draw_E(fr, save=True):
    base.apply_style()
    reg_df, sv = fr["reg_df"], fr["sv"]
    fig, ax = style.subplots_mm(W, H)

    # NFKB1 regulon in monocytes/macrophages, both timepoints. The selection is
    # the original's, verbatim.
    reg = sv[(sv["source"] == "regulon") & (sv["key"] == "NFKB1(+)")
             & (sv["stratum"] == "MoMac")]
    series, colors = [], []
    for phase in ("Pre", "Post"):
        for grp, c in (("R", A.COLOR_R), ("NR", A.COLOR_NR)):
            v = reg.loc[(reg["phase"] == phase) & (reg["group"] == grp),
                        "value"].values
            if len(v):
                series.append((f"{phase}\n{grp}", v))
                colors.append(c)
    post = reg_df[(reg_df["regulon"] == "NFKB1(+)")
                  & (reg_df["cell_type"] == "MoMac")
                  & (reg_df["phase"] == "Post")]
    title = (f"post-treatment P = {post['p_two_tailed'].iloc[0]:.3f}"
             if len(post) else None)
    _box(ax, series, colors,
         "NFKB1 regulon activity (AUCell)\nmonocytes/macrophages", title)

    base.fit(fig)
    if save:
        base.save(fig, FIG, "S9_E", "S9_E_nfkb_regulon_momac")
    return fig


def main():
    fr = frames()
    return base.run({
        # Axis 0 of the original is the block this panel keeps; the comparison
        # is made over that axis and says nothing about the one withdrawn.
        "S9_E": (lambda: draw_E(fr),
                 lambda: A._panel(fr["reg_df"], fr["fb"], fr["sv"]),
                 None, [0]),
    })


if __name__ == "__main__":
    sys.exit(main())
