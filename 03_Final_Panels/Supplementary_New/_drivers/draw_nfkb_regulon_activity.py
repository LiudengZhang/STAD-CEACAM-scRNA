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

`_box` is copied in below as a local helper. It is drawing code only. Since
2026-09-16 (the author's fifth reading) it sets its boxes with
`cnsfig.boxes.draw_boxes` - the one box of the paper - and no longer jitters
the individual samples on top; the data reach ax.boxplot untouched, so every
box statistic is the original's, and the comparator sees the removed point
collections and the added 0.79 mm fliers and nothing else.

    python draw_nfkb_regulon_activity.py --check
    python draw_nfkb_regulon_activity.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

import _driver_base as base
import panel_style_cns as style
from cnsfig import boxes
from cnsfig.rich import rich_ylabel

base.apply_style()

A = base.analysis("07_R1.8_NFkB_Specificity/scripts/nfkb_regulon_activity.py")
OUT = base.outputs_of(A)
FIG = "S10_MoMac_Identity_NFkB"

# The printed box, millimetres. Four sample-level boxes, each with a two-line
# tick label under it, a two-line y axis label and a P value as the title. Half
# the width of the two-block original, because half of it is what is drawn.
# W 64 -> 37 on 2026-09-16 (the author's fifth reading): E shares a row with C
# and D under B, so the page reads A, B, C-D-E in order; four boxes at a 0.55
# width still clear the spines (cnsfig.boxes.box_xlim).
W, H = 37.0, 60.0

# Previously 9.0 x 4.6 cm drawn at four times print size and fitted into a
# 171 x 44 mm box at 0.2391; the smallest type, `_box`'s x tick labels, was set
# at 5.0 x that scale and printed at 4.782 pt. Until 2026-09-16 every mark and
# stroke was rescaled by the factor that keeps its size relative to the type;
# the boxes now come from cnsfig.boxes at the paper's one weight, and MARK is
# kept for the record of what the marks were rescaled by until then.
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
    """`nfkb_regulon_activity._box`, redrawn in the one box style.

    ONE BOX, NO POINTS (2026-09-16, the author's fifth reading): the boxes are
    cnsfig.boxes.draw_boxes (Figure 2 H/I's: RULE_PT black lines, 0.79 mm open
    fliers, solid R/NR faces as Figure 5 J/L print them) and the jittered
    samples the original drew on top are no longer drawn. Same numbers: the
    data go to ax.boxplot untouched; only the marks changed.
    """
    labels = [lab for lab, _ in series]
    data = [vals for _, vals in series]
    bxp = boxes.draw_boxes(ax, data, range(1, len(data) + 1), colors, width=0.55)
    ax.set_xlim(*boxes.box_xlim(range(1, len(data) + 1), 0.55))
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels([f"{lab}\nn = {len(v)}" for lab, v in zip(labels, data)])
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="both", width=0.6, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    if title:
        ax.set_title(title)
    boxes.assert_no_points(ax, bxp)


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
    _box(ax, series, colors, "", title)
    rich_ylabel(ax, "*NFKB1* regulon activity (AUCell)\nmonocytes/macrophages")

    base.fit(fig)
    if save:
        base.save(fig, FIG, "S10_E", "S10_E_nfkb_regulon_momac")
    return fig


def main():
    fr = frames()
    return base.run({
        # Axis 0 of the original is the block this panel keeps; the comparison
        # is made over that axis and says nothing about the one withdrawn.
        "S10_E": (lambda: draw_E(fr),
                 lambda: A._panel(fr["reg_df"], fr["fb"], fr["sv"]),
                 None, [0]),
    })


if __name__ == "__main__":
    sys.exit(main())
