#!/usr/bin/env python3
"""Draw Supplementary Figure S9 E: NFKB1 regulon activity (AUCell) in
monocytes/macrophages, per sample, by phase and response.

Reads nfkb_sample_values.csv and nfkb_regulon_activity.csv written by
04_Revision_Analyses/07_R1.8_NFkB_Specificity/scripts/nfkb_regulon_activity.py
--check compares the drawing with the left axes of that script's panel.


    python draw_nfkb_regulon_activity.py [--check]"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

import _driver_base as base
import panel_style_cns as style
from cnsfig import boxes
from cnsfig.rich import rich_ylabel
from shared.figure_config import GROUP_COLORS

base.apply_style()

A = base.analysis("07_R1.8_NFkB_Specificity/scripts/nfkb_regulon_activity.py")
OUT = base.outputs_of(A)
FIG = "S9_MoMac_Identity_NFkB"

# Panel size, mm.
W, H = 37.0, 60.0

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
    """Sample-level boxes without individual points."""
    labels = [lab for lab, _ in series]
    data = [vals for _, vals in series]
    bxp = boxes.draw_boxes(ax, data, range(1, len(data) + 1), colors, width=0.55)
    ax.set_xlim(*boxes.box_xlim(range(1, len(data) + 1), 0.55))
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(list(labels))
    ax.set_ylabel(ylabel)
    base.finish(ax)
    if title:
        ax.set_title(title)
    boxes.assert_no_points(ax, bxp)

def draw_E(fr, save=True):
    base.apply_style()
    reg_df, sv = fr["reg_df"], fr["sv"]
    fig, ax = style.subplots_mm(W, H)

    # NFKB1(+) regulon in the MoMac stratum.
    reg = sv[(sv["source"] == "regulon") & (sv["key"] == "NFKB1(+)")
             & (sv["stratum"] == "MoMac")]
    series, colors = [], []
    for phase in ("Pre", "Post"):
        for grp in ("R", "NR"):
            c = GROUP_COLORS[f"{phase}-{grp}"]
            v = reg.loc[(reg["phase"] == phase) & (reg["group"] == grp),
                        "value"].values
            if len(v):
                series.append((f"{phase}-{grp}", v))
                colors.append(c)
    post = reg_df[(reg_df["regulon"] == "NFKB1(+)")
                  & (reg_df["cell_type"] == "MoMac")
                  & (reg_df["phase"] == "Post")]
    # Post-treatment two-tailed P as the title.
    title = (f"{style.p_exact(post['p_two_tailed'].iloc[0])}"
             if len(post) else None)
    _box(ax, series, colors, "", title)
    rich_ylabel(ax, "*NFKB1* regulon activity (AUCell)\nmonocytes/macrophages")

    base.fit(fig)
    if save:
        base.save(fig, FIG, "S9_E", "S9_E_nfkb_regulon_momac")
    return fig

def main():
    fr = frames()
    return base.run({
        # --check compares against axis 0 of the analysis panel.
        "S9_E": (lambda: draw_E(fr),
                 lambda: A._panel(fr["reg_df"], fr["fb"], fr["sv"]),
                 None, [0]),
    })

if __name__ == "__main__":
    sys.exit(main())
