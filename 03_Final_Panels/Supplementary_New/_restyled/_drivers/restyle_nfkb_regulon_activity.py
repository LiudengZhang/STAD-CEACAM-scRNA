#!/usr/bin/env python3
"""
Restyle S11C - NF-kB regulon activity and the feedback programme
(Reviewer 1, point R1.8, affirmative evidence).

Analysis: 02_New_Analyses/07_R1.8_NFkB_Specificity/scripts/nfkb_regulon_activity.py

NOTHING IS RECOMPUTED HERE. The driver does not call `load_regulons()` and so
never opens `02_Preparation_for_Panels/SCENIC/aucell_matrix.csv`. That matters
beyond speed: the deposited AUCell matrix is a 12k-cell subset, so a driver that
re-derived sample-level means from it could quietly draw a different panel from
the one `main()` drew. It also never opens a h5ad, so `score_genes` is not
re-run and the feedback-target scores are not re-derived.

Instead it reads the three tables `main()` already wrote and hands them to the
same panel code:

    outputs/nfkb_regulon_activity.csv   -> reg_df, the regulon summary; the
                                          panel takes only the post-treatment
                                          NFKB1(+)/MoMac P value from it
    outputs/nfkb_feedback_targets.csv   -> fb, which `_panel` accepts and does
                                          not use; passed for the check anyway
    outputs/nfkb_sample_values.csv      -> sv, the per-sample values that are
                                          what the boxes and points actually
                                          plot

`_box` is copied in below as a restyled local helper. It is drawing code only -
it computes nothing except the jitter, from its own seeded generator - and its
behaviour is unchanged.

    conda run -n Liudeng_Python_310 python restyle_nfkb_regulon_activity.py --check
    conda run -n Liudeng_Python_310 python restyle_nfkb_regulon_activity.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _driver_base as base
import panel_style_cns as style

A = base.analysis("07_R1.8_NFkB_Specificity/scripts/nfkb_regulon_activity.py")
OUT = base.outputs_of(A)
FIG = "S11_Affirmative_Analyses"

# Printed box, millimetres. Version A: 171 x 44 at fit 0.2391.
W, H = 171.0, 52.0

# Version A: 9.0 x 4.6 cm at SCALE 4, fitted into 171 x 44 mm at 0.2391. Its
# smallest type, `_box`'s x tick labels, was set at 5.0 * SCALE and printed at
# 5.0 * 4 * 0.2391 = 4.782 pt. Every mark and stroke below is rescaled by the
# factor that keeps its size relative to the type unchanged.
A_FIT = 0.2391
A_TYPE = 5.0 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2


def frames():
    """The three tables the analysis already wrote. No computation."""
    return dict(
        reg_df=pd.read_csv(base.require(
            OUT / "nfkb_regulon_activity.csv", "S11C regulon activity")),
        fb=pd.read_csv(base.require(
            OUT / "nfkb_feedback_targets.csv", "S11C feedback targets")),
        sv=pd.read_csv(base.require(
            OUT / "nfkb_sample_values.csv", "S11C per-sample values")),
    )


def _box(ax, series, colors, ylabel, title=None):
    """`nfkb_regulon_activity._box`, restyled. Same marks, same numbers.

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


def draw_C(fr, save=True):
    style.apply()
    reg_df, sv = fr["reg_df"], fr["sv"]
    fig, axes = style.subplots_mm(W, H, 1, 2)

    # left: NFKB1 regulon in monocytes/macrophages, both timepoints
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
    _box(axes[0], series, colors,
         "NFKB1 regulon activity (AUCell)\nmonocytes/macrophages", title)

    # right: NF-kB negative-feedback target score after treatment, by compartment
    pan = sv[(sv["source"] == "panel")
             & (sv["key"] == "NF-kB negative-feedback targets")
             & (sv["phase"] == "Post")]
    order = [c for c in ("Monocytes/Macrophages", "Epithelial", "Fibroblast")
             if c in set(pan["stratum"])]
    series, colors = [], []
    for comp in order:
        short = "Mono/Mac" if comp == "Monocytes/Macrophages" else comp
        for grp, c in (("R", A.COLOR_R), ("NR", A.COLOR_NR)):
            v = pan.loc[(pan["stratum"] == comp) & (pan["group"] == grp),
                        "value"].values
            if len(v):
                series.append((f"{short}\n{grp}", v))
                colors.append(c)
    _box(axes[1], series, colors,
         "NF-$\\kappa$B feedback target score\n(post-treatment)")
    axes[1].axhline(0, color="#666666", linewidth=0.5, zorder=0)

    style.margins_mm(fig, left=15, right=3, top=6, bottom=14, wspace=0.30)
    if save:
        base.save(fig, FIG, "S11_C", "S11_C_nfkb_regulon_and_feedback")
    return fig


def main():
    fr = frames()
    return base.run({
        "S11_C": (lambda: draw_C(fr),
                  lambda: A._panel(fr["reg_df"], fr["fb"], fr["sv"])),
    })


if __name__ == "__main__":
    sys.exit(main())
