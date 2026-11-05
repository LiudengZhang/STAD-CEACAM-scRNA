#!/usr/bin/env python3
"""Supplementary Figure S4 A: S-MP4 and S-MP5 scores across the four
response/phase groups, with two-sided Mann-Whitney P for four contrasts.

Reads the stored per-sample scores (mp4_permutation_results.json) and
outputs/mp_group_comparisons.csv of
04_Revision_Analyses/03_R1.4_MP_Direction_PrePost/scripts/mp_direction_and_prepost.py;
writes S4_CEACAM5_vs_CEACAM6/S4_A. `_draw` is the analysis's own drawing code,
used only by --check.

Run: python draw_mp_direction_and_prepost.py [--check]"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _driver_base as base
import panel_style_cns as style
from cnsfig import boxes
from shared.figure_config import GROUP_COLORS

base.apply_style()

A = base.analysis("03_R1.4_MP_Direction_PrePost/scripts/mp_direction_and_prepost.py")
OUT = base.outputs_of(A)
FIG = "S4_CEACAM5_vs_CEACAM6"

# Panel size in mm.
W, H = 171.0, 60.0

PROGRAMS = ("S-MP4", "S-MP5")
LADDER = 1.5          # bracket-to-bracket step, in units of 0.075 * top

def frames():
    """The stored per-group values and the comparison table."""
    mp_json = base.require(
        A.PREPARATION / "Metaprogram_Permutation" / "mp4_permutation_results.json",
        "S7A stored metaprogram permutation results")
    mp = json.loads(mp_json.read_text())
    vals = {prog: {g: np.asarray(mp[prog][f"{g}_values"], float)
                   for g in A.GROUPS}
            for prog in PROGRAMS}
    df = pd.read_csv(base.require(OUT / "mp_group_comparisons.csv",
                                  "S7A group comparisons"))
    return dict(vals=vals, df=df)

def _draw(fr, axes, scale, fontsize):
    """The analysis's own drawing code, used by `_original()`."""
    vals, df = fr["vals"], fr["df"]
    for ax, prog in zip(axes, PROGRAMS):
        v = vals[prog]
        data = [v[g] for g in A.GROUPS]
        bp = ax.boxplot(data, positions=range(4), widths=0.6, patch_artist=True,
                        showfliers=False,
                        boxprops=dict(linewidth=0.8 * scale),
                        whiskerprops=dict(linewidth=0.8 * scale),
                        capprops=dict(linewidth=0.8 * scale),
                        medianprops=dict(color="black", linewidth=1.2 * scale))
        for patch, g in zip(bp["boxes"], A.GROUPS):
            patch.set_facecolor(GROUP_COLORS[g])
            patch.set_edgecolor("#444444")
        rng = np.random.default_rng(0)
        for i, g in enumerate(A.GROUPS):
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(v[g])), v[g],
                       s=8 * 4 * scale ** 2, c="#333333", zorder=3, alpha=0.8,
                       edgecolors="white", linewidths=0.3 * 4 * scale)

        top = max(x.max() for x in data)
        step = 0.075 * top
        y = top + step
        for (ga, gb, _), _ in zip(A.CONTRASTS, range(4)):
            ia, ib = A.GROUPS.index(ga), A.GROUPS.index(gb)
            p = df[(df["program"] == prog) & (df["group_a"] == ga)
                   & (df["group_b"] == gb)]["p_two_tailed"].iloc[0]
            lo, hi = sorted((ia, ib))
            ax.plot([lo, lo, hi, hi], [y, y + step * 0.25, y + step * 0.25, y],
                    color="#444444", linewidth=0.8 * scale)
            p_str, p_kw = style.p_text_kw(p)
            ax.text((lo + hi) / 2, y + step * 0.3, p_str,
                    ha="center", va="bottom", **p_kw)
            y += step * 1.15

        ax.set_xticks(range(4))
        ax.set_xticklabels(A.GROUPS)
        ax.set_ylabel(f"{prog} score")
        ax.set_title(prog)
        ax.set_ylim(0, y + step)
        base.finish(ax)

def _draw_boxes(fr, axes):
    """Box plots without points; contrast brackets stacked in the analysis's
    order."""
    vals, df = fr["vals"], fr["df"]
    for ax, prog in zip(axes, PROGRAMS):
        v = vals[prog]
        data = [v[g] for g in A.GROUPS]
        bxp = boxes.draw_boxes(ax, data, range(4),
                               [GROUP_COLORS[g] for g in A.GROUPS], width=0.6)
        top = max(x.max() for x in data)
        rng = top                     # brackets step in fractions of the top value
        text = None
        for k, ((ga, gb, _), _) in enumerate(zip(A.CONTRASTS, range(4))):
            ia, ib = A.GROUPS.index(ga), A.GROUPS.index(gb)
            p = df[(df["program"] == prog) & (df["group_a"] == ga)
                   & (df["group_b"] == gb)]["p_two_tailed"].iloc[0]
            lo, hi = sorted((ia, ib))
            # First bracket 0.075 * top above the tallest box, each next one
            # LADDER steps higher.
            _, text, _ = boxes.bracket(ax.figure, ax, lo, hi, top, rng, p,
                                       kind="pair",
                                       lift=0.075 * (1 + LADDER * k),
                                       arm=0.075 * 0.25)
        ax.set_xticks(range(4))
        ax.set_xticklabels(A.GROUPS)
        ax.set_ylabel(f"{prog} score")
        ax.set_title(prog)
        ax.set_xlim(*boxes.box_xlim(range(4), 0.6))
        ax.set_ylim(bottom=0)
        boxes.ylim_above(ax, text, pad_mm=0.6)
        base.finish(ax)
        boxes.assert_no_points(ax, bxp)

def draw_A(fr, save=True):
    base.apply_style()
    fig, axes = style.subplots_mm(W, H, 1, 2)
    _draw_boxes(fr, axes)
    base.fit(fig, wspace=0.30)
    if save:
        base.save(fig, FIG, "S4_A", "S4_A_metaprogram_four_groups")
    return fig

def _original(fr):
    """The analysis's own panel, with its sizes and rcParams, for --check."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 2, figsize=(9.0 * A.SCALE * A.CM,
                                            4.6 * A.SCALE * A.CM))
    _draw(fr, axes, 1.0, 4.5 * A.SCALE)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.10, wspace=0.28)
    fig.savefig("/dev/null")          # intercepted by --check
    return fig

def main():
    fr = frames()
    return base.run({
        "S4_A": (lambda: draw_A(fr), lambda: _original(fr), base.aliases),
    })

if __name__ == "__main__":
    sys.exit(main())
