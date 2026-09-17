#!/usr/bin/env python3
"""
Draw S7A - metaprogram direction and pre/post contrasts (Reviewer 1, R1.4).

Analysis: 04_Revision_Analyses/03_R1.4_MP_Direction_PrePost/scripts/mp_direction_and_prepost.py

NOTHING IS RECOMPUTED HERE. This is the one analysis in the set whose panel is
drawn inline in `main()` rather than in a `_panel*` function, so there is no
plotting entry point to call. Separating it needs no edit to the analysis all
the same, because neither of the two things the panel draws is computed at draw
time:

  `vals`  is read verbatim out of a stored JSON of permutation results, with
          the same two lines `main()` uses. No test is re-run.
  `df`    is outputs/mp_group_comparisons.csv, which `main()` wrote. The exact
          permutation test and the four Mann-Whitney contrasts are NOT
          recomputed here; their P values are read.

So `--check` compares this drawing against the analysis's own inline drawing by
running that inline block through a copy of the loop, which is why the check for
this panel is written differently from the others: there is no `A._panel*` to
call. See `_original()`.

ONE BOX, NO POINTS (2026-09-16, the author's fifth reading). The redraw
(`_draw_boxes`) sets its boxes with `cnsfig.boxes.draw_boxes` - the Figure 2
H/I box: RULE_PT black lines, 0.79 mm open fliers, faces from the analysis's
four group colours - and its three stacked pairwise brackets with
`cnsfig.boxes.bracket`, each label's ink 0.4 mm above its own line, in black;
the jittered individual samples are no longer drawn ("box plots should look
alike throughout; don't show every point"). `_draw` is kept as it was, because
it is the analysis's own inline drawing code and `_original()` still runs it
for the check; the comparator then sees the removed point collections, the
added flier markers and the moved brackets, and not one box statistic moved.

    python draw_mp_direction_and_prepost.py --check
    python draw_mp_direction_and_prepost.py
"""

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

base.apply_style()

A = base.analysis("03_R1.4_MP_Direction_PrePost/scripts/mp_direction_and_prepost.py")
OUT = base.outputs_of(A)
FIG = "S8_CEACAM5_vs_CEACAM6"

# The printed box, millimetres. Two boxplot panels side by side across the
# 183 mm page less its 6 mm margins; the depth is what four annotated contrast
# brackets above the tallest box need at this type size.
W, H = 171.0, 60.0

# The panel as previously laid out: 9.0 x 4.6 cm drawn at four times print size
# and fitted into a 171 x 52 mm box at 0.2826. Its smallest type, the bracket P
# values, was set at 4.5 x that scale and printed at 5.087 pt. Marks are
# rescaled by the factor that keeps their size relative to the type unchanged.
A_FIT = 0.2826
A_TYPE = 4.5 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2

PROGRAMS = ("S-MP4", "S-MP5")
LADDER = 1.5          # bracket-to-bracket step, in units of 0.075 * top (see _draw_boxes)


def frames():
    """The stored values and the comparison table. Nothing is computed."""
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
    """The panel, drawn once. `scale` is 1 here, the draw scale for the original.

    This is the analysis's own inline drawing code. It is parameterised only in
    the two places the redraw touches - mark sizes and type - so that calling it
    with the original's numbers reproduces the original and calling it with
    these reproduces this panel. Every value plotted is the same in both.
    """
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
            patch.set_facecolor(A.COLORS[g])
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
            # P through panel_style_cns.p_text_kw (2026-09-15): two decimals,
            # a star below 0.05, as every main-figure bracket prints it. The
            # analysis's three-decimal strings are declared as renames in
            # labels.py RENAMES_SUPPLEMENTARY.
            p_str, p_kw = style.p_text_kw(p)
            ax.text((lo + hi) / 2, y + step * 0.3, p_str,
                    ha="center", va="bottom", **p_kw)
            y += step * 1.15

        ax.set_xticks(range(4))
        ax.set_xticklabels(A.GROUPS)
        ax.set_ylabel(f"{prog} score")
        ax.set_title(prog)
        ax.set_ylim(0, y + step)
        ax.tick_params(axis="both", width=0.6, length=2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)


def _draw_boxes(fr, axes):
    """The redraw: the same data and contrasts as `_draw`, in the one box style.

    Boxes through cnsfig.boxes.draw_boxes, brackets through
    cnsfig.boxes.bracket stacked in the analysis's contrast order, each one a
    step higher (`lift`), no individual points.
    """
    vals, df = fr["vals"], fr["df"]
    for ax, prog in zip(axes, PROGRAMS):
        v = vals[prog]
        data = [v[g] for g in A.GROUPS]
        bxp = boxes.draw_boxes(ax, data, range(4),
                               [A.COLORS[g] for g in A.GROUPS], width=0.6)
        top = max(x.max() for x in data)
        rng = top                     # the analysis steps its brackets by the top value
        text = None
        for k, ((ga, gb, _), _) in enumerate(zip(A.CONTRASTS, range(4))):
            ia, ib = A.GROUPS.index(ga), A.GROUPS.index(gb)
            p = df[(df["program"] == prog) & (df["group_a"] == ga)
                   & (df["group_b"] == gb)]["p_two_tailed"].iloc[0]
            lo, hi = sorted((ia, ib))
            # The analysis's stacking: the first bracket's foot one step
            # (0.075 * top) above the tallest box, arms a quarter of a step.
            # Each next bracket is LADDER steps higher: the analysis used
            # 1.15, which at 6 pt left a label 0.46 mm under the arm foot of
            # the bracket above it (measured 2026-09-16); the round holds any
            # label >= 0.5 mm clear of any arm (sweep_pages.check_bracket_
            # clearance), so the ladder is 1.5 steps.
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
        ax.tick_params(axis="both", width=0.6, length=2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        boxes.assert_no_points(ax, bxp)


def draw_A(fr, save=True):
    base.apply_style()
    fig, axes = style.subplots_mm(W, H, 1, 2)
    _draw_boxes(fr, axes)
    base.fit(fig, wspace=0.30)
    if save:
        base.save(fig, FIG, "S8_A", "S8_A_metaprogram_four_groups")
    return fig


def _original(fr):
    """The analysis's own panel, for the check. Same code, its own sizes."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 2, figsize=(9.0 * A.SCALE * A.CM,
                                            4.6 * A.SCALE * A.CM))
    _draw(fr, axes, 1.0, 4.5 * A.SCALE)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.10, wspace=0.28)
    fig.savefig("/dev/null")          # intercepted by the capture harness
    return fig


def main():
    fr = frames()
    return base.run({
        "S8_A": (lambda: draw_A(fr), lambda: _original(fr), base.aliases),
    })


if __name__ == "__main__":
    sys.exit(main())
