#!/usr/bin/env python3
"""
Restyle S8A - metaprogram direction and pre/post contrasts (Reviewer 1, R1.4).

Analysis: 02_New_Analyses/03_R1.4_MP_Direction_PrePost/scripts/mp_direction_and_prepost.py

This is the one analysis in the set whose panel is drawn inline in `main()`
rather than in a `_panel*` function, so there was no plotting entry point to
call. Separating it needed no edit to the analysis all the same, because
neither of the two things the panel draws is computed at draw time:

  `vals`  is read verbatim out of a stored JSON,
          02_Preparation_for_Panels/Metaprogram_Permutation/mp4_permutation_results.json,
          with the same two lines `main()` uses. No test is re-run.
  `df`    is outputs/mp_group_comparisons.csv, which `main()` wrote. The exact
          permutation test and the four Mann-Whitney contrasts are NOT
          recomputed here; their P values are read.

So `--check` compares this drawing against the analysis's own inline drawing by
running that inline block through a copy of the loop, which is why the check for
this panel is written differently from the others: there is no `A._panel*` to
call. See `_original()`.

    conda run -n Liudeng_Python_310 python restyle_mp_direction_and_prepost.py --check
    conda run -n Liudeng_Python_310 python restyle_mp_direction_and_prepost.py
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

A = base.analysis("03_R1.4_MP_Direction_PrePost/scripts/mp_direction_and_prepost.py")
OUT = base.outputs_of(A)
FIG = "S8_CEACAM_Metaprogram"

W, H = 171.0, 66.0

# Version A: 9.0 x 4.6 cm at SCALE 4, fitted into 171 x 52 mm at 0.2826.
# Its smallest type, the bracket P values, was set at 4.5 * SCALE and printed
# at 4.5 * 4 * 0.2826 = 5.087 pt.
A_FIT = 0.2826
A_TYPE = 4.5 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2

PROGRAMS = ("S-MP4", "S-MP5")


def frames():
    """The stored values and the comparison table. Nothing is computed."""
    mp_json = base.require(
        A.PREPARATION / "Metaprogram_Permutation" / "mp4_permutation_results.json",
        "S8A stored metaprogram permutation results")
    mp = json.loads(mp_json.read_text())
    vals = {prog: {g: np.asarray(mp[prog][f"{g}_values"], float)
                   for g in A.GROUPS}
            for prog in PROGRAMS}
    df = pd.read_csv(base.require(OUT / "mp_group_comparisons.csv",
                                  "S8A group comparisons"))
    return dict(vals=vals, df=df)


def _draw(fr, axes, scale, fontsize):
    """The panel, drawn once. `scale` is 1 for the restyle, SCALE for Version A.

    This is the analysis's own inline drawing code. It is parameterised only in
    the two places the restyle touches - mark sizes and type - so that calling
    it with Version A's numbers reproduces Version A exactly and calling it with
    the restyle's numbers reproduces the restyle. Every value plotted is the
    same in both.
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
            ax.text((lo + hi) / 2, y + step * 0.3, f"P = {p:.3f}",
                    ha="center", va="bottom", fontsize=fontsize)
            y += step * 1.15

        ax.set_xticks(range(4))
        ax.set_xticklabels(A.GROUPS)
        ax.set_ylabel(f"{prog} score")
        ax.set_title(prog)
        ax.set_ylim(0, y + step)
        ax.tick_params(axis="both", width=0.6, length=2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)


def draw_A(fr, save=True):
    style.apply()
    fig, axes = style.subplots_mm(W, H, 1, 2)
    _draw(fr, axes, MARK, style.tick_pt())
    style.margins_mm(fig, left=14, right=3, top=6, bottom=11, wspace=0.30)
    if save:
        base.save(fig, FIG, "S8_A", "S8_A_metaprogram_four_groups")
    return fig


def _original(fr):
    """Version A's panel, for the check. Same code, Version A's sizes."""
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
        "S8_A": (lambda: draw_A(fr), lambda: _original(fr)),
    })


if __name__ == "__main__":
    sys.exit(main())
