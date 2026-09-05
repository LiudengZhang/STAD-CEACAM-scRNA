#!/usr/bin/env python3
"""
Restyle S10A, S10B and S10C - is the inflammatory programme already present
before treatment? (Reviewer 2, point R2.1).

Analysis: 04_Revision_Analyses/08_R2.1_PreTx_Inflammatory/scripts/pretreatment_inflammatory.py

NOTHING IS RECOMPUTED HERE. The driver never opens MoMac.h5ad, so
`rank_genes_groups` and `score_genes` are not re-run and the IL-1b+ signature is
not re-derived; and it never walks `RECOMPUTE_GSEA.glob(...)` at line 136, so
the module-12 Hallmark tables are not re-read and no cell type can be added to
or dropped from S10C by a file appearing or disappearing on disk. It reads the
three tables `main()` already wrote:

    outputs/pretx_state_abundance.csv   -> frac, for S10A
    outputs/pretx_signature_scores.csv  -> score, for S10B
    outputs/nfkb_pre_vs_post.csv        -> the pivoted NES/FDR table for S10C

S10C is drawn from `n`, not from that table as read: `main()` line 198 drops the
cell types missing either timepoint and sorts by post-treatment NES before
handing it to `_panel_nfkb_shift`. That one line is copied verbatim into
`frames()` and the same object is given to both the restyled panel and the
original, so `--check` compares drawing against drawing.

`_panel_four_group` is called twice by `main()`, with `frac`/"fraction"/pct=True
and with `score`/"score", so there is one restyled implementation here and two
build functions over it. The ylabels, the `sub` directory names, the file stems
and the `pct` flag are exactly the ones `main()` passes.

S10A, S10B and S10C sit side by side on one row of the assembled figure, so all
three are drawn at the same height; their three widths tile the 171 mm page.

    conda run -n Liudeng_Python_310 python restyle_pretreatment_inflammatory.py --check
    conda run -n Liudeng_Python_310 python restyle_pretreatment_inflammatory.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _driver_base as base
import panel_style_cns as style

A = base.analysis("08_R2.1_PreTx_Inflammatory/scripts/pretreatment_inflammatory.py")
OUT = base.outputs_of(A)
FIG = "S10_PreTx_and_Adaptive"

# Printed boxes, millimetres. Version A: S10A 54 x 44, S10B 54 x 44,
# S10C 51 x 46. The three share a row of the assembled figure, so they share a
# height, and the three widths tile the 171 mm page exactly.
#
# The height is set by S10C, which stacks thirteen cell-type labels at 7 pt and
# carries a legend under a two-line x label. The widths are set by S10C too:
# its x label, "NES, TNFa/NF-kB / (positive = enriched in non-responders)",
# measures 49.2 mm at 8 pt and is centred on the axes, whose centre sits at
# 7 + W/2 once the 17 mm label gutter is reserved, so the panel needs
# W >= 63.2 mm or the label is clipped on the right (measured at 57 mm: over by
# 3.1 mm). 65 mm for S10C leaves 53 mm each for S10A and S10B.
AB_W = 53.0
C_W = 65.0
ROW_H = 60.0

# Version A, S10A/S10B: 5.5 x 4.2 cm at SCALE 4, fitted into 54 x 44 mm at
# 0.2455. Smallest type, the bracket P values, was 5 * SCALE, printed at
# 5 * 4 * 0.2455 = 4.91 pt.
AB_FIT = 0.2455
AB_TYPE = 5.0 * 4 * AB_FIT
AB_MARK = (style.tick_pt() / AB_TYPE) * AB_FIT
AB_AREA = AB_MARK ** 2

# Version A, S10C: 6.5 x 4.6 cm at SCALE 4, fitted into 51 x 46 mm at 0.1962.
# Smallest type, the cell-type labels, was 5.5 * SCALE, printed at
# 5.5 * 4 * 0.1962 = 4.316 pt.
C_FIT = 0.1962
C_TYPE = 5.5 * 4 * C_FIT
C_MARK = (style.tick_pt() / C_TYPE) * C_FIT
C_AREA = C_MARK ** 2


def frames():
    """The three tables the analysis already wrote. No computation."""
    nfkb = pd.read_csv(base.require(
        OUT / "nfkb_pre_vs_post.csv", "S10C NF-kB pre versus post"))
    # main() line 198, verbatim: the frame `_panel_nfkb_shift` is given.
    n = nfkb.dropna(subset=["pre_nes", "post_nes"]).sort_values(
        "post_nes", ascending=False)
    return dict(
        frac=pd.read_csv(base.require(
            OUT / "pretx_state_abundance.csv", "S10A IL-1b+ state abundance")),
        score=pd.read_csv(base.require(
            OUT / "pretx_signature_scores.csv", "S10B IL-1b+ signature scores")),
        n=n,
    )


# --------------------------------------------------------------------------
# S10A and S10B - the four-group boxplot. One implementation, two callers,
# exactly as in the analysis. Layout unchanged.
# --------------------------------------------------------------------------
def _four_group(data, col, ylabel, panel, stem_name, pct=False, save=True):
    style.apply()
    fig, ax = style.subplots_mm(AB_W, ROW_H)
    vals = [data.loc[data["group4"] == g, col].values * (100 if pct else 1)
            for g in A.GROUPS]
    bp = ax.boxplot(vals, positions=range(4), widths=0.6, patch_artist=True,
                    showfliers=False,
                    boxprops=dict(linewidth=0.8 * AB_MARK),
                    whiskerprops=dict(linewidth=0.8 * AB_MARK),
                    capprops=dict(linewidth=0.8 * AB_MARK),
                    medianprops=dict(color="black", linewidth=1.2 * AB_MARK))
    for patch, g in zip(bp["boxes"], A.GROUPS):
        patch.set_facecolor(A.COLORS[g])
        patch.set_edgecolor("#444444")
    rng = np.random.default_rng(3)
    for i, v in enumerate(vals):
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=9 * 4 * AB_AREA,
                   c="#333333", zorder=3, alpha=0.85,
                   edgecolors="white", linewidths=0.3 * 4 * AB_MARK)
    top = max(v.max() for v in vals)
    bot = min(v.min() for v in vals)
    span = top - bot
    for k, (i, j) in enumerate(((0, 1), (2, 3))):
        p, _ = A.mw(vals[i], vals[j])
        yy = top + span * (0.10 + 0.16 * k)
        ax.plot([i, i, j, j], [yy, yy + span * 0.04, yy + span * 0.04, yy],
                color="#444444", linewidth=0.8 * AB_MARK)
        ax.text((i + j) / 2, yy + span * 0.05, f"P = {p:.3f}", ha="center",
                va="bottom", fontsize=style.tick_pt())
    ax.set_ylim(bot - span * 0.10, top + span * 0.42)
    ax.set_xticks(range(4))
    ax.set_xticklabels(A.GROUPS, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="both", width=0.6, length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    style.margins_mm(fig, left=15, right=3, top=4, bottom=16)
    if save:
        base.save(fig, FIG, panel, stem_name)
    return fig


def draw_A(fr, save=True):
    return _four_group(fr["frac"], "fraction",
                       "IL-1$\\beta$+ state\n(% of mono/macrophages)",
                       "S10_A", "S10_A_il1b_state_four_groups",
                       pct=True, save=save)


def draw_B(fr, save=True):
    return _four_group(fr["score"], "score",
                       "IL-1$\\beta$+ signature score",
                       "S10_B", "S10_B_il1b_signature_four_groups",
                       save=save)


# --------------------------------------------------------------------------
# S10C - pre versus post NES per cell type, as a dumbbell. Layout unchanged.
# --------------------------------------------------------------------------
def draw_C(fr, save=True):
    style.apply()
    n = fr["n"]
    fig, ax = style.subplots_mm(C_W, ROW_H)
    y = np.arange(len(n))
    for i, (_, r) in enumerate(n.iterrows()):
        ax.plot([r["pre_nes"], r["post_nes"]], [i, i], color="#bbbbbb",
                linewidth=1.0 * C_MARK, zorder=1)
        ax.scatter(r["pre_nes"], i, s=24 * C_AREA, c="#a2d2ff",
                   edgecolors="#444444", linewidths=0.4 * C_MARK, zorder=3)
        ax.scatter(r["post_nes"], i, s=24 * C_AREA, c="#f1c0e8",
                   edgecolors="#444444", linewidths=0.4 * C_MARK, zorder=3)
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([A.LABELS[c] for c in n["cell_type"]])
    ax.set_xlabel("NES, TNF$\\alpha$/NF-$\\kappa$B\n"
                  "(positive = enriched in non-responders)")
    ax.tick_params(axis="x", width=0.6, length=2)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    handles = [plt.Line2D([], [], marker="o", linestyle="none",
                          markersize=4 * C_MARK, markerfacecolor=c,
                          markeredgecolor="#444444", label=l)
               for c, l in (("#a2d2ff", "Pre-treatment"),
                            ("#f1c0e8", "Post-treatment"))]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.30),
              ncol=2, frameon=False, handlelength=1.0, handletextpad=0.4,
              columnspacing=1.0)
    style.margins_mm(fig, left=17, right=3, top=3, bottom=18)
    if save:
        base.save(fig, FIG, "S10_C", "S10_C_nfkb_pre_vs_post")
    return fig


def main():
    fr = frames()
    return base.run({
        "S10_A": (lambda: draw_A(fr),
                  lambda: A._panel_four_group(
                      fr["frac"], "fraction",
                      "IL-1$\\beta$+ state\n(% of mono/macrophages)",
                      "S10_A", "S10_A_il1b_state_four_groups", pct=True)),
        "S10_B": (lambda: draw_B(fr),
                  lambda: A._panel_four_group(
                      fr["score"], "score", "IL-1$\\beta$+ signature score",
                      "S10_B", "S10_B_il1b_signature_four_groups")),
        "S10_C": (lambda: draw_C(fr), lambda: A._panel_nfkb_shift(fr["n"])),
    })


if __name__ == "__main__":
    sys.exit(main())
