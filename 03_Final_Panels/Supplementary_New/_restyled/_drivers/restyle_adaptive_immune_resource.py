#!/usr/bin/env python3
"""
Restyle S10D and S10E - adaptive immune resource (Reviewer 2, point R2.2).

Analysis: 04_Revision_Analyses/09_R2.2_Adaptive_Immune/scripts/adaptive_immune_resource.py

NOTHING IS RECOMPUTED HERE. The driver imports that module for its constants
and reads the three tables its `main()` already wrote. It never opens a h5ad and
it never reaches line 146, which is the `PREPARATION/"GSEA"` read that
SUPPLEMENTARY_AUDIT.md fault 1 raises and that is awaiting the author's ruling,
nor line 147, the `if not f.exists(): continue` that would skip a lineage in
silence. Restyling this panel must not, and does not, decide either question.
S10E draws exactly the table that is on disk today.

    conda run -n Liudeng_Python_310 python restyle_adaptive_immune_resource.py
    conda run -n Liudeng_Python_310 python restyle_adaptive_immune_resource.py --check
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _driver_base as base
import panel_style_cns as style

A = base.analysis("09_R2.2_Adaptive_Immune/scripts/adaptive_immune_resource.py")
OUT = base.outputs_of(A)
FIG = "S10_PreTx_and_Adaptive"

# Printed boxes, millimetres. See RESTYLE_REPORT.md for why S10E is re-laid out.
D_W, D_H = 171.0, 90.0
E_W, E_H = 171.0, 88.0


def frames():
    """The three tables the analysis already wrote. No computation."""
    return dict(
        fractions=pd.read_csv(base.require(
            OUT / "adaptive_state_fractions.csv", "S10D state fractions")),
        tests=pd.read_csv(base.require(
            OUT / "adaptive_state_tests.csv", "S10D state tests")),
        hallmark=pd.read_csv(base.require(
            OUT / "adaptive_hallmark_top.csv", "S10E Hallmark table")),
    )


# --------------------------------------------------------------------------
# S10D - minor state composition, 1 x 4 stacked bars. Layout unchanged.
# --------------------------------------------------------------------------
def draw_D(fr, save=True):
    style.apply()
    lineages = list(A.LINEAGES)
    fig, axes = style.subplots_mm(D_W, D_H, 1, len(lineages))
    for ax, lineage in zip(axes, lineages):
        sub = fr["fractions"][fr["fractions"]["lineage"] == lineage]
        states = sorted(sub["minor_cell_state"].unique())
        means = (sub.groupby(["group4", "minor_cell_state"], observed=True)["fraction"]
                 .mean().unstack(fill_value=0).reindex(A.GROUPS).fillna(0))
        means = means.reindex(columns=states, fill_value=0)
        bottom = np.zeros(len(A.GROUPS))
        cmap = plt.get_cmap("tab20")
        for k, state in enumerate(states):
            ax.bar(range(len(A.GROUPS)), means[state], bottom=bottom, width=0.7,
                   color=cmap(k % 20), edgecolor="white", linewidth=0.4,
                   label=state.replace("_", " "))
            bottom += means[state].values
        ax.set_xticks(range(len(A.GROUPS)))
        ax.set_xticklabels(A.GROUPS, rotation=45, ha="right")
        ax.set_title(lineage)
        ax.set_ylim(0, 1)
        ax.tick_params(axis="both", width=0.6, length=2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.legend(loc="upper left", bbox_to_anchor=(0, -0.22), ncol=1,
                  handlelength=1.0, handletextpad=0.4, labelspacing=0.25)
    axes[0].set_ylabel("Fraction of lineage")
    style.margins_mm(fig, left=13, right=2, top=6, bottom=42, wspace=0.34)
    if save:
        base.save(fig, FIG, "S10_D", "S10_D_adaptive_composition")
    return fig


# --------------------------------------------------------------------------
# S10E - top Hallmark programmes, post-treatment.
#
# LAYOUT CHANGED, content identical. Version A set four panels side by side and
# the term labels at 5.13 pt as printed. The longest label here is 29
# characters ("TNF-alpha Signaling via NF-kB"); at 7 pt that is about 36 mm of
# paper, and four of them plus their bars do not fit across a 171 mm page - the
# labels alone would need 144 mm. Two rows of two gives each cell 85 mm, which
# holds a 36 mm label and a 45 mm bar track. The same eight terms per lineage,
# the same NES values, the same order, the same colours.
# --------------------------------------------------------------------------
def _selection(hallmark, lineage):
    """Version A's selection, verbatim, so the same eight terms are drawn."""
    post = hallmark[hallmark["phase"] == "post"]
    sub = post[post["lineage"] == lineage].drop_duplicates("term")
    sub = pd.concat([sub[sub["direction"] == "NR-enriched"].head(4),
                     sub[sub["direction"] == "R-enriched"].head(4)])
    return sub.sort_values("nes")


def draw_E(fr, save=True):
    style.apply()
    lineages = list(A.LINEAGES)
    fig, axes = style.subplots_mm(E_W, E_H, 2, 2)
    for ax, lineage in zip(axes.ravel(), lineages):
        sub = _selection(fr["hallmark"], lineage)
        y = np.arange(len(sub))
        colors = ["#B2182B" if v > 0 else "#2166AC" for v in sub["nes"]]
        ax.barh(y, sub["nes"], color=colors, edgecolor="#444444", linewidth=0.4,
                height=0.7)
        ax.axvline(0, color="#666666", linewidth=0.5)
        ax.set_yticks(y)
        # Version A's truncation rule, kept exactly: Hallmark names run to 34
        # characters and "Interferon Gamma Response" is quoted by name in the
        # response letter, so the cut must not fall shorter than that.
        ax.set_yticklabels([t if len(t) < 34 else t[:31] + "..."
                            for t in sub["term"]])
        ax.set_title(lineage)
        ax.tick_params(axis="x", width=0.6, length=2)
        ax.tick_params(axis="y", length=0)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
    # supxlabel takes matplotlib's `figure.labelsize`, which is "large" - 1.2x
    # font.size, so 9.6 pt - and cnsplots has no opinion about it. Set it to the
    # body size so this one string does not sit 1.6 pt above every other axis
    # label in the set.
    fig.supxlabel("NES, post-treatment (positive = enriched in non-responders)",
                  fontsize=style.body_pt())
    style.margins_mm(fig, left=38, right=2, top=5, bottom=10,
                     wspace=0.62, hspace=0.55)
    if save:
        base.save(fig, FIG, "S10_E", "S10_E_adaptive_hallmark")
    return fig


def main():
    fr = frames()
    return base.run({
        "S10_D": (lambda: draw_D(fr),
                  lambda: A._panel_composition(fr["fractions"], fr["tests"])),
        "S10_E": (lambda: draw_E(fr),
                  lambda: A._panel_hallmark(fr["hallmark"])),
    })


if __name__ == "__main__":
    sys.exit(main())
