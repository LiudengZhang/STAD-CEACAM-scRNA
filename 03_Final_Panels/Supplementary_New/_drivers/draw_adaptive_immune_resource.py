#!/usr/bin/env python3
"""
Draw S3E - adaptive immune minor-state composition (Reviewer 2, point R2.2).

Analysis: 04_Revision_Analyses/09_R2.2_Adaptive_Immune/scripts/adaptive_immune_resource.py

NOTHING IS RECOMPUTED HERE. The driver imports that module for its constants
and reads the two tables its `main()` already wrote. It never opens a h5ad and
it never reaches the Hallmark read that only the withdrawn companion panel
needed, nor the `if not f.exists(): continue` beside it that would skip a
lineage in silence.

Until 2026-09-15 (evening) this panel was spliced below the submitted S3
page by Supplementary_Fixes/patch_S3_add_adaptive.py and so was written into
Supplementary_Fixes/S3_E/. Since the S1-S6 redraw it is a panel of the
assembled figure (S3 then, S4 since the renumbering of 2026-09-16) like the other four, written into S4_CD8_TCells/S4_E/.

    python draw_adaptive_immune_resource.py --check
    python draw_adaptive_immune_resource.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import _driver_base as base
import panel_style_cns as style

base.apply_style()

A = base.analysis("09_R2.2_Adaptive_Immune/scripts/adaptive_immune_resource.py")
OUT = base.outputs_of(A)

# The printed box, millimetres. Full width across the 183 mm S3 page less its
# 6 mm margins, and deep enough for four stacked-bar cells whose minor-state
# legends run below them.
W, H = 171.0, 90.0


def frames():
    """The two tables the analysis already wrote. No computation."""
    return dict(
        fractions=pd.read_csv(base.require(
            OUT / "adaptive_state_fractions.csv", "S3E state fractions")),
        tests=pd.read_csv(base.require(
            OUT / "adaptive_state_tests.csv", "S3E state tests")),
    )


# --------------------------------------------------------------------------
# S3E - minor state composition, 1 x 4 stacked bars.
# --------------------------------------------------------------------------
def draw_E(fr, save=True):
    base.apply_style()
    lineages = list(A.LINEAGES)
    fig, axes = style.subplots_mm(W, H, 1, len(lineages))
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
    base.fit(fig, wspace=0.34)
    if save:
        base.save(fig, "S4_CD8_TCells", "S4_E", "S4_E_adaptive_composition")
    return fig


def main():
    fr = frames()
    return base.run({
        "S4_E": (lambda: draw_E(fr),
                 lambda: A._panel_composition(fr["fractions"], fr["tests"])),
    })


if __name__ == "__main__":
    sys.exit(main())
