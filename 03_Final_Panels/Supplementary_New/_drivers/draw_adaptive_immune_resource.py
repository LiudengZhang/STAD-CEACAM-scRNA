#!/usr/bin/env python3
"""Supplementary Figure S5 E: adaptive immune minor-state composition (mean
fraction per group, stacked bars, one box per lineage).

Reads adaptive_state_fractions.csv and adaptive_state_tests.csv written by
04_Revision_Analyses/09_R2.2_Adaptive_Immune/scripts/adaptive_immune_resource.py;
writes S5_CD8_TCells/S5_E/.
Run: python draw_adaptive_immune_resource.py [--check]
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

# Panel size in mm.
W, H = 171.0, 90.0

def frames():
    """The two tables the analysis wrote."""
    return dict(
        fractions=pd.read_csv(base.require(
            OUT / "adaptive_state_fractions.csv", "S3E state fractions")),
        tests=pd.read_csv(base.require(
            OUT / "adaptive_state_tests.csv", "S3E state tests")),
    )

def draw_E(fr, save=True):
    """S5 E: minor-state composition, one stacked-bar plot per lineage."""
    base.apply_style()
    lineages = list(A.LINEAGES)
    fig, axes = style.subplots_mm(W, H, 1, len(lineages))
    for ax, lineage in zip(axes, lineages):
        sub = fr["fractions"][fr["fractions"]["lineage"] == lineage]
        states = sorted(sub["minor_cell_state"].unique())
        # A state absent from a sample counts as 0, as in the tests.
        per = sub.pivot_table(index=["sample", "group4"], columns="minor_cell_state",
                              values="fraction", fill_value=0.0)
        means = per.groupby(level="group4").mean().reindex(A.GROUPS).fillna(0)
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
        base.finish(ax)
        ax.legend(loc="upper left", bbox_to_anchor=(0, -0.22), ncol=1,
                  handlelength=1.0, handletextpad=0.4, labelspacing=0.25)
    axes[0].set_ylabel("Fraction of lineage")
    base.fit(fig, wspace=0.34)
    if save:
        base.save(fig, "S5_CD8_TCells", "S5_E", "S5_E_adaptive_composition")
    return fig

def main():
    fr = frames()
    return base.run({
        "S5_E": (lambda: draw_E(fr),
                 lambda: A._panel_composition(fr["fractions"], fr["tests"])),
    })

if __name__ == "__main__":
    sys.exit(main())
