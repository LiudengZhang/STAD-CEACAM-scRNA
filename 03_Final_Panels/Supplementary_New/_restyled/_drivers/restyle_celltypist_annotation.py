#!/usr/bin/env python3
"""
Restyle S11D - independent CellTypist annotation (Reviewer 1 point R1.7).

Analysis: 02_New_Analyses/06_R1.7_MoMac_Lineage_Markers/scripts/celltypist_annotation.py

NOTHING IS RECOMPUTED HERE. `celltypist` is never imported, no model is
downloaded, and MoMac.h5ad is never opened. The driver reads the per-cell label
table `main()` already wrote and rebuilds `tab` with the analysis's own
one-liner (celltypist_annotation.py:75), so the crosstab the panel draws is the
same object the report was written from.

    conda run -n Liudeng_Python_310 python restyle_celltypist_annotation.py
    conda run -n Liudeng_Python_310 python restyle_celltypist_annotation.py --check
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _driver_base as base
import panel_style_cns as style

A = base.analysis("06_R1.7_MoMac_Lineage_Markers/scripts/celltypist_annotation.py")
OUT = base.outputs_of(A)
FIG = "S11_Affirmative_Analyses"

# Printed box, millimetres. Version A's box was 171 x 44. At 7 pt the seven
# cell-state names on the left need 42 mm and the CellTypist legend on the
# right needs 44 mm, which leaves 85 mm of bar track; the height grows so the
# legend (one line per label at 7 pt) still clears the canvas.
W, H = 171.0, 52.0


def frames():
    """The label table the analysis already wrote, crosstabbed as it does."""
    df = pd.read_csv(base.require(
        OUT / "celltypist_labels.csv", "S11D CellTypist labels"))
    for col in ("minor_cell_state", "celltypist"):
        if col not in df.columns:
            raise RuntimeError(f"celltypist_labels.csv has no `{col}` column: "
                               f"columns are {list(df.columns)}")
    # celltypist_annotation.py:75, verbatim.
    tab = pd.crosstab(df["minor_cell_state"], df["celltypist"],
                      normalize="index") * 100
    return dict(tab=tab)


# --------------------------------------------------------------------------
# S11D - percentage of each cluster assigned to each CellTypist label.
# Layout unchanged: one stacked horizontal bar per cluster.
# --------------------------------------------------------------------------
def draw_D(fr, save=True):
    style.apply()
    tab = fr["tab"]
    keep = tab.loc[:, tab.max(axis=0) >= 5]
    keep = keep.reindex(sorted(keep.index))
    fig, ax = style.subplots_mm(W, H)
    bottom = np.zeros(len(keep))
    cmap = plt.get_cmap("tab20")
    for i, c in enumerate(keep.columns):
        ax.barh(np.arange(len(keep)), keep[c], left=bottom, height=0.7,
                color=cmap(i % 20), edgecolor="white", linewidth=0.4,
                label=str(c)[:28])
        bottom += keep[c].values
    ax.set_yticks(np.arange(len(keep)))
    ax.set_yticklabels([A.DISPLAY.get(s, s).replace("_", " ")
                        for s in keep.index])
    ax.set_xlabel("% of cells assigned by CellTypist")
    ax.set_xlim(0, 100)
    ax.tick_params(axis="both", width=0.6, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0),
              handlelength=1.0, handletextpad=0.4, labelspacing=0.35)
    style.margins_mm(fig, left=42, right=44, top=3, bottom=12)
    if save:
        base.save(fig, FIG, "S11_D", "S11_D_celltypist_annotation")
    return fig


def main():
    fr = frames()
    return base.run({
        "S11_D": (lambda: draw_D(fr), lambda: A._panel(fr["tab"])),
    })


if __name__ == "__main__":
    sys.exit(main())
