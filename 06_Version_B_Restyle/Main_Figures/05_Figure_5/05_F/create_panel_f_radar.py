#!/usr/bin/env python3
"""
Figure 5 panel H, RESTYLED (Version B) - NF-kB NES radar across thirteen cell
types, pre versus post treatment.

Version A is
`03_Final_Panels/05_Figure_5/05_F/create_panel_f_radar.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

The source table, the `method == "ttest"` selection, the cell-type ORDER, the
LABELS and the radial limits are Version A's, unchanged. Version A's own header
records why the source moved this round (the submitted panel read a table built
on the doubly normalised .X); that decision is Version A's and is not revisited
here. A restyle changes how a panel is drawn and nothing else.

  printed panel  Figure 5 H        (PROVENANCE.csv - the directory is "05_F";
                 do NOT read the directory as the letter)
  printed rect   44.6 x 33.0 mm    (panel_rects.csv)
  Version B box  95.0 x 82.0 mm

    Thirteen cell-type names are set around the circle, so the label ring costs
    roughly twice the width of the longest name ("Endothelial", ~10 mm at 7 pt)
    on top of the plotting circle, and the legend sits outside the axes at the
    lower right. 44.6 mm holds none of that at 7 pt.

MARK
    Version A drew 20 x 20 cm at SCALE = 4 and set its smallest body type - the
    radial tick labels and the legend - at `4 * SCALE`, so SMALL_PT = 4 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 16 = 0.4375

    The two series line widths, the marker sizes, the dashed zero circle, the
    angular tick pad and the title pad are scaled by it. Grid and spine widths
    are not: those are style, and cnsplots sets them.

Judgement calls, stated plainly:
  - `handlelength`, `borderpad` and `labelspacing` are left alone: matplotlib
    expresses all three in units of the legend font size, so they already track
    the type and multiplying them by MARK would shrink them twice.
  - The legend keeps its explicit `edgecolor`/`framealpha`; cnsplots turns
    legend frames off by default, and Version A deliberately draws one here to
    lift the legend off the grid. Restoring `frameon=True` is a visual decision
    that keeps Version A's, and no plotted value depends on it.

NOTE FOR THE GATE
    This script writes `nfkb_rankings_used.csv` beside itself, as Version A
    does. `compare_panel_content.py` blocks disk writes, so this panel must be
    gated with `--divert-writes`, which swallows and lists that write.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *  # noqa: E402,F403
import panel_style_cns as style  # noqa: E402

BASE_DIR = Path(__file__).parent
SRC = (Path(__file__).resolve().parents[4] / "02_New_Analyses"
       / "07_R1.8_NFkB_Specificity" / "outputs" / "nfkb_per_celltype.csv")
METHOD = "ttest"

SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 4.0                  # Version A's smallest body type, before * SCALE

PRINTED_MM = (44.6, 33.0)       # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 95.0, 82.0
MARGIN = dict(left=15.0, right=24.0, top=9.0, bottom=9.0)

COLOR_PRE = '#2166AC'
COLOR_POST = '#B2182B'

# The names the submitted panel used, so the axis labels do not move.
LABELS = {"B_cells": "B cells", "DC_cells": "DC",
          "Endothelial_cells": "Endothelial", "Epithelial": "Epithelial",
          "Fibroblast": "Fibroblast", "Mast_cells": "Mast", "MoMac": "MoMac",
          "Neutrophils": "Neutrophils", "NK_cells": "NK",
          "Pericyte": "Pericyte", "Plasma_cells": "Plasma",
          "TCD4_cells": "CD4+ T", "TCD8_cells": "CD8+ T"}
ORDER = ["B_cells", "DC_cells", "Endothelial_cells", "Epithelial", "Fibroblast",
         "Mast_cells", "MoMac", "Neutrophils", "NK_cells", "Pericyte",
         "Plasma_cells", "TCD4_cells", "TCD8_cells"]


def load():
    if not SRC.exists():
        raise SystemExit(f"missing {SRC} - run 07_R1.8_NFkB_Specificity first")
    df = pd.read_csv(SRC)
    df = df[df["method"] == METHOD]
    wide = df.pivot(index="cell_type", columns="phase",
                    values=["nes", "fdr_q", "rank", "n_sets"])
    missing = [c for c in ORDER if c not in wide.index]
    if missing:
        raise SystemExit(f"cell types absent from {SRC.name}: {missing}")
    out = pd.DataFrame({
        "cell_type": ORDER,
        "label": [LABELS[c] for c in ORDER],
        "pre_nes": [wide.loc[c, ("nes", "pre")] for c in ORDER],
        "post_nes": [wide.loc[c, ("nes", "post")] for c in ORDER],
        "pre_fdr": [wide.loc[c, ("fdr_q", "pre")] for c in ORDER],
        "post_fdr": [wide.loc[c, ("fdr_q", "post")] for c in ORDER],
        "pre_rank": [wide.loc[c, ("rank", "pre")] for c in ORDER],
        "post_rank": [wide.loc[c, ("rank", "post")] for c in ORDER],
        "pre_n_sets": [wide.loc[c, ("n_sets", "pre")] for c in ORDER],
        "post_n_sets": [wide.loc[c, ("n_sets", "post")] for c in ORDER],
    })
    if out[["pre_nes", "post_nes"]].isna().any().any():
        raise SystemExit("NES missing for at least one cell type or timepoint")
    return out


def main():
    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

    df = load()

    labels = df['label'].tolist()
    n_vars = len(labels)
    angles = np.linspace(0, 2 * np.pi, n_vars, endpoint=False).tolist()
    angles_closed = angles + [angles[0]]

    pre_nes = df['pre_nes'].tolist()
    post_nes = df['post_nes'].tolist()
    pre_nes_closed = pre_nes + [pre_nes[0]]
    post_nes_closed = post_nes + [post_nes[0]]

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM,
                                subplot_kw=dict(polar=True))

    ax.fill(angles_closed, pre_nes_closed, color=COLOR_PRE, alpha=0.10)
    ax.fill(angles_closed, post_nes_closed, color=COLOR_POST, alpha=0.10)
    ax.plot(angles_closed, pre_nes_closed, 'o--', linewidth=2 * MARK,
            color=COLOR_PRE, label='Pre-treatment', markersize=6 * MARK)

    ax.plot(angles_closed, post_nes_closed, 'o-', linewidth=2 * MARK,
            color=COLOR_POST, label='Post-treatment', markersize=6 * MARK)

    ax.set_xticks(angles)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='x', pad=2 * SCALE * MARK)
    ax.set_ylim(-2, 2.4)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_yticklabels(['-2', '-1', '0', '1', '2'], color='gray')
    ax.set_rlabel_position(90)
    ax.legend(loc='upper left', bbox_to_anchor=(0.99, 0.04),
              handlelength=1.6, borderpad=0.35, labelspacing=0.35,
              edgecolor='0.6', framealpha=1.0, frameon=True)
    ax.grid(True, linestyle='-', alpha=0.3)
    ax.set_title('NES', pad=15 * MARK)

    # Dashed zero circle
    theta_circle = np.linspace(0, 2 * np.pi, 100)
    ax.plot(theta_circle, [0] * 100, 'k--', linewidth=1 * MARK, alpha=0.7)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, BASE_DIR / 'nfkb_radar_celltype_enrichment')
    print(f"Saved: nfkb_radar_celltype_enrichment.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")

    df.to_csv(BASE_DIR / 'nfkb_rankings_used.csv', index=False)
    print(df[['label', 'pre_nes', 'post_nes']].to_string(index=False))


if __name__ == "__main__":
    main()
