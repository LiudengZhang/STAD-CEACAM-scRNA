#!/usr/bin/env python3
"""
Figure 5 panel H - NF-kB NES radar across thirteen cell types, pre versus post
treatment.

The source table, the `method == "ttest"` selection, the cell-type ORDER and
the radial limits are unchanged. The panel reads the adopted sound-input table,
04_Revision_Analyses/13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/
nfkb_per_celltype_sound13.csv - the same table the Results, the response letter
and the two supplementary NF-kB panels read - so the three NF-kB panels and the
text draw from one table. Do not repoint it at
07_R1.8_NFkB_Specificity/outputs/nfkb_per_celltype.csv, the earlier run on the
doubly normalised .X: against that table one of the twenty-six spokes changes
sign.

  printed panel  Figure 5 H       (PROVENANCE.csv - the directory is "05_F";
                                   do NOT read the directory as the letter)

Thirteen cell-type names are set around the circle, so the label ring costs
roughly twice the width of the longest name on top of the plotting circle, and
the legend sits outside the axes at the lower right. The names are therefore
printed in a short form, declared in 00_Config/shared/labels.py; the spokes,
their order and the values on them do not move.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the radial tick labels and the legend - at 4 * SCALE,
    so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The two series line widths, the marker sizes, the dashed zero circle, the
    angular tick pad and the title pad are scaled by it. Grid and spine widths
    are not: those are style, and cnsplots sets them.

Judgement calls, stated plainly:
  - `handlelength`, `borderpad` and `labelspacing` are left alone: matplotlib
    expresses all three in units of the legend font size, so they already track
    the type and multiplying them by MARK would shrink them twice.
  - The legend keeps its explicit `edgecolor`/`framealpha`; cnsplots turns
    legend frames off by default, and the frame is drawn here deliberately to
    lift the legend off the grid. No plotted value depends on it.

NOTE FOR THE GATE
    This script writes `nfkb_rankings_used.csv` beside itself.
    `compare_panel_content.py` blocks disk writes, so this panel must be gated
    with `--divert-writes`, which swallows and lists that write.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *  # noqa: E402,F403
import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402

BASE_DIR = Path(__file__).parent
# NEW_ANALYSES comes from 00_Config/paths.py, which names this tree correctly
# in both layouts: 04_Revision_Analyses here, 04_Revision_Analyses in the release.
SRC = (NEW_ANALYSES
       / "13_R1.8_Neutrophil_Rebuilt_Recompute" / "outputs"
       / "nfkb_per_celltype_sound13.csv")
METHOD = "ttest"

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 4.0                  # the earlier smallest body type, before * SCALE

PANEL_LETTER = "H"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(5, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(5, PANEL_LETTER)

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

# Thirteen names set around a 33 mm circle collide at full length: the label
# ring costs twice the longest name on top of the plotting circle. These are
# the forms printed on the spokes. LABELS above is untouched, so the table this
# script writes beside itself still carries the full names.
SPOKE_LABEL = {"Endothelial": "Endo", "Epithelial": "Epi",
               "Fibroblast": "Fibro", "Neutrophils": "Neut",
               "Pericyte": "Peri"}

# The two series, named as the group labels of panels D and E name them. At
# full length the legend is more than half the width of this panel and prints
# over the lower-right spoke names.
SERIES_LABEL = {"Pre-treatment": "Pre", "Post-treatment": "Post"}


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
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
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
            color=COLOR_PRE, label=SERIES_LABEL['Pre-treatment'],
            markersize=6 * MARK)

    ax.plot(angles_closed, post_nes_closed, 'o-', linewidth=2 * MARK,
            color=COLOR_POST, label=SERIES_LABEL['Post-treatment'],
            markersize=6 * MARK)

    ax.set_xticks(angles)
    ax.set_xticklabels([SPOKE_LABEL.get(t, t) for t in labels])
    ax.tick_params(axis='x', pad=2 * SCALE * MARK)
    ax.set_ylim(-2, 2.4)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_yticklabels(['-2', '-1', '0', '1', '2'], color='gray')
    ax.set_rlabel_position(90)
    # Pinned to the bottom-right corner of the canvas, where the published
    # panel prints it and where the label ring has nothing: anchored to the
    # axes it lands on the lower-right spoke names, because at 1:1 the legend
    # is a much larger fraction of the panel than it was at four times the size.
    ax.legend(loc='lower right', bbox_to_anchor=(1.0, 0.0),
              bbox_transform=fig.transFigure,
              handlelength=1.6, borderpad=0.35, labelspacing=0.35,
              edgecolor='0.6', framealpha=1.0, frameon=True)
    ax.grid(True, linestyle='-', alpha=0.3)
    ax.set_title('NES', pad=15 * MARK)

    # Dashed zero circle
    theta_circle = np.linspace(0, 2 * np.pi, 100)
    ax.plot(theta_circle, [0] * 100, 'k--', linewidth=1 * MARK, alpha=0.7)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    style.save_panel(fig, BASE_DIR / 'nfkb_radar_celltype_enrichment')
    print(f"Saved: nfkb_radar_celltype_enrichment.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")

    df.to_csv(BASE_DIR / 'nfkb_rankings_used.csv', index=False)
    print(df[['label', 'pre_nes', 'post_nes']].to_string(index=False))


if __name__ == "__main__":
    main()
