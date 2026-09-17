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

#: Paper left round the ring for the thirteen full spoke names (2026-09-14).
RING_MARGIN_MM = 12.0           # unused since the evening of 2026-09-14
#: The ring's radius and centre on the canvas, in mm. r_max (2.4 NES) is the
#: ring's edge; the names sit SPOKE_PAD outside it.
#: 2026-09-15 (the author's third reading): "NES" is the title, centred over
#: the ring as the page prints it, so the ring moves down to leave it a line;
#: the legend goes under the names at the bottom right instead of over
#: "Plasma"; the radius gives up 1 mm for both.
RING_R_MM = 9.8
RING_CX_MM = 23.0
RING_CY_MM = 17.3
#: Radial distance, in NES units, from the outer ring to the start of a name.
SPOKE_PAD = 0.4
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

# FIVE SPOKES KEEP A SHORT FORM, AND THIS IS WHERE THE ARITHMETIC STOPS
#   The 2026-09-11 round put the published wording back across Figures 2, 3 and
#   5 by giving the panels the blank bottom of the page. That works wherever
#   the crowding is vertical. A radar's labels ring a circle, so its crowding
#   is radial, and the extra height buys almost nothing: rebuilt with all
#   thirteen names in full, the panel gate convicted three pairs - Epithelial
#   against Fibroblast at 1.93 mm2 the worst - and clipped the C off CD4+ T.
#
#   The panel is 44.6 mm wide and that number did not change; widening it would
#   move a column of the page. So these five keep the short form, and the five
#   full names are carried in LABELS above and in the table this script writes
#   beside itself. The series names in the legend are NOT abbreviated - see
#   SERIES_LABEL - because that legend is a box, not a ring.
SPOKE_LABEL = {"Endothelial": "Endo", "Epithelial": "Epi",
               "Fibroblast": "Fibro", "Neutrophils": "Neut",
               "Pericyte": "Peri"}

# The two series keep a short form, for the same reason the five spokes above
# do. Both placements were built and measured on 2026-09-11: at full length in
# the published bottom-right corner the box reaches into the label ring and the
# gate convicted it against the CD4+ T spoke, clipping two glyphs; laid across
# the foot of the panel it touches the Plasma spoke at 0.00 pt. The radar fills
# its 44.6 x 38.0 mm and there is no corner left for a 32 mm key.
#
# Pre and Post are not an invention here. Panels D and E of this same figure
# print their four groups as Pre NR, Pre R, Post NR and Post R, so the reader
# meets the short form twice before reaching this panel.
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
    # Markers 0.41 mm across, counted off the published page (2026-09-14).
    ax.plot(angles_closed, pre_nes_closed, 'o--', linewidth=style.RULE_PT,
            color=COLOR_PRE, label='Pre-treatment',
            markersize=0.41 * style.PT_PER_MM)

    ax.plot(angles_closed, post_nes_closed, 'o-', linewidth=style.RULE_PT,
            color=COLOR_POST, label='Post-treatment',
            markersize=0.41 * style.PT_PER_MM)

    # THE THIRTEEN NAMES IN FULL, ROUND A SMALLER RING  (2026-09-14)
    #   The five short spokes and the two-word series key were the arithmetic
    #   of a ring that filled the canvas. The author asked for the page's
    #   names; the ring is drawn inside RING_MARGIN_MM of paper on each side
    #   and the names have the margin. Same thirteen spokes, same values.
    ax.set_xticks(angles)
    ax.set_xticklabels([])
    # Each name reads along its own spoke, outward, so neighbours near the
    # top and bottom of the ring - where horizontal names collided - are
    # side by side instead of on top of each other. Drawn as text: a polar
    # axis re-lays its tick labels at draw time and drops a rotation set on
    # them.
    # HORIZONTAL NAMES ROUND A 12 mm RING  (2026-09-14, evening)
    #   The radial names of the morning were rejected; the names are set
    #   horizontally, as the page sets them, each anchored on its spoke just
    #   outside the ring and aligned away from the centre (ha by the cosine,
    #   va by the sine). The ring is placed at millimetres (RING_R_MM,
    #   RING_CX_MM, RING_CY_MM) so the widest names - Neutrophils to the left,
    #   B cells to the right - end inside the canvas, and the legend sits
    #   under the ring at the right.
    r_label = 2.4 + SPOKE_PAD
    for name, ang in zip(labels, angles):
        c, s_ = np.cos(ang), np.sin(ang)
        ha = 'center' if abs(c) < 0.2 else ('left' if c > 0 else 'right')
        va = 'center' if abs(s_) < 0.2 else ('bottom' if s_ > 0 else 'top')
        ax.text(ang, r_label, name, ha=ha, va=va, fontsize=style.tick_pt())
    ax.tick_params(axis='x', pad=2 * SCALE * MARK)
    ax.set_ylim(-2, 2.4)
    # Three ring numbers, not five: at this ring size five stand 1.5 mm
    # apart and overprint.
    ax.set_yticks([-2, 0, 2])
    ax.set_yticklabels(['-2', '0', '2'], color='gray')
    # Between the B cells and DC spokes at the right, where the published
    # panel sets its ring numbers (2026-09-15; they had been at 180, along
    # the left horizontal, which the author read as "a strange place").
    ax.set_rlabel_position(10.0)   # 22.5 put "2" 0.17 pt from the DC name
    # Under the ring's names at the bottom right, clear of "Plasma"
    # (2026-09-15). See SERIES_LABEL for why the two names are the short ones.
    ax.legend(loc='lower right', bbox_to_anchor=(1.0 - 0.4 / PANEL_W_MM, 0.4 / PANEL_H_MM),
              bbox_transform=fig.transFigure,
              handlelength=1.4, borderpad=0.2, labelspacing=0.3,
              edgecolor='0.6', framealpha=1.0, frameon=True)
    ax.grid(True, linestyle='-', alpha=0.3)
    # The title, centred over the ring as the page prints it (2026-09-15).
    ax.set_title('')
    fig.text(RING_CX_MM / PANEL_W_MM, 1.0 - 0.4 / PANEL_H_MM, 'NES',
             ha='center', va='top', fontsize=style.body_pt())

    # Dashed zero circle
    theta_circle = np.linspace(0, 2 * np.pi, 100)
    ax.plot(theta_circle, [0] * 100, 'k--', linewidth=style.RULE_PT, alpha=0.7)

    ax.set_position([(RING_CX_MM - RING_R_MM) / PANEL_W_MM,
                     1.0 - (RING_CY_MM + RING_R_MM) / PANEL_H_MM,
                     2 * RING_R_MM / PANEL_W_MM, 2 * RING_R_MM / PANEL_H_MM])
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
