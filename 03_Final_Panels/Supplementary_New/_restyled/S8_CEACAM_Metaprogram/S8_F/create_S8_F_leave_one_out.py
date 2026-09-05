"""
Supplementary Figure S8F, RESTYLED (Version B) - leave-one-patient-out
stability of the four-versus-four pre-treatment comparison.

Version A of this panel is
`Supplementary_New/S8_CEACAM_Metaprogram/S8_F/create_S8_F_leave_one_out.py`
and is frozen. This file is a copy of it with three changes and no others:

  1. the type comes from `00_Config/panel_style_cns.py` (cnsplots), not from a
     local rcParams block and not from `5 * SCALE`;
  2. the canvas is the millimetre box the panel prints in, not four times it;
  3. everything specified in points that is *not* type - marker areas, marker
     edge widths, the zero rule - is rescaled by the factor below so that its
     size relative to the type is unchanged.

Not one value read, filtered, computed or plotted differs from Version A. The
drawing code is the same code.

THE MARK RESCALE FACTOR
-----------------------
Version A drew on a 40 x 12 cm canvas (SCALE = 4) and the assembler then fitted
that 1133.9 x 340.2 pt SVG into a 100 x 26 mm box, a fit of 0.2167. So a mark
set to `s = 22` pt^2 printed at sqrt(22) x 0.2167 = 1.016 pt across, and type
set to `5 * SCALE` = 20 pt printed at 4.334 pt.

Version B draws 1:1, so a size set here is the size printed. Type is now 7 pt,
which is 7 / 4.334 = 1.615x its Version A printed size. Marks are multiplied by
the same 1.615 so the panel keeps its proportions:

    printed_B_length = printed_A_length x 1.615
    L_B = L_A x 0.2167 x 1.615 = L_A x 0.35          (lengths: MARK)
    s_B = s_A x 0.35^2          = s_A x 0.1224       (areas:   AREA)

Input : 04_Revision_Analyses/10_R1.3_CrossCohort_Convergence/outputs/loo_stability.csv
Output: _restyled/S8_CEACAM_Metaprogram/S8_F/S8_F_leave_one_out.{svg,pdf,png}
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[5] / "00_Config"))
from paths import NEW_ANALYSES  # noqa: E402
import panel_style_cns as style  # noqa: E402

OUT_DIR = Path(__file__).parent
LOO = (NEW_ANALYSES / "10_R1.3_CrossCohort_Convergence" / "outputs"
       / "loo_stability.csv")

# The printed box. Version A's box was 100 x 26 mm and its panel's natural
# aspect was 10:3, so the assembler letterboxed 4 mm of dead space at the
# bottom. Drawing 100 x 32 mm removes the letterbox; the extra 2 mm over the
# 10:3 aspect is the room the legend needs now that it is set at 7 pt instead
# of printing at 4.3 pt. The plotted content is unchanged.
PANEL_W_MM = 100.0
PANEL_H_MM = 32.0

# Margins in millimetres of paper. Version A's fractions (left=0.12,
# bottom=0.34) were fractions of a 400 x 120 mm canvas; carried onto this one
# unchanged they clip the tick labels and the legend. See style.margins_mm.
MARGIN = dict(left=14.0, right=6.0, top=1.5, bottom=13.0)

# See THE MARK RESCALE FACTOR above.
A_FIT = 0.2167                     # assembler fit of Version A
A_TYPE_PT = 5.0 * 4 * A_FIT        # 4.334 pt, what Version A's small type printed at
MARK = (style.tick_pt() / A_TYPE_PT) * A_FIT   # length multiplier, 0.35
AREA = MARK ** 2                               # area multiplier, 0.1224

GENE_COLOR = {"CEACAM6": "#2166AC", "CEACAM5": "#4393C3"}


def main():
    family = style.apply()
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt")

    df = pd.read_csv(LOO)
    genes = ["CEACAM6", "CEACAM5"]

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    ax.axvline(0, color="#999999", linewidth=0.8 * MARK, linestyle="--",
               zorder=1)

    for row, gene in enumerate(genes):
        y = len(genes) - 1 - row
        d = df[df.Gene == gene]
        drops = d[d.Dropped != "none (as published)"]
        published = d[d.Dropped == "none (as published)"].iloc[0]
        c = GENE_COLOR[gene]
        ax.scatter(drops["Hedges g"], np.full(len(drops), y), s=22 * AREA, c=c,
                   alpha=0.75, edgecolors="white", linewidths=0.3 * MARK,
                   zorder=3)
        ax.scatter(published["Hedges g"], y, s=70 * AREA, marker="D", c="white",
                   edgecolors=c, linewidths=1.0 * MARK, zorder=4)
        ax.text(3.9, y + 0.30,
                f"{int((drops['Hedges g'] > 0).sum())} of {len(drops)} refits keep "
                f"the direction; $g$ {drops['Hedges g'].min():.2f}"
                f"–{drops['Hedges g'].max():.2f}",
                fontsize=style.tick_pt(), ha="right", va="center",
                color="#555555")

    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes[::-1])
    ax.set_ylim(-0.6, len(genes) - 0.4)
    ax.set_xlim(-0.4, 4.0)
    ax.set_xlabel("Effect size after dropping one patient (Hedges' $g$)")
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    # cnsplots draws legend keys at legend.markerscale (0.5). Version A had no
    # such reduction, so carrying its handle sizes over would halve the keys
    # and the open diamond would close up into a dot. Size each key so that,
    # after the scale, it prints at the same diameter as the mark it stands
    # for: sqrt(s) for a scatter, divided by the markerscale.
    key = 1.0 / plt.rcParams["legend.markerscale"]
    handles = [
        plt.Line2D([], [], marker="D", color="white", markeredgecolor="#444444",
                   markersize=(70 * AREA) ** 0.5 * key, linestyle="none",
                   label="all eight patients"),
        plt.Line2D([], [], marker="o", color="#444444",
                   markersize=(22 * AREA) ** 0.5 * key, linestyle="none",
                   label="one patient dropped"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.34),
              ncol=2, handletextpad=0.4, columnspacing=1.6)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUT_DIR / "S8_F_leave_one_out")
    print(f"Saved {OUT_DIR / 'S8_F_leave_one_out'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
