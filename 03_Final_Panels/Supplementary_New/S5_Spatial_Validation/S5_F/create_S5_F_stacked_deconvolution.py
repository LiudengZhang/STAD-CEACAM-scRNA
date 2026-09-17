#!/usr/bin/env python3
"""
S4 panel F - GraphST deconvolution: mean cell-type proportions per GSE251950
section, drawn at the size it prints at.

  printed panel  Supplementary Figure S5 F   (PROVENANCE.csv)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/S4_Spatial_Validation/
                 S5_F/create_S4_F_stacked_deconvolution.py

Reads SPATIAL_SPOT_DATA, the per-spot table the predecessor read, and takes
the same per-section means. The key stands at the right in two columns
(cnsfig.legend.group_key); the predecessor's axes legend had the same
fourteen names. GC6-PM is printed in the same face as the other sections
(the page it is on carries no footnote to point at).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import SPATIAL_SPOT_DATA                       # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import group_key                              # noqa: E402
from cnsfig.boxes import box_xlim                         # noqa: E402
import _driver_base as base                               # noqa: E402
from _spatial_gallery import SAMPLE_INFO                  # noqa: E402

FIG, PANEL = "S5_Spatial_Validation", "S5_F"
CELL_TYPE_COLORS = {
    'Epi CEACAM-high': '#d95f02', 'Epi CEACAM-low': '#e5c494',
    'Monocytes/Macrophages': '#fc8d62', 'CD4+ T cells': '#66c2a5',
    'CD8+ T cells': '#1b9e77', 'NK cells': '#7570b3', 'B cells': '#ffd92f',
    'Plasma cells': '#8da0cb', 'DC cells': '#80b1d3', 'Endothelial cells': '#a6d854',
    'Fibroblast': '#e78ac3', 'Pericyte': '#bebada', 'Neutrophils': '#b3b3b3',
    'Mast cells': '#fb8072',
}
CELL_TYPE_ORDER = [
    'Epi CEACAM-high', 'Epi CEACAM-low', 'CD4+ T cells', 'CD8+ T cells', 'NK cells',
    'B cells', 'Plasma cells', 'Monocytes/Macrophages', 'DC cells', 'Neutrophils',
    'Mast cells', 'Endothelial cells', 'Fibroblast', 'Pericyte',
]

# The printed box, millimetres: the right half of the last row since
# 2026-09-16 (the author's fifth reading: "see whether F can be squeezed up
# so S4 becomes 3 x 2"). Until then it had the page width to itself with a
# two-column key; at 80 mm the twelve names stand in ONE column at the right
# (two columns of "Monocytes/Macrophages" and "Endothelial cells" are 47 mm
# at 6 pt, more than the key can have), the ten bars keep a 4.6 mm pitch,
# and the 75 mm title breaks into two lines (labels.py RENAMES_S1_S6).
W, H = 80.0, 40.0
KEY_W_MM = 27.0
KEY_TOP_MM = 2.5              # the key starts at the canvas top, right of the title
MARGIN = dict(left=8.0, right=KEY_W_MM + 1.5, top=7.0, bottom=9.0)
TITLE = 'GraphST Deconvolution:\nCell Type Proportions per Sample'


def draw(df):
    base.apply_style()
    cols = [c for c in CELL_TYPE_ORDER if c in df.columns]
    means = df.groupby('sample')[cols].mean().loc[[s for s, _ in SAMPLE_INFO]]
    labels = [l for _, l in SAMPLE_INFO]
    fig, ax = style.subplots_mm(W, H)
    style.margins_mm(fig, **MARGIN)
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    for ct in cols:
        vals = means[ct].to_numpy()
        ax.bar(x, vals, bottom=bottom, color=CELL_TYPE_COLORS.get(ct, '#999999'),
               label=ct, width=0.75, edgecolor='white', linewidth=style.EDGE_PT)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('Mean proportion')
    ax.set_ylim(0, 1)
    # At a 4.6 mm bar pitch matplotlib's default x margin leaves 0.6 mm
    # between the first bar and the y spine; the page check wants 1.0.
    ax.set_xlim(*box_xlim(x, 0.75, clear=0.5))
    ax.set_title(TITLE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    box = ax.get_position()
    group_key(fig, [(ct, CELL_TYPE_COLORS.get(ct, '#999999')) for ct in cols],
              x_mm=W - KEY_W_MM, y_mm=KEY_TOP_MM, ncol=1, linespacing=0.15)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"ink outside the {W} x {H} mm canvas: {over}")
    return fig


def main():
    base.save(draw(pd.read_csv(SPATIAL_SPOT_DATA)), FIG, PANEL, "panel_S5_F")
    return 0


if __name__ == "__main__":
    sys.exit(main())
