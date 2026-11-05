#!/usr/bin/env python3
"""Supplementary Figure S6 F: GraphST deconvolution, mean proportions of
fourteen cell types per GSE251950 section (stacked bars).

Reads SPATIAL_SPOT_DATA (per-spot table); writes the panel through
_driver_base.save. Run: python create_S6_F_stacked_deconvolution.py"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import SPATIAL_SPOT_DATA                       # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import group_key                              # noqa: E402
from cnsfig.boxes import box_xlim                         # noqa: E402
import _driver_base as base                               # noqa: E402
from _spatial_gallery import SAMPLE_INFO                  # noqa: E402

FIG, PANEL = "S6_Spatial_Validation", "S6_F"
CELL_TYPE_COLORS = {
    'Epi CEACAM-high': '#d95f02', 'Epi CEACAM-low': '#e5c494',
    'Monocytes/Macrophages': '#fc8d62', 'CD4+ T cells': '#66c2a5',
    'CD8+ T cells': '#1b9e77', 'NK cells': '#7570b3', 'B cells': '#ffd92f',
    'Plasma cells': '#8da0cb', 'DC cells': '#80b1d3', 'Endothelial cells': '#a6d854',
    'Fibroblast': '#e78ac3', 'Pericyte': '#bebada', 'Neutrophils': '#b3b3b3',
    'Mast cells': '#fb8072',
}
# Names shared with the major-cell-type palette take its colours.
from shared.figure_config import CELL_TYPE_COLORS as _MAJOR  # noqa: E402
CELL_TYPE_COLORS.update({k: v for k, v in _MAJOR.items() if k in CELL_TYPE_COLORS})
CELL_TYPE_ORDER = [
    'Epi CEACAM-high', 'Epi CEACAM-low', 'CD4+ T cells', 'CD8+ T cells', 'NK cells',
    'B cells', 'Plasma cells', 'Monocytes/Macrophages', 'DC cells', 'Neutrophils',
    'Mast cells', 'Endothelial cells', 'Fibroblast', 'Pericyte',
]

# Panel size in mm; the key is one column at the right.
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
    # Keep the first bar clear of the y spine.
    ax.set_xlim(*box_xlim(x, 0.75, clear=0.5))
    ax.set_title(TITLE)
    base.despine(ax)
    group_key(fig, [(ct.replace('CEACAM-', 'CEACAM5/6-'),
                     CELL_TYPE_COLORS.get(ct, '#999999')) for ct in cols],
              x_mm=W - KEY_W_MM, y_mm=KEY_TOP_MM, ncol=1, linespacing=0.15)
    return fig

if __name__ == "__main__":
    base.save(draw(pd.read_csv(SPATIAL_SPOT_DATA)), FIG, PANEL, "panel_S6_F")
