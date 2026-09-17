#!/usr/bin/env python3
"""
One 2 x 5 gallery of the ten GSE251950 spatial sections coloured by a per-spot
value, drawn at print size.

Supplementary panels S4 A-E were five copies of one submission-tree script differing
in the column, the colour map and the colour-bar label; this is that script's
drawing half once, with the sizes in millimetres. The spots are read from
SPATIAL_SPOT_DATA as before (a table already; nothing is computed here beyond
the 2nd/98th percentile colour limits the predecessors also took).

The section labels are the printed page's: "GC6-PM" without the asterisk and
without the footnote the predecessors set - the submitted S4 was hand-edited
to drop both (00_GROUND_TRUTH/README.md, amendment 2026-09-10), and the page
is the ground truth. Both declared in labels.py RENAMES/REMOVALS_S1_S6.
"""

import matplotlib.colors as mcolors
import pandas as pd

import panel_style_cns as style

SAMPLE_INFO = [
    ('sample_01', 'GC1'), ('sample_02', 'GC2'), ('sample_03', 'GC3'),
    ('sample_04', 'GC4'), ('sample_05', 'GC5'), ('sample_06', 'GC6'),
    ('sample_07', 'GC6-PM'), ('sample_08', 'GC7'), ('sample_09', 'GC8'),
    ('sample_10', 'GC9'),
]
SPOT_MM = 0.25            # a spot's diameter on the page


def draw(df, column, *, cmap, cbar_label, w_mm, h_mm, map_mm=12.5, left_mm=1.0,
         top_mm=3.8, gap_mm=1.0, row_gap_mm=3.6, cbar_w_mm=1.6, cbar_gap_mm=2.0,
         panel=""):
    """Returns fig. Axes: the ten maps in reading order, then the colour bar."""
    import matplotlib.pyplot as plt
    fig = style.figure_mm(w_mm, h_mm)
    vmin = df[column].quantile(0.02)
    vmax = df[column].quantile(0.98)
    for idx, (sample, label) in enumerate(SAMPLE_INFO):
        row, col = divmod(idx, 5)
        x = left_mm + col * (map_mm + gap_mm)
        y_top = top_mm + row * (map_mm + row_gap_mm)
        ax = fig.add_axes([x / w_mm, 1 - (y_top + map_mm) / h_mm, map_mm / w_mm, map_mm / h_mm])
        sub = df[df['sample'] == sample]
        if len(sub) == 0:
            ax.text(0.5, 0.5, f'{label}\nNo data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=style.tick_pt())
            ax.set_xticks([]); ax.set_yticks([])
            continue
        ax.scatter(sub['x'], sub['y'], c=sub[column], cmap=cmap,
                   s=(SPOT_MM * style.PT_PER_MM) ** 2, alpha=0.8, vmin=vmin, vmax=vmax,
                   rasterized=True, edgecolors='none', linewidths=0)
        ax.set_title(label, pad=1.5, fontsize=style.tick_pt())
        ax.set_aspect('equal')
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():           # a full frame round each section
            sp.set_visible(True)
            sp.set_linewidth(style.EDGE_PT)
    x_cb = left_mm + 5 * map_mm + 4 * gap_mm + cbar_gap_mm
    cax = fig.add_axes([x_cb / w_mm, 1 - (top_mm + map_mm + row_gap_mm + 0.8 * map_mm) / h_mm,
                        cbar_w_mm / w_mm, (map_mm + row_gap_mm + 0.6 * map_mm) / h_mm])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(cbar_label, fontsize=style.tick_pt())
    cb.ax.tick_params(labelsize=style.tick_pt(), width=style.RULE_PT, length=1.5)
    cb.outline.set_linewidth(style.EDGE_PT)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"{panel}: ink outside the {w_mm} x {h_mm} mm canvas: {over}")
    return fig
