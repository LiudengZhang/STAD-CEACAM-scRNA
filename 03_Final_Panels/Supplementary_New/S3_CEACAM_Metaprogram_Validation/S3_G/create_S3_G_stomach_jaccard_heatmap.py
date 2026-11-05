#!/usr/bin/env python3
"""Supplementary Figure S3 G: Jaccard similarity of all stomach NMF programs
(lower triangle), ordered by metaprogram size.

Reads the Jaccard matrix and metaprogram assignments from NMF_INTERMEDIATE;
writes the panel through _driver_base.save.
Run: python create_S3_G_stomach_jaccard_heatmap.py"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import NMF_INTERMEDIATE                        # noqa: E402
import panel_style_cns as style                           # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_G"
MP_COLORS = {'MP1': '#E41A1C', 'MP2': '#377EB8', 'MP3': '#4DAF4A',
             'MP4': '#984EA3', 'MP5': '#FF7F00', 'MP6': '#FFFF33', 'Other': '#999999'}
CMAP_COLORS = ['#FFFFCC', '#FFEDA0', '#FED976', '#FEB24C', '#FD8D3C',
               '#FC4E2A', '#E31A1C', '#BD0026', '#800026']

# Panel size in mm.
W, H = 61.0, 48.0
MAP_MM = 36.0
LEFT_MM, TOP_MM = 6.0, 8.5
CBAR_W_MM = 1.6
LABEL_GAP = 2.6 / (36.0 / 116)     # 2.6 mm between block names, in cells
LABEL_TOP = 1.5                    # the first name's centre row, in cells

def load():
    matrix = pd.read_csv(NMF_INTERMEDIATE / "panel_A2_stomach_all_jaccard_matrix.csv", index_col=0)
    assignments = pd.read_csv(NMF_INTERMEDIATE / "panel_A2_stomach_all_mp_assignments.csv")
    sizes = assignments.groupby('metaprogram_id').size().sort_values(ascending=False)
    ordered, info = [], []
    for mp in sizes.index:
        progs = [p for p in assignments[assignments['metaprogram_id'] == mp]['program_id']
                 if p in matrix.index]
        ordered.extend(progs)
        if progs:
            info.append({'mp_id': mp, 'size': len(progs)})
    unassigned = [p for p in matrix.index if p not in set(ordered)]
    if unassigned:
        ordered.extend(unassigned)
        info.append({'mp_id': 'Other', 'size': len(unassigned)})
    return matrix.loc[ordered, ordered].values, info

def draw(data, info):
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle
    base.apply_style()
    fig = style.figure_mm(W, H)
    ax = fig.add_axes([LEFT_MM / W, 1 - (TOP_MM + MAP_MM) / H, MAP_MM / W, MAP_MM / H])
    cmap = LinearSegmentedColormap.from_list('stomach', CMAP_COLORS)
    n = data.shape[0]
    mask = np.triu(np.ones((n, n), dtype=bool), k=0)
    im = ax.imshow(np.ma.masked_where(mask, data), cmap=cmap, aspect='equal',
                   vmin=0, vmax=0.5, rasterized=False)
    bar_px = 4.375                          # colour-bar width, in cells
    bounds = [0]
    for mp in info:
        bounds.append(bounds[-1] + mp['size'])
    for b in bounds[1:-1]:
        ax.plot([-0.5, b - 0.5], [b - 0.5, b - 0.5], color='black', linestyle='--',
                linewidth=style.RULE_PT, clip_on=True)
        ax.plot([b - 0.5, b - 0.5], [b - 0.5, n - 0.5], color='black', linestyle='--',
                linewidth=style.RULE_PT, clip_on=True)
    for i, mp in enumerate(info):
        start, end = bounds[i], bounds[i + 1]
        color = MP_COLORS.get(mp['mp_id'], '#999999')
        ax.add_patch(Rectangle((-bar_px - 0.5, start - 0.5), bar_px, end - start,
                               facecolor=color, edgecolor='black', linewidth=style.EDGE_PT))
        ax.add_patch(Rectangle((start - 0.5, n - 0.5), end - start, bar_px,
                               facecolor=color, edgecolor='black', linewidth=style.EDGE_PT))
    # Block names in the empty upper triangle; raise if one reaches the diagonal.
    labels = []
    for i, mp in enumerate(info):
        labels.append(ax.text(n - 1.0, LABEL_TOP + i * LABEL_GAP, mp['mp_id'],
                              ha='right', va='center', fontsize=style.tick_pt(),
                              color=MP_COLORS.get(mp['mp_id'], '#999999')))
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for t in labels:
        bb = t.get_window_extent(r)
        (x0, y_low), (x1, y_high) = ax.transData.inverted().transform(
            [(bb.x0, bb.y0), (bb.x1, bb.y1)])
        y_bottom = max(y_low, y_high)                # image rows grow downward
        if x0 <= y_bottom + 1.0:
            raise RuntimeError(f"block name {t.get_text()!r} reaches the "
                               f"diagonal: left edge {x0:.1f}, bottom row "
                               f"{y_bottom:.1f} (cells)")
        if x1 > n - 0.5 + 1e-6 or y_low < -0.5:
            raise RuntimeError(f"block name {t.get_text()!r} leaves the matrix")
    cax = fig.add_axes([(LEFT_MM + MAP_MM + 4.5) / W, 1 - (TOP_MM + MAP_MM * 0.8) / H,
                        CBAR_W_MM / W, MAP_MM * 0.6 / H])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Jaccard similarity', fontsize=style.tick_pt())
    cb.ax.tick_params(labelsize=style.tick_pt(), width=style.RULE_PT, length=1.5)
    cb.outline.set_linewidth(style.EDGE_PT)
    ax.set_xlim(-bar_px - 1, n + bar_px + 0.5)
    ax.set_ylim(n + bar_px + 0.5, -bar_px - 1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect('equal')
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title('Epithelial Meta-Programs\n(Stomach, all samples)', pad=2)
    return fig

if __name__ == "__main__":
    data, info = load()
    base.save(draw(data, info), FIG, PANEL, "panel_S3_G")
