#!/usr/bin/env python3
"""
Panel 2P: Stomach-only NMF Meta-Program Jaccard Similarity Heatmap.
Lower-left triangle of all stomach NMF programs (32 samples), ordered by MP1–MP5.
4× scaling method for Nature Cancer.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

BASE_DIR = Path(__file__).parent
INTERMEDIATE = NMF_INTERMEDIATE

# 4× scaling — canvas + fonts only
DPI = 300
SCALE = 4
CM_TO_INCH = 1 / 2.54

# MP color palette (consistent with S03.1 bowtie)
MP_COLORS = {
    'MP1': '#E41A1C', 'MP2': '#377EB8', 'MP3': '#4DAF4A',
    'MP4': '#984EA3', 'MP5': '#FF7F00', 'MP6': '#FFFF33',
    'Other': '#999999',
}

# YlOrRd colormap (same as stomach triangle in bowtie)
CMAP = LinearSegmentedColormap.from_list('stomach',
    ['#FFFFCC', '#FFEDA0', '#FED976', '#FEB24C', '#FD8D3C',
     '#FC4E2A', '#E31A1C', '#BD0026', '#800026'])


def load_data():
    matrix = pd.read_csv(INTERMEDIATE / "panel_A2_stomach_all_jaccard_matrix.csv", index_col=0)
    assignments = pd.read_csv(INTERMEDIATE / "panel_A2_stomach_all_mp_assignments.csv")
    return matrix, assignments


def order_by_mp(matrix, assignments):
    mp_sizes = assignments.groupby('metaprogram_id').size().sort_values(ascending=False)

    ordered_ids = []
    mp_info = []
    for mp_id in mp_sizes.index:
        mp_programs = assignments[assignments['metaprogram_id'] == mp_id]['program_id'].tolist()
        valid_programs = [p for p in mp_programs if p in matrix.index]
        ordered_ids.extend(valid_programs)
        if valid_programs:
            mp_info.append({'mp_id': mp_id, 'size': len(valid_programs)})

    all_ids = matrix.index.tolist()
    assigned_ids = set(ordered_ids)
    unassigned = [p for p in all_ids if p not in assigned_ids]

    if unassigned:
        ordered_ids.extend(unassigned)
        mp_info.append({'mp_id': 'Other', 'size': len(unassigned)})

    matrix_ordered = matrix.loc[ordered_ids, ordered_ids].values
    return matrix_ordered, mp_info


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("=" * 60)
    print("Panel 2P: Stomach Jaccard Heatmap (lower-left triangle)")
    print("=" * 60)

    matrix, assignments = load_data()
    data, mp_info = order_by_mp(matrix, assignments)
    n = data.shape[0]

    print(f"  {n} programs, MPs: {[(m['mp_id'], m['size']) for m in mp_info]}")

    # Mask upper-right triangle (including diagonal)
    mask = np.triu(np.ones((n, n), dtype=bool), k=0)
    data_masked = np.ma.masked_where(mask, data)

    # Figure sizing — cell-based
    cell_size = 0.08 * SCALE
    bar_width = 0.35 * SCALE
    bar_px = bar_width / cell_size
    fig_size = n * cell_size + 2.5 * SCALE

    fig, ax = plt.subplots(figsize=(fig_size * CM_TO_INCH, fig_size * CM_TO_INCH))

    # Lower-left triangle heatmap
    im = ax.imshow(data_masked, cmap=CMAP, aspect='equal', vmin=0, vmax=0.5,
                   rasterized=False)

    # MP boundary lines (dashed) — only in lower-left triangle
    mp_bounds = [0]
    for mp in mp_info:
        mp_bounds.append(mp_bounds[-1] + mp['size'])

    for b in mp_bounds[1:-1]:
        # Horizontal line: from x=0 to x=b (within lower-left)
        ax.plot([-0.5, b - 0.5], [b - 0.5, b - 0.5],
                color='black', linestyle='--', linewidth=0.6, clip_on=True)
        # Vertical line: from y=b to y=n (within lower-left)
        ax.plot([b - 0.5, b - 0.5], [b - 0.5, n - 0.5],
                color='black', linestyle='--', linewidth=0.6, clip_on=True)

    # MP color bars — LEFT and BOTTOM
    for i, mp in enumerate(mp_info):
        start = mp_bounds[i]
        end = mp_bounds[i + 1]
        color = MP_COLORS.get(mp['mp_id'], '#999999')

        # Left bar
        rect_l = Rectangle((-bar_px - 0.5, start - 0.5), bar_px, end - start,
                            facecolor=color, edgecolor='black', linewidth=0.5)
        ax.add_patch(rect_l)

        # Bottom bar
        rect_b = Rectangle((start - 0.5, n - 0.5), end - start, bar_px,
                            facecolor=color, edgecolor='black', linewidth=0.5)
        ax.add_patch(rect_b)

        # Labels
        mid = (start + end) / 2
        text_color = 'white' if mp['mp_id'] in ['MP1', 'MP2'] else 'black'
        ax.text(-bar_px / 2 - 0.5, mid, mp['mp_id'], ha='center', va='center',
                fontsize=5 * SCALE, rotation=90, color=text_color, fontweight='bold')
        ax.text(mid, n + bar_px / 2, mp['mp_id'], ha='center', va='center',
                fontsize=5 * SCALE, color=text_color, fontweight='bold')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04, shrink=0.6)
    cbar.set_label('Jaccard similarity', fontsize=6 * SCALE)
    cbar.ax.tick_params(labelsize=5 * SCALE, width=1.0, length=4)

    # Axis settings
    ax.set_xlim(-bar_px - 1, n + bar_px + 0.5)
    ax.set_ylim(n + bar_px + 0.5, -bar_px - 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')

    ax.set_title('Epithelial Meta-Programs\n(Stomach, all samples)',
                 fontsize=7 * SCALE, pad=12)

    plt.tight_layout()

    # Save
    out_png = BASE_DIR / 'stomach_jaccard_heatmap.png'
    fig.savefig(out_png, dpi=DPI, bbox_inches='tight', facecolor='white')
    out_svg = BASE_DIR / 'stomach_jaccard_heatmap.svg'
    fig.savefig(out_svg, format='svg', bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\nSaved to {BASE_DIR}")


if __name__ == '__main__':
    main()
