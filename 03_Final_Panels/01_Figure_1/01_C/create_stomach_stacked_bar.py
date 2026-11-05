#!/usr/bin/env python3
"""
Figure 1C: cell-type composition of each stomach sample as stacked bars, one
sample per column, with the cell-type key below the bars.

Reads data/proportions.csv (sample, cell_type, proportion); the table is
computed from the full h5ad when absent or with --recompute.
Writes 01_C_sample_balance_stomach.[svg|pdf|png] beside this script.

Run: python create_stomach_stacked_bar.py [--recompute]

Samples are labelled with their Supplementary Table 1 study IDs and ordered
alphabetically by specimen; cell types are stacked most abundant first.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FULL_DATASET_H5AD                       # noqa: E402
from cnsfig import cache, group_key, Panel, style          # noqa: E402

PANEL = Panel(1, "C", __file__, small_pt=6.0)
# Layout in mm: room for the y axis (LABEL_W_MM) and the rotated sample labels
# (LABEL_H_MM), the gap above the key, and the key's columns. The bars take the
# remaining height. SLOT_CLEAR is the space beside the outer bars, in bar slots.
LABEL_W_MM = 9.0
LABEL_H_MM = 8.5
KEY_GAP_MM = 1.5
KEY_NCOL = 4
KEY_LINESPACING = 0.5
KEY_BOTTOM_MM = 0.5
SLOT_CLEAR = 0.6
MARGIN = dict(left=LABEL_W_MM, right=1.0, top=4.5)

# Major-cell-type palette shared with Figure 1B.
from shared.figure_config import CELL_TYPE_COLORS  # noqa: E402,F401

MERGE = {
    'CD4+ T cells': 'T/NK cells',
    'CD8+ T cells': 'T/NK cells',
    'NK cells': 'T/NK cells',
    'DC cells': 'Dendritic cells',
    'Fibroblast': 'Fibroblasts',
    'Pericyte': 'Pericytes',
    'Hepatocyte': 'Hepatocytes',
}


def compute_proportions():
    """Per-sample cell-type proportions, labelled by study ID and ordered
    alphabetically by specimen."""
    import scanpy as sc
    adata = sc.read_h5ad(FULL_DATASET_H5AD)
    obs = adata.obs[adata.obs['Sample site'] == 'Stomach'].copy()
    obs['cell_type'] = obs['major_cell_type'].astype(str).replace(MERGE)

    proportions = obs.groupby(['sample', 'cell_type'], observed=True).size().unstack(fill_value=0)
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    if 'Sample ID' not in obs.columns:
        raise SystemExit("the h5ad carries no 'Sample ID' (study ID) column; "
                         "the x labels would print specimen numbers")
    orig_to_sample = (obs.groupby('sample', observed=True)['Sample ID']
                      .agg(lambda v: v.astype(str).iloc[0]).to_dict())
    sample_order = sorted(proportions.index.tolist())
    rows = []
    for i, s in enumerate(sample_order):
        for ct in proportions.columns:
            rows.append({"order": i, "sample": orig_to_sample[s], "cell_type": ct,
                         "proportion": float(proportions.loc[s, ct])})
    return pd.DataFrame(rows)


def main():
    long = cache.table(PANEL.dir, "proportions", compute_proportions)
    wide = long.pivot(index="order", columns="cell_type", values="proportion").sort_index()
    labels = long.drop_duplicates("order").sort_values("order")["sample"].tolist()
    # Most abundant first.
    cell_type_order = wide.mean().sort_values(ascending=False).index.tolist()
    wide = wide[cell_type_order]

    fig, ax = PANEL.subplots()
    key_entries = [(ct, CELL_TYPE_COLORS.get(ct, '#999999')) for ct in cell_type_order]

    # Draw and measure the key first; the bars take the space above it.
    probe = group_key(fig, key_entries, x_mm=LABEL_W_MM, y_mm=0.0,
                      ncol=KEY_NCOL, columnspacing=1.2, linespacing=KEY_LINESPACING)
    fig.canvas.draw()
    key_h_mm = probe.get_window_extent().height / fig.dpi * 25.4
    probe.remove()

    key_top_mm = PANEL.h_mm - KEY_BOTTOM_MM - key_h_mm       # from the canvas top
    axes_bottom_mm = key_top_mm - KEY_GAP_MM - LABEL_H_MM
    axes_h = axes_bottom_mm - MARGIN["top"]
    if axes_h < 25.0:
        raise RuntimeError(f"the bars would be {axes_h:.1f} mm tall under a "
                           f"{key_h_mm:.1f} mm key; the slot is {PANEL.h_mm} mm")
    style.margins_mm(fig, **dict(MARGIN, bottom=PANEL.h_mm - axes_bottom_mm))

    positions = np.arange(len(wide))
    bottom = np.zeros(len(wide))
    for ct in cell_type_order:
        ax.bar(positions, wide[ct].to_numpy(), bottom=bottom,
               color=CELL_TYPE_COLORS.get(ct, '#999999'), label=ct, width=0.8,
               edgecolor='white', linewidth=style.EDGE_PT)
        bottom += wide[ct].to_numpy()

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=90, ha='center')
    ax.set_xlim(-SLOT_CLEAR, len(wide) - 1 + SLOT_CLEAR)
    ax.set_ylabel('Proportion')
    ax.set_ylim(0, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # The key, under the sample labels, aligned with the plotting box.
    leg = group_key(fig, key_entries, x_mm=LABEL_W_MM, y_mm=key_top_mm,
                    ncol=KEY_NCOL, columnspacing=1.2, linespacing=KEY_LINESPACING)

    # Check the drawn layout: sample pitch, key width, label-to-key gap and
    # spacing between key columns.
    fig.canvas.draw()
    px_mm = 25.4 / fig.dpi
    pitch_mm = (ax.get_window_extent().width * px_mm) / (len(wide) - 1 + 2 * SLOT_CLEAR)
    if pitch_mm < 2.45:
        raise RuntimeError(f"sample pitch {pitch_mm:.2f} mm < 2.45 mm: rotated 6 pt "
                           f"labels would touch")
    leg_box = leg.get_window_extent()
    if leg_box.x1 * px_mm > PANEL.w_mm - 1.0:
        raise RuntimeError(f"a {KEY_NCOL}-column key is {leg_box.width * px_mm:.1f} mm "
                           f"wide; use fewer columns")
    label_bottom = min(t.get_window_extent().y0 for t in ax.get_xticklabels())
    gap = (label_bottom - leg_box.y1) * px_mm
    if gap < KEY_GAP_MM - 0.05:
        raise RuntimeError(f"{gap:.2f} mm between the sample labels and the key")
    handles = getattr(leg, "legend_handles", None) or leg.legendHandles
    tb = [t.get_window_extent() for t in leg.get_texts()]
    hb = [h.get_window_extent() for h in handles]
    # A label's right edge to the next column's circle on the same row.
    clear = min((h.x0 - t.x1) * px_mm for t in tb for h in hb
                if h.x0 > t.x1 and abs(h.y0 - t.y0) * px_mm < 1.5)
    if clear < 1.0:
        raise RuntimeError(f"{clear:.2f} mm between a key label and the next column")

    PANEL.finish(fig, '01_C_sample_balance_stomach')


if __name__ == '__main__':
    main()
