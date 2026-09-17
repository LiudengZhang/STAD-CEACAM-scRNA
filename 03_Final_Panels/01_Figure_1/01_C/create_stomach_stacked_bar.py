#!/usr/bin/env python3
"""
Figure 1 panel C - cell-type composition of each stomach sample, as upright
stacked bars, one sample a column, with the twelve types keyed BELOW the bars.

  printed panel  Figure 1 C       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

DRAWN AT THE SIZE IT PRINTS, since 2026-09-15 - see 01_B/create_stomach_umap.py
for why Figure 1 is now assembled into slots on the 171.10 mm page. The type
is set by 00_Config/panel_style_cns.py; the sample labels are the study IDs
of Supplementary Table 1, read from the h5ad's own 'Sample ID' column (the
earlier script joined them from a submission-tree copy of ST1 through a column that
table no longer carries).

THE KEY STANDS UNDER THE BARS, since 2026-09-16 (the author's fifth reading:
"1C: put the key below after all, it doesn't fit; 1B and 1C should be the
same height"). The history, so it is not re-argued: the fourth reading
(2026-09-15) asked for B and C on one row and for C's key at the right in two
columns of six; at 6 pt that 55 mm key, the 9 mm y axis and 32 rotated
sample labels at a 2.64 mm pitch left the UMAP 15 mm, so for one night the
bars ran across the panel, one sample a row, and the row was 96 mm tall
against a 56 mm cloud. With the key below, the upright bars fit: 32 columns
at a 2.6 mm pitch in an 83 mm plotting box, the 6 pt sample labels rotated
90 degrees under them, the key in KEY_NCOL columns under the labels, and the
panel is 62 mm - the same height as B (build_grid_v2.REPAGED). The same 32
samples in the same order left to right, the same twelve types in the same
order bottom to top, the same proportions (data/proportions.csv, held by
freeze_baseline.py, is untouched); the round-39 upright script (fcc9192) is
the drawing this returns to, and compare_panel_content.py holds them equal.

DRAWING READS A TABLE (cnsfig.cache). data/proportions.csv holds one row per
sample and cell type - sample (the ST1 study ID), cell_type, proportion -
written from the full-dataset h5ad only when the table is absent or
--recompute is passed. The merge of T/NK, the renamings, the alphabetical
sample order and the abundance order of the types are the earlier script's.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FULL_DATASET_H5AD                       # noqa: E402
import panel_style_cns as style                           # noqa: E402
import slots                                              # noqa: E402
from cnsfig import cache, group_key                       # noqa: E402

PANEL_LETTER = "C"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(1, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(1, PANEL_LETTER)
OUTPUT_DIR = Path(__file__).resolve().parent
#: The key stands UNDER the bars, since 2026-09-16 (the author's fifth
#: reading). LABEL_W_MM is the y axis at the left ("Proportion" and its
#: ticks); LABEL_H_MM the rotated 6 pt sample labels under the bars ("P15-P1"
#: is 6.9 mm long, plus the tick and its pad); KEY_GAP_MM the paper between
#: the lowest label ink and the key's top; KEY_NCOL the key's columns -
#: the most that fit 93 mm at 6 pt with 1 mm between the longest label and
#: the next column's circle (measured after drawing, below). The bars take
#: what is left of the 62 mm. SLOT_CLEAR is the paper between the outer bars
#: and the y spine, in slots (0.6 x 2.6 mm = 1.56 mm; check_box_clearance
#: holds every bar 1 mm clear of a vertical rule).
LABEL_W_MM = 9.0
LABEL_H_MM = 8.5
KEY_GAP_MM = 1.5
KEY_NCOL = 4
KEY_LINESPACING = 0.5
KEY_BOTTOM_MM = 0.5
SLOT_CLEAR = 0.6
MARGIN = dict(left=LABEL_W_MM, right=1.0, top=4.5)

# Set2 + Set3 palette (12 major cell types, consistent across Fig 1B & 1C)
CELL_TYPE_COLORS = {
    'Epithelial cells': '#e5c494',
    'T/NK cells': '#66c2a5',
    'Monocytes/Macrophages': '#fc8d62',
    'Plasma cells': '#8da0cb',
    'B cells': '#ffd92f',
    'Endothelial cells': '#a6d854',
    'Fibroblasts': '#e78ac3',
    'Neutrophils': '#b3b3b3',
    'Mast cells': '#fb8072',
    'Pericytes': '#bebada',
    'Dendritic cells': '#80b1d3',
    'Hepatocytes': '#bc80bd',
}

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
    """The computing half: per-sample cell-type proportions, samples named
    by their ST1 study ID and ordered as the earlier drawing ordered them
    (alphabetically by specimen)."""
    import scanpy as sc
    print(f"Loading data from: {FULL_DATASET_H5AD}")
    adata = sc.read_h5ad(FULL_DATASET_H5AD)
    print(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    obs = adata.obs[adata.obs['Sample site'] == 'Stomach'].copy()
    print(f"Stomach samples: {len(obs):,} cells, {obs['sample'].nunique()} samples")
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
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt")

    long = cache.table(OUTPUT_DIR, "proportions", compute_proportions)
    wide = long.pivot(index="order", columns="cell_type", values="proportion").sort_index()
    labels = long.drop_duplicates("order").sort_values("order")["sample"].tolist()
    # Most abundant first, as before.
    cell_type_order = wide.mean().sort_values(ascending=False).index.tolist()
    wide = wide[cell_type_order]

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    key_entries = [(ct, CELL_TYPE_COLORS.get(ct, '#999999')) for ct in cell_type_order]

    # The key is drawn first, at the canvas's bottom-left, and measured; the
    # bars take the paper above it. A key is a fixed thing (12 names at 6 pt
    # in KEY_NCOL columns) and the bars are not, so the bars give.
    probe = group_key(fig, key_entries, x_mm=LABEL_W_MM, y_mm=0.0,
                      ncol=KEY_NCOL, columnspacing=1.2, linespacing=KEY_LINESPACING)
    fig.canvas.draw()
    key_h_mm = probe.get_window_extent().height / fig.dpi * 25.4
    probe.remove()

    key_top_mm = PANEL_H_MM - KEY_BOTTOM_MM - key_h_mm       # from the canvas top
    axes_bottom_mm = key_top_mm - KEY_GAP_MM - LABEL_H_MM
    axes_h = axes_bottom_mm - MARGIN["top"]
    if axes_h < 25.0:
        raise RuntimeError(f"the bars would be {axes_h:.1f} mm tall under a "
                           f"{key_h_mm:.1f} mm key; the slot is {PANEL_H_MM} mm")
    style.margins_mm(fig, **dict(MARGIN, bottom=PANEL_H_MM - axes_bottom_mm))

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

    # The key, at fixed 1.3 mm circles (cnsfig.legend.group_key), under the
    # sample labels, its left on the plotting box's left.
    leg = group_key(fig, key_entries, x_mm=LABEL_W_MM, y_mm=key_top_mm,
                    ncol=KEY_NCOL, columnspacing=1.2, linespacing=KEY_LINESPACING)

    # Measure what was declared: the sample pitch, the key's right edge, the
    # paper between the lowest label and the key, and the paper between each
    # label and the next column's circle.
    fig.canvas.draw()
    px_mm = 25.4 / fig.dpi
    pitch_mm = (ax.get_window_extent().width * px_mm) / (len(wide) - 1 + 2 * SLOT_CLEAR)
    if pitch_mm < 2.45:
        raise RuntimeError(f"sample pitch {pitch_mm:.2f} mm < 2.45 mm: rotated 6 pt "
                           f"labels would touch")
    leg_box = leg.get_window_extent()
    if leg_box.x1 * px_mm > PANEL_W_MM - 1.0:
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
    print(f"  bars {axes_h:.1f} mm tall, pitch {pitch_mm:.2f} mm; key {KEY_NCOL} columns, "
          f"{leg_box.width * px_mm:.1f} x {key_h_mm:.1f} mm, {gap:.2f} mm under the "
          f"labels, {clear:.2f} mm between columns")

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, OUTPUT_DIR / '01_C_sample_balance_stomach')
    print(f"Saved: 01_C_sample_balance_stomach.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
