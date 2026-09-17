#!/usr/bin/env python3
"""
S2 panel C (S1 H until 2026-09-16, when the author split S1) - per-sample cell-type composition of the 38 non-stomach samples,
drawn at the size it prints at.

  printed panel  Supplementary Figure S2 C   (PROVENANCE.csv)
  predecessors   submission-tree/03_Final_Panels/10_Supplementaries/S1_QC_Annotation/
                 S1_H/create_S1_H_non_stomach_stacked_bar.py and its
                 de-identified copy in Supplementary_Fixes/S1_H (both a
                 12 x 4 in canvas cropped to its ink and fitted into the page)

Figure 1 C's drawing, copied and pointed at the non-stomach samples: the same
merge of T/NK, the same palette, the samples in alphabetical specimen order
named by their study IDs, the key in two columns at the right of the bars.

DRAWING READS A TABLE (cnsfig.cache): data/proportions.csv, written from the
full-dataset h5ad only when absent or with --recompute.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import FULL_DATASET_H5AD                       # noqa: E402
import panel_style_cns as style                           # noqa: E402
import _driver_base as base                               # noqa: E402
from cnsfig import cache, group_key                       # noqa: E402

FIG, PANEL = "S2_Cell_Annotation", "S2_C"
OUTPUT_DIR = HERE
# The printed box, millimetres: the page width; 38 bars at a 2.7 mm pitch
# and the two-column key.
PANEL_W_MM, PANEL_H_MM = 171.0, 58.0
#: The key stands at the RIGHT of the bars in two columns of six, since
#: 2026-09-15 (the author's fourth reading). It could not on 2026-09-15,
#: when C shared its row with the UMAP: beside the bars in a 96 mm slot it
#: left the 32 sample labels a 1.9 mm pitch, at which 6 pt labels
#: overprint, so it went under them in three columns. C now has a row to
#: itself (build_grid_v2.REPAGED): 161 mm, of which the key takes KEY_W_MM
#: ("Monocytes/Macrophages" 27.6 mm + "Endothelial cells" 17.9 mm + two
#: 1.3 mm circles, their pads and the column gap) and the bars 95 mm, a
#: 3.0 mm pitch. LABEL_H_MM is the rotated sample labels under the bars.
KEY_W_MM = 55.0
KEY_GAP_MM = 2.0
LABEL_H_MM = 12.0    # "P28-LN1_2" is 14 mm of rotated type
MARGIN = dict(left=9.0, right=KEY_W_MM + KEY_GAP_MM, top=4.5, bottom=LABEL_H_MM)

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
    obs = adata.obs[adata.obs['Sample site'] != 'Stomach'].copy()
    print(f"Non-stomach samples: {len(obs):,} cells, {obs['sample'].nunique()} samples")
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
    style.margins_mm(fig, **MARGIN)

    positions = np.arange(len(wide))
    bottom = np.zeros(len(wide))
    for ct in cell_type_order:
        ax.bar(positions, wide[ct].to_numpy(), bottom=bottom,
               color=CELL_TYPE_COLORS.get(ct, '#999999'), label=ct, width=0.8,
               edgecolor='white', linewidth=style.EDGE_PT)
        bottom += wide[ct].to_numpy()

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=90, ha='center')
    ax.set_xlim(-0.6, len(wide) - 0.4)
    ax.set_ylabel('Proportion')
    ax.set_ylim(0, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # The key, at fixed 1.3 mm circles (cnsfig.legend.group_key), in two
    # columns of six at the right of the bars, its top on the bars' top.
    box = ax.get_position()
    group_key(fig, [(ct, CELL_TYPE_COLORS.get(ct, '#999999')) for ct in cell_type_order],
              x_mm=PANEL_W_MM - KEY_W_MM,
              y_mm=(1.0 - box.y1) * PANEL_H_MM, ncol=2, columnspacing=1.2,
              linespacing=1.0)

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    base.save(fig, FIG, PANEL, "panel_S2_C")


if __name__ == '__main__':
    main()
