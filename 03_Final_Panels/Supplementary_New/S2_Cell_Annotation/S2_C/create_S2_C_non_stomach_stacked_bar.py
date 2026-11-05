#!/usr/bin/env python3
"""Supplementary Figure S2 C: per-sample cell-type composition of the 38
non-stomach samples (stacked bars, T/NK merged, samples labelled by study ID).

Reads data/proportions.csv; the full-dataset h5ad is opened only if that table
is missing or with --recompute. Writes the panel through _driver_base.save.
Run: python create_S2_C_non_stomach_stacked_bar.py [--recompute]"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import FULL_DATASET_H5AD                       # noqa: E402
import panel_style_cns as style                           # noqa: E402
import _driver_base as base                               # noqa: E402
from cnsfig import cache, group_key                       # noqa: E402

FIG, PANEL = "S2_Cell_Annotation", "S2_C"
# Panel size in mm; the key sits right of the bars in two columns.
PANEL_W_MM, PANEL_H_MM = 171.0, 58.0
KEY_W_MM = 55.0
KEY_GAP_MM = 2.0
LABEL_H_MM = 12.0    # room for the rotated sample labels
MARGIN = dict(left=9.0, right=KEY_W_MM + KEY_GAP_MM, top=4.5, bottom=LABEL_H_MM)

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
    """Per-sample cell-type proportions, samples named by their ST1 study ID
    and ordered alphabetically by specimen."""
    import scanpy as sc
    adata = sc.read_h5ad(FULL_DATASET_H5AD)
    obs = adata.obs[adata.obs['Sample site'] != 'Stomach'].copy()
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
    base.apply_style()
    long = cache.table(HERE, "proportions", compute_proportions)
    wide = long.pivot(index="order", columns="cell_type", values="proportion").sort_index()
    labels = long.drop_duplicates("order").sort_values("order")["sample"].tolist()
    # Most abundant first.
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
    base.despine(ax)

    # Key top aligned with the top of the bars.
    box = ax.get_position()
    group_key(fig, [(ct, CELL_TYPE_COLORS.get(ct, '#999999')) for ct in cell_type_order],
              x_mm=PANEL_W_MM - KEY_W_MM,
              y_mm=(1.0 - box.y1) * PANEL_H_MM, ncol=2, columnspacing=1.2,
              linespacing=1.0)

    base.save(fig, FIG, PANEL, "panel_S2_C")

if __name__ == '__main__':
    main()
