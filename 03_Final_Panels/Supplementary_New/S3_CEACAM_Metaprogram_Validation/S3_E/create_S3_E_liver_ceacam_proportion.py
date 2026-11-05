#!/usr/bin/env python3
"""Supplementary Figure S3 E: CEACAM5/6+ epithelial proportion in liver
metastases, pre-treatment R against NR (two-sided Mann-Whitney U).

Reads data/proportions.csv; EPITHELIAL_H5AD is opened only if that table is
missing or with --recompute. Run: python create_S3_E_liver_ceacam_proportion.py"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import EPITHELIAL_H5AD                         # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
from cnsfig import boxes                                  # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_E"
from shared.figure_config import PRE_PAIR  # noqa: E402
COLORS = {'R': PRE_PAIR[0], 'NR': PRE_PAIR[1]}
TITLE = 'Liver Met., Epi_CEACAM5/6'

# Panel size in mm.
W, H = 40.0, 44.0
MARGIN = dict(left=9.0, right=1.0, top=8.5, bottom=4.5)

def compute_proportions():
    import scanpy as sc
    from scipy.stats import mannwhitneyu
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    adata = adata[adata.obs['Sample site'] == 'Liver'].copy()
    adata.obs['response'] = adata.obs['liver_pre_grouping'].map(
        {'Responsed': 'R', 'No-response': 'NR'}).astype(str)
    adata = adata[adata.obs['response'].isin(['R', 'NR'])].copy()
    st = adata.obs['minor_cell_state']
    c2 = st.str.contains('C2_', na=False) & st.str.contains('CEACAM', na=False)
    if c2.sum() == 0:
        c2 = st.str.contains('CEACAM', na=False)
    adata.obs['is_c2'] = c2.astype(int)
    df = adata.obs.groupby(['sample', 'response'], observed=True).agg(
        n_total=('is_c2', 'size'), n_c2=('is_c2', 'sum')).reset_index()
    df['proportion'] = df['n_c2'] / df['n_total'] * 100
    R = df[df['response'] == 'R']['proportion'].to_numpy(dtype=float)
    NR = df[df['response'] == 'NR']['proportion'].to_numpy(dtype=float)
    _, p = mannwhitneyu(R, NR, alternative='two-sided')
    return pd.DataFrame({'group': df['response'].astype(str).to_numpy(),
                         'proportion': df['proportion'].to_numpy(dtype=float),
                         'p_two_sided': float(p)})

def draw(tab):
    base.apply_style()
    fig, ax = style.subplots_mm(W, H)
    style.margins_mm(fig, **MARGIN)
    R = tab[tab['group'] == 'R']['proportion'].to_numpy()
    NR = tab[tab['group'] == 'NR']['proportion'].to_numpy()
    p = float(tab['p_two_sided'][0])
    bxp = boxes.draw_boxes(ax, [R, NR], [1, 2], [COLORS['R'], COLORS['NR']],
                           width=0.6)
    y_max = float(np.concatenate([R, NR]).max())
    _, p_text, _ = boxes.bracket(fig, ax, 1, 2, y_max, y_max, p, kind="pair")
    ax.set_title(TITLE, pad=3)
    ax.set_ylabel('Proportion (%)')
    ax.set_xticks([1, 2]); ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_xlim(*boxes.box_xlim([1, 2], 0.6))
    ax.set_ylim(0, y_max * 1.40)
    base.despine(ax)
    boxes.ylim_above(ax, p_text)
    boxes.assert_no_points(ax, bxp)
    return fig

if __name__ == "__main__":
    base.save(draw(cache.table(HERE, "proportions", compute_proportions)),
              FIG, PANEL, "panel_S3_E")
