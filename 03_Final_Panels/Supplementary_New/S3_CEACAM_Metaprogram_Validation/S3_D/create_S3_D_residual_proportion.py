#!/usr/bin/env python3
"""Supplementary Figure S3 D: tumour-score-adjusted CEACAM5/6+ epithelial
proportion in pre-treatment primary tumours, R against NR.

Each sample's residual from regressing its CEACAM5/6+ proportion on its mean
tumour score; Mann-Whitney U, the panel prints the one-sided P (NR > R).
Reads data/residuals.csv; EPITHELIAL_DEPOSIT_H5AD is opened only if that table
is missing or with --recompute. Run: python create_S3_D_residual_proportion.py"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import EPITHELIAL_DEPOSIT_H5AD                 # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
from cnsfig import boxes                                  # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_D"
from shared.figure_config import PRE_PAIR  # noqa: E402
COLORS = {'Responsed': PRE_PAIR[0], 'No-response': PRE_PAIR[1]}
TITLE = 'Primary Tumor, Epi_CEACAM5/6\n(Tumor-score adjusted)'

# Panel size in mm.
W, H = 44.0, 44.0
MARGIN = dict(left=9.0, right=2.5, top=8.5, bottom=4.5)

def compute_residuals():
    import scanpy as sc
    from scipy import stats
    adata = sc.read_h5ad(EPITHELIAL_DEPOSIT_H5AD)
    pre = adata.obs['Treatment phase'] == 'Pre'
    valid = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
    sub = adata[pre & valid].copy()
    sub.obs['is_C2'] = (sub.obs['minor_cell_state'] == 'C2_Epi_CEACAM6').astype(int)
    df = sub.obs.groupby('sample', observed=True).agg(
        C2_proportion=('is_C2', lambda x: x.mean() * 100),
        mean_tumor_score=('tumor_score', 'mean'),
        group=('stomach_pre_grouping', 'first')).reset_index().dropna()
    x = df['mean_tumor_score'].to_numpy(dtype=float)
    y = df['C2_proportion'].to_numpy(dtype=float)
    slope, intercept, _, _, _ = stats.linregress(x, y)
    df['residual'] = y - (slope * x + intercept)
    R = df[df['group'] == 'Responsed']['residual'].to_numpy()
    NR = df[df['group'] == 'No-response']['residual'].to_numpy()
    # Both P values are cached; the panel prints the one-sided one.
    _, p = stats.mannwhitneyu(NR, R, alternative='two-sided')
    _, p_one = stats.mannwhitneyu(NR, R, alternative='greater')
    return pd.DataFrame({'group': df['group'].astype(str).to_numpy(),
                         'residual': df['residual'].to_numpy(dtype=float),
                         'p_two_sided': float(p),
                         'p_one_tailed': float(p_one)})

def draw(tab):
    base.apply_style()
    fig, ax = style.subplots_mm(W, H)
    style.margins_mm(fig, **MARGIN)
    R = tab[tab['group'] == 'Responsed']['residual'].to_numpy()
    NR = tab[tab['group'] == 'No-response']['residual'].to_numpy()
    p = float(tab['p_one_tailed'][0])      # one-sided, NR > R
    data = [R, NR]
    bxp = boxes.draw_boxes(ax, data, [1, 2],
                           [COLORS['Responsed'], COLORS['No-response']], width=0.6)
    y_max, y_min = max(R.max(), NR.max()), min(R.min(), NR.min())
    span = y_max - y_min
    _, p_text, _ = boxes.bracket(fig, ax, 1, 2, y_max, span, p, kind="pair",
                                 exact=True)
    ax.set_title(TITLE, pad=3)
    ax.set_ylabel('Residual proportion (%)')
    ax.set_xticks([1, 2]); ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_xlim(*boxes.box_xlim([1, 2], 0.6))
    ax.set_ylim(y_min - span * 0.15, y_max + span * 0.35)
    base.despine(ax)
    boxes.ylim_above(ax, p_text)
    boxes.assert_no_points(ax, bxp)
    return fig

if __name__ == "__main__":
    base.save(draw(cache.table(HERE, "residuals", compute_residuals)),
              FIG, PANEL, "panel_S3_D")
