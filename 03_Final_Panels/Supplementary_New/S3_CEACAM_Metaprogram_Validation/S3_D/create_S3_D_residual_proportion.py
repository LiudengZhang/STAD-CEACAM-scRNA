#!/usr/bin/env python3
"""
S2 panel D - the purity-adjusted primary-tumour CEACAM5/6+ epithelial
proportion, R against NR (pre-treatment), drawn at the size it prints at.

  printed panel  Supplementary Figure S3 D   (PROVENANCE.csv)
  predecessor    03_Final_Panels/02_Figure_2/02_E2/
                 create_c2_proportion_adjusted.py, whose second figure is the
                 printed panel (its first, the tumour-score scatter, was never
                 printed). The submitted panel printed the one-sided P = 0.03;
                 the two-sided P = 0.057 (R1.3c) has shipped since 2026-09-08
                 (Supplementary_Fixes/patch_S2_two_sided.py, now retired) and
                 is what this panel prints.

The title is the printed page's ("Primary Tumor, Epi_CEACAM5/6 (Tumor-score
adjusted)"), declared in labels.py RENAMES_SUPPLEMENTARY against the
predecessor's; beside panel E ("Liver Met.") the tissue has to be named.

DRAWING READS A TABLE (cnsfig.cache): data/residuals.csv holds each
pre-treatment sample's residual after regressing the CEACAM5/6+ proportion
on the mean tumour score, its group, and the Mann-Whitney P; the epithelial
h5ad is opened only when the table is absent or with --recompute.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import EPITHELIAL_DEPOSIT_H5AD                 # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
from cnsfig import boxes                                  # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_D"
COLORS = {'Responsed': '#0072B2', 'No-response': '#D55E00'}
TITLE = 'Primary Tumor, Epi_CEACAM5/6\n(Tumor-score adjusted)'

# The printed box, millimetres: one of four in row 2.
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
    _, p = stats.mannwhitneyu(NR, R, alternative='two-sided')
    # The specimen code is not written; the order within a group is the
    # predecessor's (groupby order), which the jitter seed depends on.
    return pd.DataFrame({'group': df['group'].astype(str).to_numpy(),
                         'residual': df['residual'].to_numpy(dtype=float),
                         'p_two_sided': float(p)})


def draw(tab):
    base.apply_style()
    fig, ax = style.subplots_mm(W, H)
    style.margins_mm(fig, **MARGIN)
    R = tab[tab['group'] == 'Responsed']['residual'].to_numpy()
    NR = tab[tab['group'] == 'No-response']['residual'].to_numpy()
    p = float(tab['p_two_sided'][0])
    data = [R, NR]
    # ONE BOX, ONE BRACKET, NO POINTS (cnsfig.boxes, 2026-09-16): the
    # predecessor's group-coloured medians, 2.2 pt fliers and the jittered
    # sample points (rng 0) are gone; the boxes are Figure 2 H/I's and the
    # label's ink sits 0.4 mm above the bracket.
    bxp = boxes.draw_boxes(ax, data, [1, 2],
                           [COLORS['Responsed'], COLORS['No-response']], width=0.6)
    y_max, y_min = max(R.max(), NR.max()), min(R.min(), NR.min())
    span = y_max - y_min
    _, p_text, _ = boxes.bracket(fig, ax, 1, 2, y_max, span, p, kind="pair")
    ax.set_title(TITLE, pad=3)
    ax.set_ylabel('Residual proportion (%)')
    ax.set_xticks([1, 2]); ax.set_xticklabels(['R', 'NR'])
    ax.set_xlim(*boxes.box_xlim([1, 2], 0.6))
    ax.set_ylim(y_min - span * 0.15, y_max + span * 0.35)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    boxes.ylim_above(ax, p_text)
    boxes.assert_no_points(ax, bxp)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"ink outside the {W} x {H} mm canvas: {over}")
    return fig


def main():
    tab = cache.table(HERE, "residuals", compute_residuals)
    base.save(draw(tab), FIG, PANEL, "panel_S3_D")
    return 0


if __name__ == "__main__":
    sys.exit(main())
