#!/usr/bin/env python3
"""
S2 panel E - the CEACAM5/6+ epithelial proportion in the liver metastases,
R against NR, drawn at the size it prints at.

  printed panel  Supplementary Figure S2 E   (PROVENANCE.csv; the submission-tree
                 directory was S3_I - the letter is looked up, not inferred)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/
                 S3_CEACAM_Metaprogram_Validation/S3_I/
                 create_S2_I_liver_c2ceacam_proportion.py

The title is the printed page's ("Liver Met., Epi_CEACAM5/6") and the P is
printed as a value through panel_style_cns.p_text_kw (the page prints
"p = 0.19"; the predecessor set "ns"); both declared in labels.py
RENAMES_SUPPLEMENTARY.

DRAWING READS A TABLE (cnsfig.cache): data/proportions.csv holds each liver
sample's CEACAM5/6+ proportion and group and the Mann-Whitney P; the
epithelial h5ad is opened only when the table is absent or with --recompute.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import EPITHELIAL_H5AD                         # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
from cnsfig import boxes                                  # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_E"
COLORS = {'R': '#0072B2', 'NR': '#D55E00'}
TITLE = 'Liver Met., Epi_CEACAM5/6'

# The printed box, millimetres: the last of four in row 2.
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
    # ONE BOX, ONE BRACKET (cnsfig.boxes, 2026-09-16): Figure 2 H/I's box
    # replaces the group-coloured medians and 2.2 pt fliers; the label's ink
    # sits 0.4 mm above the bracket. The y range is [0, y_max], so the
    # bracket is lifted by a tenth of y_max.
    bxp = boxes.draw_boxes(ax, [R, NR], [1, 2], [COLORS['R'], COLORS['NR']],
                           width=0.6)
    y_max = float(np.concatenate([R, NR]).max())
    _, p_text, _ = boxes.bracket(fig, ax, 1, 2, y_max, y_max, p, kind="pair")
    ax.set_title(TITLE, pad=3)
    ax.set_ylabel('Proportion (%)')
    ax.set_xticks([1, 2]); ax.set_xticklabels(['R', 'NR'])
    ax.set_xlim(*boxes.box_xlim([1, 2], 0.6))
    ax.set_ylim(0, y_max * 1.40)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    boxes.ylim_above(ax, p_text)
    boxes.assert_no_points(ax, bxp)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"ink outside the {W} x {H} mm canvas: {over}")
    return fig


def main():
    tab = cache.table(HERE, "proportions", compute_proportions)
    base.save(draw(tab), FIG, PANEL, "panel_S3_E")
    return 0


if __name__ == "__main__":
    sys.exit(main())
