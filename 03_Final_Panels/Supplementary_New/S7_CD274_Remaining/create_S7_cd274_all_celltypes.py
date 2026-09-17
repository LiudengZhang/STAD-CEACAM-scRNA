#!/usr/bin/env python3
"""
S6 panels A-I - sample-level PD-L1 (CD274) expression, post-treatment R
against NR, in the nine cell types not shown in Figure 5, drawn at the size
they print at.

  printed panels  Supplementary Figure S7 A-I   (PROVENANCE.csv; one row per
                  panel, each pointing at this script)
  predecessor     submission-tree/03_Final_Panels/10_Supplementaries/S6_CD274_Remaining/
                  create_S6_cd274_all_celltypes.py  (nine 3.5 x 5 cm canvases
                  at 4x)

One script draws the nine panels, as before, each into its own directory
S7_<letter>/. Since 2026-09-16 the boxes and brackets are cnsfig.boxes' (one
box style through the paper, no individual points - the author's ruling of
that day) and the P labels go through panel_style_cns.p_text_kw like every
other pairwise bracket: a star below 0.05, the two-decimal value at or above
it. The predecessor and the printed page set "ns" there; the legend says
which.

DRAWING READS A TABLE (cnsfig.cache): data/sample_means.csv holds, per cell
type, each post-treatment sample's mean CD274 (from .raw, samples with at
least MIN_CELLS cells) and its group, and data/tests.csv the Mann-Whitney P
per cell type; the lineage h5ads and the full dataset are opened only when
the tables are absent or with --recompute.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parent / "_drivers"))
from paths import (FULL_DATASET_H5AD, TCD8_H5AD, TCD4_H5AD, NK_CELLS_H5AD,   # noqa: E402
                   B_CELLS_H5AD, NEUTROPHILS_H5AD, ENDOTHELIAL_H5AD, PERICYTE_H5AD)
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
from cnsfig import boxes                                  # noqa: E402
from cnsfig.boxes import box_xlim                         # noqa: E402
import _driver_base as base                               # noqa: E402

FIG = "S7_CD274_Remaining"
MIN_CELLS = 20
COLOR_R, COLOR_NR = '#ffcfd2', '#f1c0e8'
PANELS = [
    ('A', TCD8_H5AD, 'CD8+ T cells', 'CD8+ T cells'),
    ('B', TCD4_H5AD, 'CD4+ T cells', 'CD4+ T cells'),
    ('C', NK_CELLS_H5AD, 'NK cells', 'NK cells'),
    ('D', B_CELLS_H5AD, 'B cells', 'B cells'),
    ('E', None, 'Plasma cells', 'Plasma cells'),
    ('F', None, 'Mast cells', 'Mast cells'),
    ('G', NEUTROPHILS_H5AD, 'Neutrophils', 'Neutrophils'),
    ('H', ENDOTHELIAL_H5AD, 'Endothelial cells', 'Endothelial cells'),
    ('I', PERICYTE_H5AD, 'Pericyte', 'Pericytes'),
]

# The printed box, millimetres: five to a row.
W, H = 30.0, 46.0
MARGIN = dict(left=12.0, right=1.0, top=8.0, bottom=7.0)   # left 11 until 2026-09-16 (D's 0.15 tick)


def _sample_means(adata):
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata = adata[adata.obs['Treatment phase'] == 'Post'].copy()
    adata.obs['response'] = adata.obs['stomach_post_grouping'].map(
        {'Responsed': 'R', 'No-response': 'NR'}).astype(str)
    adata = adata[adata.obs['response'].isin(['R', 'NR'])].copy()
    gene = 'CD274'
    if adata.raw is not None and gene in adata.raw.var_names:
        expr = adata.raw.X[:, list(adata.raw.var_names).index(gene)]
    elif gene in adata.var_names:
        expr = adata.X[:, list(adata.var_names).index(gene)]
    else:
        return pd.DataFrame(columns=['response', 'mean'])
    expr = expr.toarray().flatten() if hasattr(expr, 'toarray') else np.asarray(expr).flatten()
    df = pd.DataFrame({'sample': adata.obs['sample'].astype(str).to_numpy(),
                       'response': adata.obs['response'].to_numpy(), 'cd274_expr': expr})
    counts = df.groupby('sample', observed=True).size()
    df = df[df['sample'].isin(counts[counts >= MIN_CELLS].index)]
    out = df.groupby(['sample', 'response'], observed=True)['cd274_expr'].mean().reset_index()
    return out.rename(columns={'cd274_expr': 'mean'})[['response', 'mean']]


def compute_means():
    import scanpy as sc
    full = None
    rows = []
    for letter, path, cell_type, _ in PANELS:
        if path is not None and path.exists():
            adata = sc.read_h5ad(path)
        else:
            if full is None:
                full = sc.read_h5ad(FULL_DATASET_H5AD)
            adata = full[full.obs['major_cell_type'] == cell_type].copy()
        means = _sample_means(adata)
        for r in means.itertuples():
            rows.append({'panel': letter, 'group': r.response, 'mean': float(r.mean)})
    return pd.DataFrame(rows)


def compute_tests():
    from scipy.stats import mannwhitneyu
    tab = cache.table(HERE, "sample_means", compute_means)
    rows = []
    for letter, _, _, _ in PANELS:
        sub = tab[tab['panel'] == letter]
        r = sub[sub['group'] == 'R']['mean'].to_numpy()
        nr = sub[sub['group'] == 'NR']['mean'].to_numpy()
        p = mannwhitneyu(nr, r, alternative='two-sided')[1] if len(r) >= 2 and len(nr) >= 2 else np.nan
        rows.append({'panel': letter, 'p_two_sided': float(p)})
    return pd.DataFrame(rows)


def draw(letter, display_name, means, tests):
    fig, ax = style.subplots_mm(W, H)
    style.margins_mm(fig, **MARGIN)
    sub = means[means['panel'] == letter]
    r = sub[sub['group'] == 'R']['mean'].to_numpy()
    nr = sub[sub['group'] == 'NR']['mean'].to_numpy()
    p = float(tests[tests['panel'] == letter]['p_two_sided'].iloc[0])
    # ONE BOX, ONE BRACKET, NO POINTS (cnsfig.boxes, 2026-09-16). Until then:
    # 0.7-alpha faces, the samples jittered on top (rng 42), a local star
    # ladder printing "ns".
    if np.isnan(p):
        raise SystemExit(f"S6 {letter}: no P (a group has fewer than two samples)")
    bxp = boxes.draw_boxes(ax, [r, nr], [0, 1], [COLOR_R, COLOR_NR], width=0.5)
    allv = np.concatenate([r, nr]) if len(r) and len(nr) else np.array([0.0])
    y_max = allv.max()
    y_rng = (allv.max() - allv.min()) if len(allv) > 1 else 0.1
    y_rng = y_rng or 0.1
    _, p_text, _ = boxes.bracket(fig, ax, 0, 1, y_max, y_rng, p, kind="pair",
                                 lift=0.10, arm=0.02)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Post-R\n(n={len(r)})", f"Post-NR\n(n={len(nr)})"])
    ax.set_xlim(*box_xlim([0, 1], 0.5))
    ax.set_ylabel('PD-L1 (CD274)\nExpression')
    ax.set_title(f'PD-L1 (CD274)\n{display_name}', pad=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    boxes.ylim_above(ax, p_text)
    boxes.assert_no_points(ax, bxp)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"S6 {letter}: ink outside the {W} x {H} mm canvas: {over}")
    return fig


def main():
    base.apply_style()
    means = cache.table(HERE, "sample_means", compute_means)
    tests = cache.table(HERE, "tests", compute_tests)
    for letter, _, _, display_name in PANELS:
        base.save(draw(letter, display_name, means, tests), FIG, f"S7_{letter}", f"panel_S7_{letter}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
