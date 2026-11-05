#!/usr/bin/env python3
"""Supplementary Figure S10 A-I: per-sample mean PD-L1 (CD274) expression,
post-treatment R vs NR, in nine cell types not shown in Figure 5.

Reads data/sample_means.csv and data/tests.csv (two-sided Mann-Whitney P per
cell type); with --recompute or when these are missing, rebuilds them from the
per-cell-type h5ads and the full dataset. Writes one panel per cell type,
S10_A ... S10_I.

Run: python create_S10_cd274_all_celltypes.py [--recompute]

"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / "_drivers"),
                str(Path(__file__).resolve().parents[3] / "00_Config")]
from paths import (FULL_DATASET_H5AD, TCD8_H5AD, TCD4_H5AD, NK_CELLS_H5AD,   # noqa: E402
                   B_CELLS_H5AD, NEUTROPHILS_H5AD, ENDOTHELIAL_H5AD,
                   PERICYTE_H5AD, SOUND_SOURCE)
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"
                       / "shared"))
from expression import expression_source                  # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
from cnsfig import boxes                                  # noqa: E402
from cnsfig.boxes import box_xlim                         # noqa: E402
import _driver_base as base                               # noqa: E402

FIG = "S10_CD274_Remaining"
MIN_CELLS = 20


def source(path):
    """Return the file to read for `path`, as mapped by paths.SOUND_SOURCE.

    """
    if path is None:
        return None
    return SOUND_SOURCE.get(Path(path).name, path)
from shared.figure_config import POST_PAIR  # noqa: E402
COLOR_R, COLOR_NR = POST_PAIR          # post-treatment: the saturated pair
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

# Panel size and margins in mm.
W, H = 30.0, 46.0
MARGIN = dict(left=12.0, right=1.0, top=8.0, bottom=7.0)

def _sample_means(adata, path):
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata = adata[adata.obs['Treatment phase'] == 'Post'].copy()
    adata.obs['response'] = adata.obs['stomach_post_grouping'].map(
        {'Responsed': 'R', 'No-response': 'NR'}).astype(str)
    adata = adata[adata.obs['response'].isin(['R', 'NR'])].copy()
    gene = 'CD274'
    # expression_source raises if it cannot resolve the matrix; so does a
    # missing gene.
    src = expression_source(adata, path)
    if gene not in src.var_names:
        raise RuntimeError(
            f"{gene} is not in this object's expression matrix. Refusing to "
            "return an empty table: a panel with no data must fail, not draw "
            "blank.")
    expr = src.X[:, list(src.var_names).index(gene)]
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
        # E and F have no object of their own and are cut from the full dataset.
        read = source(path)
        if read is not None:
            if not read.exists():
                raise FileNotFoundError(
                    f"panel {letter} declares {read} and it is not there. "
                    "Refusing to fall back to full_dataset.h5ad: its .raw is "
                    "not log1p CP10K, so the fallback would silently swap in "
                    "an incorrectly normalised matrix.")
            adata = sc.read_h5ad(read)
        else:
            if full is None:
                read_full = source(FULL_DATASET_H5AD)
                full = sc.read_h5ad(read_full)
            adata = full[full.obs['major_cell_type'] == cell_type].copy()
            read = source(FULL_DATASET_H5AD)
        means = _sample_means(adata, read)
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
    if np.isnan(p):
        raise SystemExit(f"S7 {letter}: no P (a group has fewer than two samples)")
    bxp = boxes.draw_boxes(ax, [r, nr], [0, 1], [COLOR_R, COLOR_NR], width=0.5)
    allv = np.concatenate([r, nr]) if len(r) and len(nr) else np.array([0.0])
    y_max = allv.max()
    y_rng = (allv.max() - allv.min()) if len(allv) > 1 else 0.1
    y_rng = y_rng or 0.1
    _, p_text, _ = boxes.bracket(fig, ax, 0, 1, y_max, y_rng, p, kind="pair",
                                 lift=0.10, arm=0.02)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Post-R", "Post-NR"])
    ax.set_xlim(*box_xlim([0, 1], 0.5))
    ax.set_ylabel('PD-L1 (CD274)\nExpression')
    ax.set_title(f'PD-L1 (CD274)\n{display_name}', pad=3)
    base.despine(ax)
    boxes.ylim_above(ax, p_text)
    boxes.assert_no_points(ax, bxp)
    return fig

if __name__ == "__main__":
    base.apply_style()
    means = cache.table(HERE, "sample_means", compute_means)
    tests = cache.table(HERE, "tests", compute_tests)
    for letter, _, _, display_name in PANELS:
        base.save(draw(letter, display_name, means, tests), FIG, f"S10_{letter}",
                  f"panel_S10_{letter}")
