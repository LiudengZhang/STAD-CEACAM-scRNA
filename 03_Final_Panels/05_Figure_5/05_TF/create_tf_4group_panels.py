#!/usr/bin/env python3
"""Figure 5D and 5E: BACH1 (D) and NFKB1 (E) SCENIC regulon scores
(sc.tl.score_genes) in C3 MoMac cells, averaged per sample, across the four
treatment/response groups.

Post-R is compared with the other three groups pooled by an exact permutation
test on the difference in means; the panel prints the one-sided P (Post-R
lower) and the two-sided P is printed to the console for the legend. A
Kruskal-Wallis P among the other three groups is shown under their comb.

Per-sample scores are cached in data/sample_scores.csv; the h5ad is read only
when that table is missing or --recompute is passed.

Input:  MOMAC_H5AD (00_Config/paths.py)
Output: bach1_tf_4group.{svg,pdf,png} and nfkb1_tf_4group.{svg,pdf,png}
        in this directory
Run:    python create_tf_4group_panels.py [--recompute]
"""
import pandas as pd
import numpy as np
from scipy.stats import kruskal
from itertools import combinations as comb
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD
from cnsfig import Panel, style
from cnsfig.rich import rich_title
from cnsfig.boxes import box_xlim, draw_boxes, bracket, ylim_above, assert_no_points, one_vs_others  # noqa: E402
from cnsfig.layout import pin_frame_mm
from cnsfig import cache
import warnings
warnings.filterwarnings('ignore')

PANEL = {'BACH1': Panel(5, 'D', __file__, small_pt=5.0),
         'NFKB1': Panel(5, 'E', __file__, small_pt=5.0)}

OUT = str(Path(__file__).resolve().parent)
H5AD = str(MOMAC_H5AD)

# SCENIC regulon gene lists
REGULONS = {
    'BACH1': ['ACSL1', 'ARHGAP26', 'ASAP1', 'AZIN1', 'BACH1', 'BTG3', 'CDC42EP3',
              'CREM', 'CSGALNACT2', 'DSE', 'ELOVL5', 'FAM102B', 'FAM210A', 'FNDC3A',
              'FNDC3B', 'GPAT4', 'GPCPD1', 'HIVEP2', 'ITGB1', 'IVNS1ABP', 'JAK1',
              'JARID2', 'KAT6A', 'KCNA3', 'KDM3A', 'KMT2C', 'KMT2E', 'NAB1', 'NAMPT',
              'NFKB1', 'NR3C1', 'PAG1', 'PCNX1', 'PHF21A', 'PIM1', 'PPP3CA', 'RAP1B',
              'RASA2', 'SEC24A', 'SLC44A1', 'SPEN', 'TP53BP2', 'TRIP12', 'USP12',
              'ZFYVE16', 'ZNF395'],
    'NFKB1': ['ABCA1', 'ACSL1', 'ACSL4', 'AFF4', 'AFTPH', 'AKT3', 'ANKRD12', 'ARAP2',
              'ASAP1', 'ATP1B3', 'ATXN1', 'AZIN1', 'B3GNT5', 'B4GALT5', 'BACH1',
              'BASP1', 'BAZ1A', 'BTG3', 'CCNI', 'CDC42EP3', 'CEP170', 'CSGALNACT2',
              'CTNNB1', 'CYLD', 'DNAJB6', 'DSE', 'DUSP16', 'ELOVL7', 'EML4', 'EPB41L3',
              'F3', 'FAM102B', 'FAM107B', 'FNBP1', 'FNDC3A', 'FNDC3B', 'FRMD6',
              'GALNTL6', 'GPBP1', 'HIVEP1', 'HIVEP2', 'IVNS1ABP', 'JAK1', 'JARID2',
              'KCNJ2', 'KDM7A', 'KMT2E', 'KPNA4', 'LDLRAD4', 'LYN', 'MAPK6',
              'MIR155HG', 'MIR3945HG', 'N4BP2', 'NABP1', 'NAMPT', 'NFAT5', 'NFE2L2',
              'NFKB1', 'NR3C1', 'PDE4B', 'PELI1', 'PIM1', 'RAP1B', 'REL', 'SNX9',
              'SRSF12', 'STK26', 'SUSD6', 'TET2', 'TP53BP2', 'UAP1', 'USP12', 'WTAP',
              'ZFYVE16', 'ZSWIM6'],
}

# Response/phase colours shared by all figures.
from shared.figure_config import GROUP_COLORS, PRE_PAIR, POST_PAIR  # noqa: E402
COLORS = dict(GROUP_COLORS)

style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)


def exact_permutation_test(x, y, alternative='two-sided'):
    """Exact permutation test: enumerate all C(n, n_x) groupings."""
    all_vals = np.concatenate([x, y])
    n_total = len(all_vals)
    n_x = len(x)
    observed_stat = np.mean(x) - np.mean(y)
    count_extreme = 0
    n_perms = 0
    for idx in comb(range(n_total), n_x):
        perm_x = all_vals[list(idx)]
        perm_rest = all_vals[np.setdiff1d(range(n_total), idx)]
        perm_stat = np.mean(perm_x) - np.mean(perm_rest)
        if alternative == 'greater':
            if perm_stat >= observed_stat:
                count_extreme += 1
        elif alternative == 'less':
            if perm_stat <= observed_stat:
                count_extreme += 1
        else:
            if abs(perm_stat) >= abs(observed_stat):
                count_extreme += 1
        n_perms += 1
    return count_extreme / n_perms


# Per-sample regulon scores, cached as data/sample_scores.csv with columns
# tf, group, score (sample IDs are not written).
def compute_sample_scores():
    adata = sc.read_h5ad(H5AD)

    if adata.raw is not None:
        adata = adata.raw.to_adata()

    c3_mask = adata.obs['minor_cell_state'].astype(str).str.contains('C3')
    adata = adata[c3_mask].copy()

    for tf, genes in REGULONS.items():
        present = [g for g in genes if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=present, score_name=f'{tf}_score',
                          ctrl_size=len(present))

    pre_grouping = adata.obs['stomach_pre_grouping'].astype(str)
    post_grouping = adata.obs['stomach_post_grouping'].astype(str)

    adata.obs['group'] = 'Other'
    adata.obs.loc[pre_grouping == 'Responsed', 'group'] = 'Pre-R'
    adata.obs.loc[pre_grouping == 'No-response', 'group'] = 'Pre-NR'
    adata.obs.loc[post_grouping == 'Responsed', 'group'] = 'Post-R'
    adata.obs.loc[post_grouping == 'No-response', 'group'] = 'Post-NR'

    c3_stomach = adata[adata.obs['group'] != 'Other'].copy()

    out = []
    for tf_name in ['BACH1', 'NFKB1']:
        df = pd.DataFrame({
            'sample': c3_stomach.obs['sample'].values,
            'group': c3_stomach.obs['group'].values,
            'score': c3_stomach.obs[f'{tf_name}_score'].values,
        })
        sample_group = df.groupby('sample')['group'].first().reset_index()
        sample_mean = df.groupby('sample')['score'].mean().reset_index()
        merged = sample_group.merge(sample_mean, on='sample')
        merged.insert(0, 'tf', tf_name)
        out.append(merged[['tf', 'group', 'score']])
    return pd.concat(out, ignore_index=True)


SAMPLE_SCORES = cache.table(OUT, 'sample_scores', compute_sample_scores)

# Post-R last, set apart from the three groups it is compared with.
groups_order = ['Pre-R', 'Pre-NR', 'Post-NR', 'Post-R']
positions = [0, 1, 2, 3.5]

for tf_name in ['BACH1', 'NFKB1']:
    sample_scores = SAMPLE_SCORES[SAMPLE_SCORES['tf'] == tf_name]

    panel = PANEL[tf_name]
    fig, ax = style.subplots_mm(panel.w_mm, panel.h_mm)

    data_list = []
    colors_list = []
    for grp in groups_order:
        scores = sample_scores[sample_scores['group'] == grp]['score'].values
        data_list.append(scores if len(scores) > 0 else [np.nan])
        colors_list.append(COLORS[grp])

    bp = draw_boxes(ax, data_list, positions, colors_list, width=0.6)

    # Post-R vs the other three groups pooled (exact permutation test)
    post_r_vals = sample_scores[sample_scores['group'] == 'Post-R']['score'].values
    others_vals = np.concatenate([
        sample_scores[sample_scores['group'] == g]['score'].values
        for g in ['Pre-NR', 'Pre-R', 'Post-NR']
    ])
    perm_p_two_sided = exact_permutation_test(post_r_vals, others_vals,
                                             alternative='two-sided')
    perm_p = exact_permutation_test(post_r_vals, others_vals,
                                    alternative='less')

    other_data = [sample_scores[sample_scores['group'] == g]['score'].values
                  for g in ['Pre-NR', 'Pre-R', 'Post-NR']]
    if all(len(d) >= 2 for d in other_data):
        kw_stat, kw_p = kruskal(*other_data)
    else:
        kw_p = 1.0

    y_all = np.concatenate([d for d in data_list if not np.all(np.isnan(d))])
    y_max = np.nanmax(y_all)

    y_rng = y_max
    # Label shifted down to clear the panel letter.
    ax.set_ylabel('Regulon activity score', y=0.40)
    rich_title(ax, f'*{tf_name}* regulon')
    ax.set_xticks(positions)
    ax.set_xticklabels(groups_order, rotation=45, ha='right')
    ax.set_xlim(*box_xlim(positions, 0.6))
    assert_no_points(ax, bp)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Comb over the three other groups, a bracket "vs. others, P = x" from
    # Post-R, and the Kruskal-Wallis P among the three.
    p_text, _ = one_vs_others(fig, ax, positions[3], positions[:3], y_max,
                              y_rng, perm_p, kw_p, comb_lift=0.36)
    ax.set_ylim(top=y_max * 1.50)
    ylim_above(ax, p_text)
    panel.fit(fig)
    # Frame top 4.3 mm below the slot top, level with panels F and G.
    pin_frame_mm(fig, ax, top_mm=4.3,
                 bottom_mm=ax.get_position().y0 * panel.h_mm)
    ylim_above(ax, p_text)

    panel.finish(fig, f'{tf_name.lower()}_tf_4group')
    print(f"  panel {panel.letter}: permutation P = {perm_p:.4g} "
          f"one-tailed, {perm_p_two_sided:.4g} two-sided; KW P = {kw_p:.4f}")
