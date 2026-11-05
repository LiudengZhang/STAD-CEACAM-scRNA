#!/usr/bin/env python3
"""
Figures 2H (MP4) and 2I (MP5): per-sample stomach metaprogram scores by
treatment phase and response (Pre-R, Pre-NR, Post-R, Post-NR).

Statistics, per metaprogram:
  - exact permutation test of Pre-R against the other three groups pooled,
    over all groupings of the samples;
  - Kruskal-Wallis test among the other three groups.
The panel prints the exact one-tailed permutation P (Pre-R lower than the
others) on a bracket from Pre-R to a comb over the other three, and the
Kruskal-Wallis P under the comb. The two-sided permutation P is in the legend.

Reads data/sample_scores.csv, data/tests.csv and data/tests_directional.csv
beside this script; with --recompute, computes them from the NMF metaprogram
outputs (NMF_INTERMEDIATE, NMF_PER_SAMPLE) and the epithelial h5ad.
Writes mp4_horizontal_boxplot.* and mp5_horizontal_boxplot.* beside this
script, one figure per printed panel.

Run: python create_mp45_horizontal_boxplot.py [--recompute]
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import kruskal
from itertools import combinations as comb
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
from collections import Counter

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import NMF_INTERMEDIATE, NMF_PER_SAMPLE, EPITHELIAL_H5AD  # noqa: E402
from shared.expression import assert_x_stands_in_for_raw             # noqa: E402
from shared.sample_ids import sample_id_map, to_study_ids            # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig.boxes import box_xlim, draw_boxes, bracket, ylim_above, assert_no_points, one_vs_others  # noqa: E402
from cnsfig import cache                                  # noqa: E402
import slots                                              # noqa: E402

BASE_DIR = Path(__file__).parent
INTERMEDIATE = NMF_INTERMEDIATE
STOMACH_NMF_DIR = NMF_PER_SAMPLE
H5AD_PATH = EPITHELIAL_H5AD

# Shared four-group palette.
from shared.figure_config import GROUP_COLORS, PRE_PAIR, POST_PAIR  # noqa: E402
COLORS = dict(GROUP_COLORS)

SCALE = 4                           # sets MARK relative to the tick size
SMALL_PT = 5.0

# metaprogram -> (printed letter, output stem), in printed order
PRINTED = {
    'S-MP4': ('H', 'mp4_horizontal_boxplot'),
    'S-MP5': ('I', 'mp5_horizontal_boxplot'),
}

# Panel titles: the stomach metaprogram columns 'S-MP4'/'S-MP5' print as
# 'MP4'/'MP5', as in Figure 2G.
TITLE_OF = {'S-MP4': 'MP4', 'S-MP5': 'MP5'}


def get_program_genes(program_id, nmf_dir):
    """The top 50 genes of one per-sample NMF program."""
    parts = program_id.rsplit('_', 2)
    sample, k, p_num = parts[0], int(parts[1][1:]), int(parts[2][1:])
    genes_file = nmf_dir / sample / f'{sample}_nmf_k{k}_genes.csv'
    if genes_file.exists():
        df = pd.read_csv(genes_file)
        return list(df[df.columns[p_num - 1]].dropna().tolist()[:50])
    return []


def get_mp_consensus_genes(assignments, mp_id, nmf_dir, top_n=50):
    """A metaprogram's consensus genes, rank-weighted over its programs.

    Also used by 02_H/create_mp45_gene_dotplot.py (Figure 2G).
    """
    programs = assignments[assignments['metaprogram_id'] == mp_id]['program_id'].tolist()
    gene_counts = Counter()
    for pid in programs:
        for i, g in enumerate(get_program_genes(pid, nmf_dir)):
            gene_counts[g] += (50 - i)
    return [g for g, c in gene_counts.most_common(top_n)], gene_counts


def calculate_mp_score(adata, genes):
    raw = adata.raw if adata.raw is not None else adata
    valid = [g for g in genes if g in raw.var_names]
    if not valid:
        return np.zeros(adata.n_obs)
    expr = raw[:, valid].X
    if hasattr(expr, 'toarray'):
        expr = expr.toarray()
    return np.nanmean(expr, axis=1)


def exact_permutation_test(x, y, alternative='two-sided'):
    """
    Exact permutation test on mean(x) - mean(y), over all C(n, n_x)
    groupings. alternative='less' tests mean(x) < mean(y), 'greater' the
    reverse, 'two-sided' compares absolute differences.
    """
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
        if alternative == 'less':
            if perm_stat <= observed_stat:
                count_extreme += 1
        elif alternative == 'greater':
            if perm_stat >= observed_stat:
                count_extreme += 1
        else:
            if abs(perm_stat) >= abs(observed_stat):
                count_extreme += 1
        n_perms += 1
    return count_extreme / n_perms


def assign_group(row):
    phase = row['Treatment phase']
    if phase == 'Pre':
        resp = row['stomach_pre_grouping']
        if resp == 'Responsed': return 'Pre-R'
        elif resp == 'No-response': return 'Pre-NR'
    elif phase == 'Post':
        resp = row['stomach_post_grouping']
        if resp == 'Responsed': return 'Post-R'
        elif resp == 'No-response': return 'Post-NR'
    return None


# Group order used across the paper; Pre-R first, so the other three sit
# together under the comb.
GROUPS_ORDER = ['Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']
MP_COLS = ['S-MP4', 'S-MP5']


def compute_sample_scores():
    """Per-sample mean MP4/MP5 score and treatment group, from the h5ad."""
    stomach_assign = pd.read_csv(INTERMEDIATE / 'panel_A2_stomach_all_mp_assignments.csv')

    mp_genes = {}
    for mp in ['MP4', 'MP5']:
        genes, _ = get_mp_consensus_genes(stomach_assign, mp, STOMACH_NMF_DIR)
        mp_genes[f'S-{mp}'] = genes

    adata = sc.read_h5ad(H5AD_PATH)
    # Raises if this file has no .raw and its .X is not known to equal it.
    assert_x_stands_in_for_raw(adata, H5AD_PATH)
    ids = sample_id_map(adata.obs)
    adata_stomach = adata[adata.obs['Sample site'] == 'Stomach'].copy()

    for mp_name, genes in mp_genes.items():
        adata_stomach.obs[mp_name] = calculate_mp_score(adata_stomach, genes)

    mp_cols = ['S-MP4', 'S-MP5']
    sample_data = adata_stomach.obs.groupby('sample', observed=True).agg({
        **{mp: 'mean' for mp in mp_cols},
        'Patient ID': 'first',
        'Treatment phase': 'first',
        'stomach_pre_grouping': 'first',
        'stomach_post_grouping': 'first',
    }).reset_index()

    sample_data['group'] = sample_data.apply(assign_group, axis=1)
    sample_data = sample_data.dropna(subset=['group'])
    sample_data['sample'] = to_study_ids(sample_data['sample'], ids)  # Supplementary Table 1 IDs
    return sample_data[['sample', 'group'] + MP_COLS].reset_index(drop=True)


def compute_directional_tests(sample_data):
    """One-tailed permutation P, mean(Pre-R) < mean(others), the value the
    panel prints. Kept in its own table beside data/tests.csv.
    """
    rows = []
    for mp_name in MP_COLS:
        pre_r_vals = sample_data[sample_data['group'] == 'Pre-R'][mp_name].values
        others_vals = np.concatenate([
            sample_data[sample_data['group'] == g][mp_name].values
            for g in ['Post-R', 'Pre-NR', 'Post-NR']
        ])
        rows.append({'metaprogram': mp_name,
                     'perm_p_one_tailed': exact_permutation_test(
                         pre_r_vals, others_vals, alternative='less')})
    return pd.DataFrame(rows)


def compute_tests(sample_data):
    """Exact permutation P (Pre-R vs the other three pooled) and the
    Kruskal-Wallis P among the other three, per metaprogram."""
    rows = []
    for mp_name in MP_COLS:
        pre_r_vals = sample_data[sample_data['group'] == 'Pre-R'][mp_name].values
        others_vals = np.concatenate([
            sample_data[sample_data['group'] == g][mp_name].values
            for g in ['Post-R', 'Pre-NR', 'Post-NR']
        ])
        perm_p = exact_permutation_test(pre_r_vals, others_vals, alternative='two-sided')
        other_data = [sample_data[sample_data['group'] == g][mp_name].values
                      for g in ['Post-R', 'Pre-NR', 'Post-NR']]
        kw_stat, kw_p = kruskal(*other_data)
        rows.append({'metaprogram': mp_name, 'perm_p_two_sided': perm_p,
                     'kw_p_others': kw_p, 'kw_stat': kw_stat})
    return pd.DataFrame(rows)


def create_panel_F():
    sample_data = cache.table(BASE_DIR, 'sample_scores', compute_sample_scores)
    tests = cache.table(BASE_DIR, 'tests',
                        lambda: compute_tests(sample_data)).set_index('metaprogram')
    directional = cache.table(
        BASE_DIR, 'tests_directional',
        lambda: compute_directional_tests(sample_data)).set_index('metaprogram')

    groups_order = GROUPS_ORDER
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)

    for mp_name, (letter, stem) in PRINTED.items():
        panel_w, panel_h = slots.size_mm(2, letter, sub=1)
        letter_cell = slots.letter_cell_mm(2, letter)
        fig, ax = style.subplots_mm(panel_w, panel_h)

        data_list = []
        colors_list = []
        for grp_name in groups_order:
            scores = sample_data[sample_data['group'] == grp_name][mp_name].values
            data_list.append(scores if len(scores) > 0 else [np.nan])
            colors_list.append(COLORS[grp_name])

        bp = draw_boxes(ax, data_list, range(4), colors_list, width=0.6)

        # The panel prints the one-tailed P; the two-sided P goes to the
        # legend and Supplementary Table 6.
        perm_p_two_sided = float(tests.loc[mp_name, 'perm_p_two_sided'])
        perm_p = float(directional.loc[mp_name, 'perm_p_one_tailed'])
        kw_p = float(tests.loc[mp_name, 'kw_p_others'])

        y_all = np.concatenate(data_list)
        y_max = np.max(y_all)
        y_min = np.min(y_all)
        y_rng = y_max - y_min
        # Bracket from Pre-R to a comb over the other three ("vs. others,
        # P = x"), with the Kruskal-Wallis P under the comb. comb_lift keeps
        # the KW label clear of the highest outlier.
        p_text, _ = one_vs_others(fig, ax, 0, [1, 2, 3], y_max, y_rng,
                                  perm_p, kw_p, comb_lift=0.42)

        ax.set_ylabel('MP Score')
        ax.set_title(TITLE_OF[mp_name], color='black')
        ax.set_xticks(range(4))
        ax.set_xticklabels(groups_order, rotation=45, ha='right')
        ax.set_xlim(*box_xlim(range(4), 0.6))   # boxes off the spines
        ax.set_ylim(top=y_max + 0.75 * y_rng)
        ylim_above(ax, p_text)
        assert_no_points(ax, bp)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        style.fit_margins(fig, pad_mm=0.6, cell_mm=letter_cell)
        over = style.overflow_mm(fig)
        if max(over) > 0:
            raise RuntimeError(
                f"panel {letter}: ink outside the {panel_w} x {panel_h} mm "
                f"canvas (l,r,b,t mm): {over}")
        intruders = style.letter_clear(fig, letter_cell)
        if intruders:
            raise RuntimeError(f"panel {letter}: ink under the panel letter "
                               f"cell: {intruders}")

        style.save_panel(fig, BASE_DIR / stem, close=False)


if __name__ == '__main__':
    create_panel_F()
