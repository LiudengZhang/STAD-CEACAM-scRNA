#!/usr/bin/env python3
"""
Panel F: MP4 and MP5 Boxplots Side by Side (horizontal layout)
- Pre-R vs Others (Post-R, Pre-NR, Post-NR) planned contrast
- Exact permutation test (Fisher's exact permutation)
- KW homogeneity test among the other 3 groups
- Sample-level analysis
- Two subplots side by side (MP4 left, MP5 right)
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import kruskal
from itertools import combinations as comb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter

# Central config
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

BASE_DIR = Path(__file__).parent
INTERMEDIATE = NMF_INTERMEDIATE
STOMACH_NMF_DIR = NMF_PER_SAMPLE
H5AD_PATH = EPITHELIAL_H5AD

# Colors — standard 4-group palette
COLORS = {
    'Pre-R':   '#bde0fe',
    'Post-R':  '#ffcfd2',
    'Pre-NR':  '#a2d2ff',
    'Post-NR': '#f1c0e8',
}

# Nature Cancer 4× scaling
DPI = 300
SCALE = 4
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 8.0 * SCALE    # wide for side-by-side
PANEL_HEIGHT_CM = 3.5 * SCALE   # shorter row 2 height


def get_program_genes(program_id, nmf_dir):
    parts = program_id.rsplit('_', 2)
    sample, k, p_num = parts[0], int(parts[1][1:]), int(parts[2][1:])
    genes_file = nmf_dir / sample / f'{sample}_nmf_k{k}_genes.csv'
    if genes_file.exists():
        df = pd.read_csv(genes_file)
        return list(df[df.columns[p_num - 1]].dropna().tolist()[:50])
    return []


def get_mp_consensus_genes(assignments, mp_id, nmf_dir, top_n=50):
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


def exact_permutation_test(x, y, alternative='less'):
    """
    Exact permutation test: enumerate all C(n, n_x) groupings.
    alternative='less': test if mean(x) < mean(y).
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


def create_panel_F():
    print("Creating Panel F: MP4, MP5 Horizontal Boxplots...")

    stomach_assign = pd.read_csv(INTERMEDIATE / 'panel_A2_stomach_all_mp_assignments.csv')

    mp_genes = {}
    for mp in ['MP4', 'MP5']:
        genes, _ = get_mp_consensus_genes(stomach_assign, mp, STOMACH_NMF_DIR)
        mp_genes[f'S-{mp}'] = genes
        print(f"  S-{mp}: {len(genes)} genes")

    adata = sc.read_h5ad(H5AD_PATH)
    adata_stomach = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    print(f"  Stomach: {adata_stomach.n_obs} cells")

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

    groups_order = ['Pre-R', 'Post-R', 'Pre-NR', 'Post-NR']
    for g in groups_order:
        n = len(sample_data[sample_data['group'] == g])
        print(f"  {g}: n={n}")

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    # 1 row, 2 columns — side by side
    fig, axes = plt.subplots(1, 2, figsize=(PANEL_WIDTH_CM * CM_TO_INCH, PANEL_HEIGHT_CM * CM_TO_INCH), sharey=False)

    focus_mps = ['S-MP4', 'S-MP5']

    for i, mp_name in enumerate(focus_mps):
        ax = axes[i]

        data_list = []
        colors_list = []
        for grp_name in groups_order:
            scores = sample_data[sample_data['group'] == grp_name][mp_name].values
            data_list.append(scores if len(scores) > 0 else [np.nan])
            colors_list.append(COLORS[grp_name])

        bp = ax.boxplot(data_list, positions=range(4), widths=0.6, patch_artist=True,
                        boxprops=dict(linewidth=0.5),
                        whiskerprops=dict(color='black', linewidth=0.5),
                        capprops=dict(color='black', linewidth=0.5),
                        flierprops=dict(marker='o', markerfacecolor='white', markersize=4,
                                       markeredgecolor='black', markeredgewidth=0.5))

        for patch, color in zip(bp['boxes'], colors_list):
            patch.set_facecolor(color)
        for median in bp['medians']:
            median.set(color='black', linewidth=0.8)

        # Stats: Exact permutation test (Pre-R < Others)
        pre_r_vals = sample_data[sample_data['group'] == 'Pre-R'][mp_name].values
        others_vals = np.concatenate([
            sample_data[sample_data['group'] == g][mp_name].values
            for g in ['Post-R', 'Pre-NR', 'Post-NR']
        ])
        perm_p = exact_permutation_test(pre_r_vals, others_vals, alternative='less')
        print(f"  {mp_name}: Exact permutation p = {perm_p:.4f} (Pre-R < Others)")

        # KW homogeneity among other 3
        other_data = [sample_data[sample_data['group'] == g][mp_name].values
                      for g in ['Post-R', 'Pre-NR', 'Post-NR']]
        kw_stat, kw_p = kruskal(*other_data)
        print(f"  {mp_name}: KW among others p = {kw_p:.4f}")

        # Bracket: Pre-R vs Others (position 0 vs midpoint of 1,2,3)
        y_all = np.concatenate(data_list)
        y_max = np.max(y_all)
        bracket_y = y_max * 1.08
        ax.plot([0, 0, 2, 2], [bracket_y, bracket_y*1.03, bracket_y*1.03, bracket_y],
                'k-', lw=0.5)
        p_str = '***' if perm_p < 0.001 else '**' if perm_p < 0.01 else '*' if perm_p < 0.05 else 'ns'
        ax.text(1.0, bracket_y*1.04, p_str, ha='center', fontsize=6 * SCALE)

        # ns bracket among others
        ns_y = y_max * 1.22
        ax.plot([1, 1, 3, 3], [ns_y, ns_y*1.02, ns_y*1.02, ns_y],
                'k-', lw=0.5, alpha=0.6)
        ax.text(2.0, ns_y*1.03, 'ns', ha='center', fontsize=5 * SCALE, color='#666666')

        ax.set_ylabel('MP Score', fontsize=6 * SCALE)
        ax.set_title(mp_name, fontsize=7 * SCALE, color='black', fontweight='normal')
        ax.set_xticks(range(4))
        ax.set_xticklabels(['Pre-R', 'Post-R', 'Pre-NR', 'Post-NR'],
                           fontsize=6 * SCALE, rotation=45, ha='right')
        ax.set_ylim(top=y_max * 1.45)
        ax.tick_params(axis='both', labelsize=6 * SCALE, width=0.5, length=4)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['bottom'].set_linewidth(0.5)

    plt.tight_layout()

    output = BASE_DIR / 'mp45_horizontal_boxplot.png'
    plt.savefig(output, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.svg'), dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(BASE_DIR / 'mp45_horizontal_boxplot.pdf', format='pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output}")


if __name__ == '__main__':
    create_panel_F()
