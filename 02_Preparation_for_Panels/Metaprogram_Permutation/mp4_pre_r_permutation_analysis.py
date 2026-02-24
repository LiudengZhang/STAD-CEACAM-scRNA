#!/usr/bin/env python3
"""
MP4 Pre-R vs Others: Exact Permutation Test + Effect Size Analysis
Hypothesis: Pre-treatment responders have uniquely low S-MP4 scores
compared to all other groups (Post-R, Pre-NR, Post-NR), which are
homogeneous among themselves.

Statistical framework:
1. Exact permutation test (a priori planned contrast, Pre-R < Others)
2. Kruskal-Wallis among other 3 groups (homogeneity confirmation)
3. Rank-biserial effect size + bootstrap 95% CI
4. Visualization: 1-vs-3 boxplot
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import kruskal, mannwhitneyu
from collections import Counter
from itertools import combinations
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import json

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from paths import *

OUTPUT_DIR = Path(__file__).parent

# 4x scaling
SCALE = 4
DPI = 300
CM_TO_INCH = 1 / 2.54

# =========================================================================
# Helper functions (from existing panel F)
# =========================================================================
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
    return [g for g, c in gene_counts.most_common(top_n)]


def calculate_mp_score(adata, genes):
    raw = adata.raw if adata.raw is not None else adata
    valid = [g for g in genes if g in raw.var_names]
    if not valid:
        return np.zeros(adata.n_obs)
    expr = raw[:, valid].X
    if hasattr(expr, 'toarray'):
        expr = expr.toarray()
    return np.nanmean(expr, axis=1)


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


# =========================================================================
# Exact permutation test
# =========================================================================
def exact_permutation_test(x, y, alternative='less'):
    """
    Exact permutation test comparing group x vs group y.
    alternative='less': test if mean(x) < mean(y)
    Returns exact p-value by enumerating all C(n, n_x) permutations.
    """
    all_vals = np.concatenate([x, y])
    n_total = len(all_vals)
    n_x = len(x)
    observed_stat = np.mean(x) - np.mean(y)

    count_extreme = 0
    n_perms = 0
    for idx in combinations(range(n_total), n_x):
        perm_x = all_vals[list(idx)]
        perm_stat = np.mean(perm_x) - np.mean(all_vals[np.setdiff1d(range(n_total), idx)])
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

    return count_extreme / n_perms, n_perms, observed_stat


# =========================================================================
# Main
# =========================================================================
def main():
    results = {}

    # --- Load data ---
    print("Loading data...")
    stomach_assign = pd.read_csv(NMF_INTERMEDIATE / 'panel_A2_stomach_all_mp_assignments.csv')

    mp_genes = {}
    for mp in ['MP4', 'MP5']:
        genes = get_mp_consensus_genes(stomach_assign, mp, NMF_PER_SAMPLE)
        mp_genes[f'S-{mp}'] = genes
        print(f"  S-{mp}: {len(genes)} consensus genes")

    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    adata_stomach = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    print(f"  Stomach epithelial cells: {adata_stomach.n_obs}")

    for mp_name, genes in mp_genes.items():
        adata_stomach.obs[mp_name] = calculate_mp_score(adata_stomach, genes)

    # --- Sample-level aggregation ---
    sample_data = adata_stomach.obs.groupby('sample', observed=True).agg({
        'S-MP4': 'mean', 'S-MP5': 'mean',
        'Treatment phase': 'first',
        'stomach_pre_grouping': 'first',
        'stomach_post_grouping': 'first',
    }).reset_index()

    sample_data['group'] = sample_data.apply(assign_group, axis=1)
    sample_data = sample_data.dropna(subset=['group'])
    sample_data['S-MP4+5'] = (sample_data['S-MP4'] + sample_data['S-MP5']) / 2

    groups_order = ['Pre-R', 'Post-R', 'Pre-NR', 'Post-NR']
    print("\nSample counts:")
    for g in groups_order:
        n = len(sample_data[sample_data['group'] == g])
        print(f"  {g}: n={n}")

    # --- Run analysis for MP4, MP5, and combined ---
    for mp in ['S-MP4', 'S-MP5', 'S-MP4+5']:
        print(f"\n{'='*60}")
        print(f"  {mp}")
        print(f"{'='*60}")

        mp_results = {}
        group_data = {}
        for g in groups_order:
            vals = sample_data[sample_data['group'] == g][mp].values
            group_data[g] = vals
            mp_results[f'{g}_n'] = int(len(vals))
            mp_results[f'{g}_mean'] = float(np.mean(vals))
            mp_results[f'{g}_median'] = float(np.median(vals))
            mp_results[f'{g}_values'] = vals.tolist()
            print(f"  {g:8s}: n={len(vals)}, mean={np.mean(vals):.4f}, "
                  f"median={np.median(vals):.4f}")

        # 1. Kruskal-Wallis (4 groups)
        kw_stat, kw_p = kruskal(*[group_data[g] for g in groups_order])
        mp_results['kw_4group_H'] = float(kw_stat)
        mp_results['kw_4group_p'] = float(kw_p)
        print(f"\n  [1] Kruskal-Wallis (4 groups): H={kw_stat:.3f}, p={kw_p:.4f}")

        # 2. Exact permutation test: Pre-R < Others
        pre_r_vals = group_data['Pre-R']
        others_vals = np.concatenate([group_data[g] for g in ['Post-R', 'Pre-NR', 'Post-NR']])
        perm_p, n_perms, obs_diff = exact_permutation_test(
            pre_r_vals, others_vals, alternative='less')
        mp_results['permutation_p'] = float(perm_p)
        mp_results['permutation_n_perms'] = int(n_perms)
        mp_results['permutation_observed_diff'] = float(obs_diff)
        print(f"  [2] Exact permutation (Pre-R < Others): p={perm_p:.4f} "
              f"({n_perms} permutations, obs_diff={obs_diff:.4f})")

        # 3. KW among other 3 groups (homogeneity)
        other_groups = ['Post-R', 'Pre-NR', 'Post-NR']
        kw_other_stat, kw_other_p = kruskal(*[group_data[g] for g in other_groups])
        mp_results['kw_others_H'] = float(kw_other_stat)
        mp_results['kw_others_p'] = float(kw_other_p)
        print(f"  [3] KW among other 3 (homogeneity): H={kw_other_stat:.3f}, "
              f"p={kw_other_p:.4f}")

        # 4. Mann-Whitney (Pre-R < Others) for comparison
        mw_stat, mw_p = mannwhitneyu(pre_r_vals, others_vals, alternative='less')
        mp_results['mannwhitney_U'] = float(mw_stat)
        mp_results['mannwhitney_p'] = float(mw_p)
        print(f"  [4] Mann-Whitney U (Pre-R < Others): U={mw_stat}, p={mw_p:.4f}")

        # 5. Effect size: rank-biserial correlation
        U_two, _ = mannwhitneyu(pre_r_vals, others_vals, alternative='two-sided')
        n1, n2 = len(pre_r_vals), len(others_vals)
        r_rb = 1 - (2 * U_two) / (n1 * n2)
        mp_results['rank_biserial_r'] = float(r_rb)
        print(f"  [5] Rank-biserial r = {r_rb:.3f}")

        # 6. Bootstrap 95% CI for mean difference
        np.random.seed(42)
        n_boot = 10000
        boot_diffs = []
        for _ in range(n_boot):
            bx = np.random.choice(pre_r_vals, n1, replace=True)
            by = np.random.choice(others_vals, n2, replace=True)
            boot_diffs.append(np.mean(bx) - np.mean(by))
        ci_lo, ci_hi = np.percentile(boot_diffs, [2.5, 97.5])
        mp_results['bootstrap_mean_diff'] = float(np.mean(pre_r_vals) - np.mean(others_vals))
        mp_results['bootstrap_ci_lo'] = float(ci_lo)
        mp_results['bootstrap_ci_hi'] = float(ci_hi)
        mp_results['bootstrap_ci_excludes_zero'] = bool(ci_lo > 0 or ci_hi < 0)
        print(f"  [6] Bootstrap 95% CI: [{ci_lo:.4f}, {ci_hi:.4f}] "
              f"(excludes 0: {mp_results['bootstrap_ci_excludes_zero']})")

        results[mp] = mp_results

    # --- Save results as JSON ---
    results_file = OUTPUT_DIR / 'mp4_permutation_results.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved: {results_file}")

    # =====================================================================
    # Visualization: 1-vs-3 boxplots for MP4, MP5, and combined
    # =====================================================================
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'font.size': 6 * SCALE,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    COLORS = {
        'Pre-R':  '#bde0fe',
        'Post-R': '#ffcfd2',
        'Pre-NR': '#a2d2ff',
        'Post-NR': '#f1c0e8',
    }

    plot_configs = [
        ('S-MP4', 'S-MP4 (CEACAM-associated)', 'mp4_pre_r_vs_others'),
        ('S-MP5', 'S-MP5 (CEACAM-associated)', 'mp5_pre_r_vs_others'),
        ('S-MP4+5', 'S-MP4+5 Combined', 'mp45_combined_pre_r_vs_others'),
    ]

    for mp, title, filename in plot_configs:
        print(f"\nCreating {mp} boxplot...")
        fig, ax = plt.subplots(figsize=(4.5 * SCALE * CM_TO_INCH, 3.5 * SCALE * CM_TO_INCH))

        data_list = []
        colors_list = []
        for g in groups_order:
            vals = sample_data[sample_data['group'] == g][mp].values
            data_list.append(vals)
            colors_list.append(COLORS[g])

        bp = ax.boxplot(data_list, positions=[0, 1, 2, 3], widths=0.6, patch_artist=True,
                        boxprops=dict(linewidth=0.8),
                        whiskerprops=dict(color='black', linewidth=0.8),
                        capprops=dict(color='black', linewidth=0.8),
                        flierprops=dict(marker='o', markerfacecolor='white', markersize=4,
                                       markeredgecolor='black', markeredgewidth=0.5),
                        medianprops=dict(color='black', linewidth=1.2))

        for patch, color in zip(bp['boxes'], colors_list):
            patch.set_facecolor(color)

        # Bracket: Pre-R vs Others (positions 0 vs {1,2,3})
        y_all = np.concatenate(data_list)
        y_max = np.max(y_all)
        bracket_y = y_max * 1.08
        ax.plot([0, 0, 2, 2], [bracket_y, bracket_y * 1.03, bracket_y * 1.03, bracket_y],
                'k-', linewidth=0.8)

        perm_p = results[mp]['permutation_p']
        p_text = f'P = {perm_p:.3f}' if perm_p >= 0.001 else f'P < 0.001'
        ax.text(1.0, bracket_y * 1.05, p_text, ha='center', va='bottom', fontsize=5 * SCALE)

        # ns bracket among others
        ns_y = y_max * 1.22
        ax.plot([1, 1, 3, 3], [ns_y, ns_y * 1.02, ns_y * 1.02, ns_y],
                'k-', linewidth=0.5, alpha=0.6)
        ax.text(2.0, ns_y * 1.03, 'ns', ha='center', va='bottom',
                fontsize=5 * SCALE, color='#666666')

        ax.set_ylabel(f'{mp} Score', fontsize=6 * SCALE)
        ax.set_title(title, fontsize=6 * SCALE, fontweight='normal')
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['Pre-R', 'Post-R', 'Pre-NR', 'Post-NR'],
                           fontsize=5 * SCALE, rotation=45, ha='right')
        ax.set_ylim(0, y_max * 1.45)
        ax.tick_params(axis='both', labelsize=5 * SCALE, width=0.8, length=3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.8)
        ax.spines['bottom'].set_linewidth(0.8)

        plt.tight_layout()

        for ext in ['png', 'svg', 'pdf']:
            plt.savefig(OUTPUT_DIR / f'{filename}.{ext}',
                        dpi=DPI, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"  Saved: {OUTPUT_DIR / f'{filename}.png'}")

    # =====================================================================
    # Print summary table
    # =====================================================================
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'Metric':<40} {'S-MP4':>10} {'S-MP5':>10} {'Combined':>10}")
    print("-" * 70)
    for key, label in [
        ('permutation_p', 'Exact permutation p (Pre-R < Others)'),
        ('kw_4group_p', 'Kruskal-Wallis p (4 groups)'),
        ('kw_others_p', 'KW p (other 3 homogeneity)'),
        ('rank_biserial_r', 'Rank-biserial effect size'),
        ('bootstrap_ci_excludes_zero', 'Bootstrap 95% CI excludes 0'),
    ]:
        vals = []
        for mp_key in ['S-MP4', 'S-MP5', 'S-MP4+5']:
            v = results[mp_key][key]
            if isinstance(v, bool):
                vals.append('Yes' if v else 'No')
            elif isinstance(v, float):
                vals.append(f'{v:.4f}')
            else:
                vals.append(str(v))
        print(f"  {label:<38} {vals[0]:>10} {vals[1]:>10} {vals[2]:>10}")
    print("=" * 70)


if __name__ == '__main__':
    main()
