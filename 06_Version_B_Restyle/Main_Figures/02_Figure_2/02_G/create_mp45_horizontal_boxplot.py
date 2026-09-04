#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# Original one-tailed version: Round_5/03_Final_Panels/02_Figure_2/02_G/create_mp45_horizontal_boxplot.py
"""
Figure 2, printed panels H AND I, RESTYLED (Version B) - MP4 and MP5 boxplots
side by side.

- Pre-R vs Others (Post-R, Pre-NR, Post-NR) planned contrast
- Exact permutation test (Fisher's exact permutation)
- KW homogeneity test among the other 3 groups
- Sample-level analysis
- Two subplots side by side (MP4 left, MP5 right)

Version A is
`03_Final_Panels/02_Figure_2/02_G/create_mp45_horizontal_boxplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every test and every string is Version A's.
The drawing code is the same code.

  printed panels Figure 2 H AND Figure 2 I     (PROVENANCE.csv; NOT inferred
                 from "02_G"). This is ONE drawing with TWO side-by-side axes,
                 which the paper letters H (MP4) and I (MP5). It stays ONE
                 file and ONE figure - the assembler places it once and puts
                 two letters on it. Do not split it.
  printed rect   H 35.1 x 28.8 mm and I 33.4 x 28.9 mm (panel_rects.csv), so
                 about 70 x 29 mm for the drawing as a whole
  Version B box  100.0 x 50.0 mm

MARK
    Version A drew at SCALE = 4 (8.0 x 3.5 cm x 4 = 320 x 140 mm) and its
    smallest body type is the "ns" label at `5 * SCALE`. So

        SCALE = 4, SMALL_PT = 5
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA  = MARK ** 2 = 0.1225

    Every non-type length - box, whisker, cap, median and the two bracket line
    widths, the flier marker and its edge - is multiplied by MARK.

    Version A's spine linewidths and `tick_params(width=..., length=...)` are
    NOT carried over: those are axes furniture, cnsplots has its own settings
    for them, and following the library rather than rescaling the old numbers
    is the standard-methods rule.

    The drawing grew from about 70 x 29 mm to 100 x 50 mm: each axes carries
    four 45-degree group labels ("Post-NR" is the longest) plus a y label and a
    title at 7/8 pt, and both stacked brackets have to clear the boxes.
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

# Central config
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402

BASE_DIR = Path(__file__).parent
INTERMEDIATE = NMF_INTERMEDIATE                           # noqa: F405
STOMACH_NMF_DIR = NMF_PER_SAMPLE                          # noqa: F405
H5AD_PATH = EPITHELIAL_H5AD                               # noqa: F405

# Colors - standard 4-group palette
COLORS = {
    'Pre-R':   '#bde0fe',
    'Post-R':  '#ffcfd2',
    'Pre-NR':  '#a2d2ff',
    'Post-NR': '#f1c0e8',
}

PRINTED_MM = ((35.1, 28.8), (33.4, 28.9))   # published rects of H and I
PANEL_W_MM, PANEL_H_MM = 100.0, 50.0
MARGIN = dict(left=11.0, right=3.0, top=5.5, bottom=11.0, wspace=0.42)

SCALE = 4                           # Version A's canvas multiplier
SMALL_PT = 5.0                      # Version A's smallest body type


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


def exact_permutation_test(x, y, alternative='two-sided'):
    """
    Exact permutation test: enumerate all C(n, n_x) groupings.
    alternative='two-sided': test if mean(x) < mean(y).
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
    print("Creating printed panels H and I: MP4, MP5 Horizontal Boxplots...")

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

    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    # 1 row, 2 columns - side by side. ONE figure, TWO axes: printed H and I.
    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, 1, 2, sharey=False)

    focus_mps = ['S-MP4', 'S-MP5']

    for i, mp_name in enumerate(focus_mps):
        ax = axes[i]

        data_list = []
        colors_list = []
        for grp_name in groups_order:
            scores = sample_data[sample_data['group'] == grp_name][mp_name].values
            data_list.append(scores if len(scores) > 0 else [np.nan])
            colors_list.append(COLORS[grp_name])

        bp = ax.boxplot(data_list, positions=range(4), widths=0.6,
                        patch_artist=True,
                        boxprops=dict(linewidth=0.5 * MARK),
                        whiskerprops=dict(color='black', linewidth=0.5 * MARK),
                        capprops=dict(color='black', linewidth=0.5 * MARK),
                        flierprops=dict(marker='o', markerfacecolor='white',
                                        markersize=4 * MARK,
                                        markeredgecolor='black',
                                        markeredgewidth=0.5 * MARK))

        for patch, color in zip(bp['boxes'], colors_list):
            patch.set_facecolor(color)
        for median in bp['medians']:
            median.set(color='black', linewidth=0.8 * MARK)

        # Stats: Exact permutation test (Pre-R < Others)
        pre_r_vals = sample_data[sample_data['group'] == 'Pre-R'][mp_name].values
        others_vals = np.concatenate([
            sample_data[sample_data['group'] == g][mp_name].values
            for g in ['Post-R', 'Pre-NR', 'Post-NR']
        ])
        perm_p = exact_permutation_test(pre_r_vals, others_vals, alternative='two-sided')
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
                'k-', lw=0.5 * MARK)
        p_str = f'P = {perm_p:.3f}' if perm_p >= 0.001 else 'P < 0.001'
        ax.text(1.0, bracket_y*1.04, p_str, ha='center',
                fontsize=style.tick_pt())

        # ns bracket among others
        ns_y = y_max * 1.22
        ax.plot([1, 1, 3, 3], [ns_y, ns_y*1.02, ns_y*1.02, ns_y],
                'k-', lw=0.5 * MARK, alpha=0.6)
        ax.text(2.0, ns_y*1.03, 'ns', ha='center', fontsize=style.tick_pt(),
                color='#666666')

        ax.set_ylabel('MP Score')
        ax.set_title(mp_name, color='black')
        ax.set_xticks(range(4))
        ax.set_xticklabels(['Pre-R', 'Post-R', 'Pre-NR', 'Post-NR'],
                           rotation=45, ha='right')
        ax.set_ylim(top=y_max * 1.45)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    style.save_panel(fig, BASE_DIR / 'mp45_horizontal_boxplot')
    print(f"  Saved: {BASE_DIR / 'mp45_horizontal_boxplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    create_panel_F()
