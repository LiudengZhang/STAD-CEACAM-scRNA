#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# The one-tailed version this replaces is the one used in the preprint,
# https://www.biorxiv.org/content/10.64898/2026.03.05.708917
"""
Figure 2 panels H and I - stomach metaprogram scores MP4 and MP5 by treatment
phase and response.

- Pre-R against the other three groups, as a planned contrast
- Exact permutation test over all groupings of the samples
- Kruskal-Wallis homogeneity test among the other three groups
- Sample-level analysis

ONE SCRIPT, TWO DRAWINGS, TWO PRINTED LETTERS
    The page letters MP4 as panel H and MP5 as panel I, and gives each its own
    printed footprint with a 2.4 mm gutter between them. A slot-preserving
    assembler places one drawing in one slot, so each letter is drawn into its
    own figure at its own sub-box, read from 03_Final_Panels/slot_subrects.csv
    through 00_Config/slots.py. The two figures share this script's data, and
    nothing else: the scores, the groups, the tests, the brackets and every
    string are computed exactly as they were when the two axes sat side by side
    in one figure.

    Figures are built in printed order - H first, then I - so that the two
    drawings line up with the two axes of the earlier single figure.

  printed panels Figure 2 H (MP4) and Figure 2 I (MP5)
                 (PROVENANCE.csv; NOT inferred from "02_G")

The panels are drawn at the millimetre rectangles they print in, so the type
size set here is the type size printed. Margins are measured from the rendered
ink rather than typed, and each panel-letter corner is left clear for the
assembler.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the "ns" label - at 5 * SCALE. MARK carries the
    non-type point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    Every non-type length - the box, whisker, cap, median and the two bracket
    line widths, the flier marker and its edge - takes it. Tick widths and
    lengths and spine widths do not: those are style, and cnsplots sets them.

Every value read, every filter, every test and every string is the earlier
drawing's. The drawing code is the same code.
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
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402
import slots                                              # noqa: E402

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

SCALE = 4                           # the earlier canvas multiplier
SMALL_PT = 5.0                      # the earlier smallest body type

# metaprogram -> (printed letter, output stem), in printed order
PRINTED = {
    'S-MP4': ('H', 'mp4_horizontal_boxplot'),
    'S-MP5': ('I', 'mp5_horizontal_boxplot'),
}


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

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

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
        print(f"  Saved: {stem}.[svg|pdf|png] at {panel_w} x {panel_h} mm")


if __name__ == '__main__':
    create_panel_F()
