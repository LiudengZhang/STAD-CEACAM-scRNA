#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
"""
Figure 5 panels D and E - BACH1 and NFKB1 regulon activity across the four
treatment/response groups, one point per sample.

ONE SCRIPT, TWO PRINTED PANELS. `bach1_tf_4group` is printed panel D and
`nfkb1_tf_4group` is printed panel E (PROVENANCE.csv, confirmed against the
published page: panel D is titled BACH1 regulon and panel E NFKB1 regulon).
Each is drawn at its own printed rectangle, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py - the two are not
the same size and are stacked, not side by side.

The two SCENIC regulon gene lists, the C3 filter, `sc.tl.score_genes`, the
per-sample aggregation, the exact permutation test, the Kruskal-Wallis
homogeneity test, the bracket heights and both P-value strings are unchanged,
and `np.random.seed(42 + i)` stays exactly where it was. Both P values are
printed, each on its own line: the Kruskal-Wallis homogeneity P over the
three-group bracket and the permutation P over the four-group bracket.

The four group names are set on two lines - "Post" over "NR" - instead of on
one line at 45 degrees. The same four groups, in the same order, at the same
tick positions; only the line break and the rotation change, and the break is
declared in 00_Config/shared/labels.py.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the y axis label, the x tick labels, both sets of tick
    labels and the grey Kruskal-Wallis annotation - at 5 * SCALE, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    Box, whisker, cap, median and bracket line widths, the flier marker size
    and edge width, and the jitter marker area and edge width are scaled by
    those. Tick widths and lengths and spine widths are not: those are style,
    and cnsplots sets them.

Judgement calls, stated plainly:
  - The earlier drawing set three annotation sizes. cnsplots offers two, so the
    ordering is kept rather than the numbers: the grey Kruskal-Wallis note
    takes `tick_pt`, the permutation P takes `body_pt`, and the title takes the
    axes-title size.
  - The title's explicit `fontweight='normal'` is dropped so cnsplots' bold
    axis title applies.
"""
import pandas as pd
import numpy as np
from scipy.stats import kruskal
from itertools import combinations as comb
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD
import panel_style_cns as style
import slots
import warnings
warnings.filterwarnings('ignore')

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier smallest body type, before * SCALE

# Which printed panel each regulon is, and therefore which slot it is drawn in.
PANEL_LETTER = {'BACH1': 'D', 'NFKB1': 'E'}

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

COLORS = {
    'Pre-R':   '#bde0fe',
    'Post-R':  '#ffcfd2',
    'Pre-NR':  '#a2d2ff',
    'Post-NR': '#f1c0e8',
}

_family = style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
MARK = style.tick_pt() / (SMALL_PT * SCALE)
AREA = MARK ** 2
print(f"  type set in {_family}; body {style.body_pt():g} pt, "
      f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")


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


# ── Load MoMac and filter to C3_Mac ──
print("Loading MoMac h5ad...")
adata = sc.read_h5ad(H5AD)
print(f"Full MoMac: {adata.n_obs} cells, X genes: {adata.n_vars}")

if adata.raw is not None:
    print(f"Using .raw layer: {adata.raw.n_vars} genes")
    adata = adata.raw.to_adata()

c3_mask = adata.obs['minor_cell_state'].astype(str).str.contains('C3')
adata = adata[c3_mask].copy()
print(f"C3_Mac cells: {adata.n_obs}")

# ── Score regulon activity ──
print("\nScoring regulon activity with sc.tl.score_genes()...")
for tf, genes in REGULONS.items():
    present = [g for g in genes if g in adata.var_names]
    print(f"  {tf}: {len(present)}/{len(genes)} genes found")
    sc.tl.score_genes(adata, gene_list=present, score_name=f'{tf}_score',
                      ctrl_size=len(present))

# ── Assign 4 groups ──
pre_grouping = adata.obs['stomach_pre_grouping'].astype(str)
post_grouping = adata.obs['stomach_post_grouping'].astype(str)

adata.obs['group'] = 'Other'
adata.obs.loc[pre_grouping == 'Responsed', 'group'] = 'Pre-R'
adata.obs.loc[pre_grouping == 'No-response', 'group'] = 'Pre-NR'
adata.obs.loc[post_grouping == 'Responsed', 'group'] = 'Post-R'
adata.obs.loc[post_grouping == 'No-response', 'group'] = 'Post-NR'

c3_stomach = adata[adata.obs['group'] != 'Other'].copy()
print(f"\nC3_Mac stomach R/NR: {c3_stomach.n_obs}")
for grp in ['Pre-R', 'Post-R', 'Pre-NR', 'Post-NR']:
    n_cells = (c3_stomach.obs['group'] == grp).sum()
    n_samples = c3_stomach.obs.loc[c3_stomach.obs['group'] == grp, 'sample'].nunique()
    print(f"  {grp}: {n_cells} cells, {n_samples} samples")

# ── Create individual panels ──
groups_order = ['Pre-NR', 'Pre-R', 'Post-NR', 'Post-R']
positions = [0, 1, 2, 3.5]

for tf_name in ['BACH1', 'NFKB1']:
    score_col = f'{tf_name}_score'

    # Per-sample aggregation
    df = pd.DataFrame({
        'sample': c3_stomach.obs['sample'].values,
        'group': c3_stomach.obs['group'].values,
        'score': c3_stomach.obs[score_col].values,
    })
    sample_group = df.groupby('sample')['group'].first().reset_index()
    sample_mean = df.groupby('sample')['score'].mean().reset_index()
    sample_scores = sample_group.merge(sample_mean, on='sample')

    letter = PANEL_LETTER[tf_name]
    panel_w_mm, panel_h_mm = slots.size_mm(5, letter)
    letter_cell = slots.letter_cell_mm(5, letter)
    fig, ax = style.subplots_mm(panel_w_mm, panel_h_mm)

    data_list = []
    colors_list = []
    for grp in groups_order:
        scores = sample_scores[sample_scores['group'] == grp]['score'].values
        data_list.append(scores if len(scores) > 0 else [np.nan])
        colors_list.append(COLORS[grp])

    bp = ax.boxplot(data_list, positions=positions, widths=0.6, patch_artist=True,
                    boxprops=dict(linewidth=1.0 * SCALE * MARK),
                    whiskerprops=dict(color='black', linewidth=1.0 * SCALE * MARK),
                    capprops=dict(color='black', linewidth=1.0 * SCALE * MARK),
                    flierprops=dict(marker='o', markerfacecolor='white',
                                    markersize=4 * SCALE * MARK,
                                    markeredgecolor='black',
                                    markeredgewidth=0.5 * SCALE * MARK))
    for patch, color in zip(bp['boxes'], colors_list):
        patch.set_facecolor(color)
    for median in bp['medians']:
        median.set(color='black', linewidth=1.5 * SCALE * MARK)

    # Jitter points
    for i, (grp, pos) in enumerate(zip(groups_order, positions)):
        grp_data = sample_scores[sample_scores['group'] == grp]['score'].values
        if len(grp_data) > 0:
            np.random.seed(42 + i)
            jitter = np.random.uniform(-0.12, 0.12, size=len(grp_data))
            ax.scatter(np.full(len(grp_data), pos) + jitter, grp_data,
                       c='black', s=20 * SCALE * AREA, zorder=3,
                       edgecolors='white',
                       linewidths=0.3 * SCALE * MARK, alpha=0.85)

    # ── Stats: Post-R vs Others (exact permutation test) ──
    post_r_vals = sample_scores[sample_scores['group'] == 'Post-R']['score'].values
    others_vals = np.concatenate([
        sample_scores[sample_scores['group'] == g]['score'].values
        for g in ['Pre-NR', 'Pre-R', 'Post-NR']
    ])
    perm_p = exact_permutation_test(post_r_vals, others_vals, alternative='two-sided')
    print(f"  {tf_name}: Exact permutation p = {perm_p:.4f} (Post-R < Others)")

    # KW homogeneity among other 3
    other_data = [sample_scores[sample_scores['group'] == g]['score'].values
                  for g in ['Pre-NR', 'Pre-R', 'Post-NR']]
    if all(len(d) >= 2 for d in other_data):
        kw_stat, kw_p = kruskal(*other_data)
        print(f"  {tf_name}: KW among others p = {kw_p:.4f}")
    else:
        kw_p = 1.0

    # ── Brackets ──
    y_all = np.concatenate([d for d in data_list if not np.all(np.isnan(d))])
    y_max = np.nanmax(y_all)

    # Bracket 1 (lower): KW ns among Others (positions 0, 1, 2)
    ns_y = y_max * 1.08
    ax.plot([0, 0, 2, 2], [ns_y, ns_y * 1.03, ns_y * 1.03, ns_y],
            'k-', lw=0.8 * SCALE * MARK, alpha=0.6)
    # Exact P rather than a star: R1.3c asks for exact values, and the test
    # above is already two-sided.
    kw_str = f'P = {kw_p:.3f}' if kw_p >= 0.001 else 'P < 0.001'
    ax.text(1.0, ns_y * 1.04, kw_str, ha='center',
            fontsize=style.tick_pt(), color='#666666')

    # Bracket 2 (upper): Post-R vs Others — spans all 4 (positions 0 to 3.5)
    bracket_y = y_max * 1.22
    ax.plot([0, 0, 3.5, 3.5], [bracket_y, bracket_y * 1.03, bracket_y * 1.03, bracket_y],
            'k-', lw=0.8 * SCALE * MARK)
    p_str = f'P = {perm_p:.3f}' if perm_p >= 0.001 else 'P < 0.001'
    ax.text(1.75, bracket_y * 1.04, p_str, ha='center',
            fontsize=style.body_pt())

    ax.set_ylabel('Regulon activity score')
    ax.set_title(f'$\\it{{{tf_name}}}$ regulon')
    ax.set_xticks(positions)
    # Two lines rather than one line at 45 degrees: the same four names in
    # the same order at the same tick positions, in half the width.
    ax.set_xticklabels([g.replace('-', '\n') for g in groups_order])
    ax.set_ylim(top=y_max * 1.50)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=letter_cell)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {panel_w_mm} x {panel_h_mm} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, letter_cell)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    stem = f'{OUT}/{tf_name.lower()}_tf_4group'
    style.save_panel(fig, stem)
    print(f"  Saved: {stem}.[svg|pdf|png] at {panel_w_mm} x {panel_h_mm} mm "
          f"(printed panel {letter})")

print("\nDone — both TF panels created.")
