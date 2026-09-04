#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# Original one-tailed version: Round_5/03_Final_Panels/05_Figure_5/05_TF/create_tf_4group_panels.py
"""
Figure 5 panels D and E, RESTYLED (Version B) - BACH1 and NFKB1 regulon
activity across the four treatment/response groups, one point per sample.

ONE SCRIPT, TWO PRINTED PANELS. `bach1_tf_4group` is printed panel D and
`nfkb1_tf_4group` is printed panel E (PROVENANCE.csv). Both are drawn by the
loop below, in that order.

Version A is
`03_Revised_Panels/Main_Figures/05_Figure_5/05_TF/create_tf_4group_panels.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panels print in, non-type point sizes are rescaled by MARK
(areas by AREA), margins are millimetres, and the save is `style.save_panel`.

The two SCENIC regulon gene lists, the C3 filter, `sc.tl.score_genes`, the
per-sample aggregation, the exact permutation test, the Kruskal-Wallis
homogeneity test, the bracket heights and both P-value strings are Version A's,
unchanged. `np.random.seed(42 + i)` stays exactly where Version A put it.

  printed panels Figure 5 D (BACH1) and Figure 5 E (NFKB1)
  printed rects  D 32.9 x 32.8 mm, E 71.6 x 33.5 mm   (panel_rects.csv)
  Version B box  56.0 x 58.0 mm, both

    Four group names set at 45 degrees, a rotated y axis label, an italic
    two-word title and two stacked significance brackets do not fit in 33 mm at
    7/8 pt. Both panels are drawn at the same size because they are the same
    plot of two different regulons and are printed side by side.

MARK
    Version A drew 6 x 5.5 cm at SCALE = 4 and set its smallest body type - the
    y axis label, the x tick labels, both sets of tick labels and the grey
    Kruskal-Wallis annotation - at `5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2                                     = 0.1225

    Box, whisker, cap, median and bracket line widths, the flier marker size
    and edge width, and the jitter marker area and edge width are scaled by
    those. Tick widths/lengths and spine widths are not: those are style, and
    cnsplots sets them.

Judgement calls, stated plainly:
  - Version A set three annotation sizes (5, 6 and 7 pt nominal). cnsplots
    offers two, so the ordering is kept rather than the numbers: the grey
    Kruskal-Wallis note takes `tick_pt` (7), the permutation P takes `body_pt`
    (8), and the title takes the axes-title size.
  - The title's explicit `fontweight='normal'` is dropped so cnsplots' bold
    axis title applies (PANEL_SPEC: follow the library).
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
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import MOMAC_H5AD
import panel_style_cns as style
import warnings
warnings.filterwarnings('ignore')

SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 5.0                  # Version A's smallest body type, before * SCALE

PRINTED_MM = {'BACH1': (32.9, 32.8), 'NFKB1': (71.6, 33.5)}  # panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 56.0, 58.0
MARGIN = dict(left=14.0, right=3.0, top=9.0, bottom=12.0)

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

_family = style.apply()
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

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

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
    ax.set_xticklabels(groups_order, rotation=45, ha='right')
    ax.set_ylim(top=y_max * 1.50)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    stem = f'{OUT}/{tf_name.lower()}_tf_4group'
    style.save_panel(fig, stem)
    print(f"  Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")

print("\nDone — both TF panels created.")
