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
homogeneity test and the bracket heights are unchanged. The permutation P
prints over the upper bracket (a star below 0.05); the Kruskal-Wallis
bracket prints "ns" (2026-09-16, see below).

ONE BOX, NO POINTS, "ns" IN BLACK (2026-09-16, the author's fifth reading):
the boxes are cnsfig.boxes.draw_boxes (the 2H/I box; 0.79 mm open fliers),
the jittered sample points (`np.random.seed(42 + i)`) are no longer drawn -
the author's ruling, a declared departure from the published panels - and
the Kruskal-Wallis bracket prints "ns" in black instead of a grey
"P = 0.92" / "P = 0.90", as 2H/I does. Both brackets' vertices are the old
ones; the labels' ink sits 0.4 mm above their lines (cnsfig.boxes.bracket).

THE UPPER BRACKET, since the evening of 2026-09-15 (the author's fourth
reading): it runs from Post-R (x = 3.5) to the centre of the three-group
bracket (x = 1), not across all four groups as the submitted panel drew it,
so the page reads "one group against three" exactly as Figure 2 H/I does
("if 5D and 5E are also one against three, the bracket ends at the midpoint
of the three"). The test, its P value and the bracket height are unchanged.

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
from cnsfig.rich import rich_title
from cnsfig.boxes import box_xlim, draw_boxes, bracket, ylim_above, assert_no_points  # noqa: E402
from cnsfig.layout import pin_frame_mm
from cnsfig import cache
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


# ── Per-sample regulon scores, from data/sample_scores.csv (cnsfig.cache,
#    2026-09-15): the MoMac h5ad is read and scored only when the table is
#    absent or --recompute is passed. Columns: tf, group, score - one row per
#    stomach R/NR sample and regulon; the specimen number is not written. ──
def compute_sample_scores():
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

# ── Create individual panels ──
groups_order = ['Pre-NR', 'Pre-R', 'Post-NR', 'Post-R']
positions = [0, 1, 2, 3.5]

for tf_name in ['BACH1', 'NFKB1']:
    sample_scores = SAMPLE_SCORES[SAMPLE_SCORES['tf'] == tf_name]

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

    # ONE BOX, NO POINTS (2026-09-16, the author's fifth reading: "box
    # plots should look alike throughout; don't show every point"). The
    # 0.82 mm published flier is the family's 0.79 mm now; the jittered
    # sample points (seeded 42 + i) are gone - a declared departure from the
    # published panels, recorded in PROVENANCE.
    bp = draw_boxes(ax, data_list, positions, colors_list, width=0.6)

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

    # Bracket 1 (lower): KW among the other three (positions 0, 1, 2).
    # "ns" IN BLACK (2026-09-16, the fifth reading: "5D/E needn't print the
    # real P, ns is enough; and not grey") - kind="omnibus" prints ns when
    # P >= 0.05, as 2H/I does. The vertices are the old ones: the line at
    # 1.08 x y_max, its arms 3% of that above (lift/arm written against
    # y_rng = y_max), so only the label and its colour change.
    _, ns_text, _ = bracket(fig, ax, 0, 2, y_max, y_max, kw_p, kind="omnibus",
                            lift=0.08, arm=0.03 * 1.08)

    # Bracket 2 (upper): Post-R vs Others - from Post-R (3.5) to the centre of
    # the three-group bracket (1), the author's ruling of the evening of
    # 2026-09-15; it spanned all four positions (0 to 3.5) until then. The
    # star's ink 0.4 mm over the line.
    # THE LEFT ARM STANDS CLEAR OF "ns" (2026-09-16, verifier defect D1): at
    #   the old height (line at 1.22 x y_max) the arm's foot, at x = 1, came
    #   down onto the "ns" centred at the same x - 0.00 mm between them,
    #   where 2H/I has 3.5 mm between its two levels. The lift is raised
    #   until the foot clears the label's box by ARM_CLEAR_MM, measured with
    #   the renderer after ylim_above has set the limits (raising the limit
    #   compresses the axes, so it is measured again after each step); the
    #   span of the bracket and both ink gaps are unchanged.
    #   sweep_pages.check_bracket_clearance holds every page to it.
    # y=0.40: centred, the label's top end ran into the 10 pt letter cell once
    # the frame was pinned to the D/F/G line (2026-09-15).
    ax.set_ylabel('Regulon activity score', y=0.40)
    # The symbol italic, digit included (cnsfig.rich; 2026-09-14 evening).
    rich_title(ax, f'*{tf_name}* regulon')
    ax.set_xticks(positions)
    # One line at 45 degrees, as the published page sets them (2026-09-14).
    ax.set_xticklabels(groups_order, rotation=45, ha='right')
    ax.set_xlim(*box_xlim(positions, 0.6))   # boxes off the spines
    assert_no_points(ax, bp)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ARM_CLEAR_MM = 0.6
    lift = 0.22
    for _ in range(12):
        lines_up, p_text, _ = bracket(fig, ax, 1, 3.5, y_max, y_max, perm_p,
                                      kind="pair", lift=lift, arm=0.03 * 1.22)
        ax.set_ylim(top=y_max * 1.50)
        ylim_above(ax, p_text)
        style.fit_margins(fig, pad_mm=0.6, cell_mm=letter_cell)
        # ONE FRAME LINE FOR D, F AND G (sweep_pages.ROW_ALIGN): the frame
        # top at 4.3 mm below the slot's top, as F and G pin theirs. The
        # 9 pt star of 2026-09-15 had pushed this one 0.6 mm down through
        # fit_margins.
        pin_frame_mm(fig, ax, top_mm=4.3,
                     bottom_mm=ax.get_position().y0 * panel_h_mm)
        fig.canvas.draw()
        per_mm = fig.dpi / 25.4
        foot_px = ax.transData.transform((1.0, y_max + lift * y_max))[1]
        ns_top_px = ns_text.get_window_extent(fig.canvas.get_renderer()).y1
        clear_mm = (foot_px - ns_top_px) / per_mm
        if clear_mm >= ARM_CLEAR_MM:
            break
        for ln in lines_up:
            ln.remove()
        p_text.remove()
        # the deficit in mm, back into data units at the current scale
        lo, hi = ax.get_ylim()
        data_per_mm = (hi - lo) / (ax.get_window_extent().height / per_mm)
        lift += (ARM_CLEAR_MM - clear_mm + 0.1) * data_per_mm / y_max
    else:
        raise RuntimeError("the upper bracket never cleared the ns label")
    print(f"  {tf_name}: upper bracket lift {lift:.3f}; arm foot "
          f"{clear_mm:.2f} mm above the ns label's box")

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
