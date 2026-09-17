#!/usr/bin/env python3
"""
Figure 3 panel E - CEACAM5 and CEACAM6 against the Tex ssGSEA score in TCGA-STAD, one point
per tumour, open or filled by vital status.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 3 E       (PROVENANCE.csv; the directory is 03_D, and
                                   the letter was looked up, not inferred)

  03_C/create_tcga_ceacam_scatter.py is a different file with the same
  name drawing printed panel B; the two are not interchangeable.

Every value read, every filter, every statistic and every string is the earlier
drawing's: the deconvolved epithelial expression, the ssGSEA scoring, the
per-subplot zero-expression filter, the Spearman test and the P-value
formatting are untouched.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the in-axes legend - at 4.5 * SCALE. MARK carries the
    non-type point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The marker area, the marker edge, the dashed regression line and the two
    legend keys are scaled by it. Tick widths and lengths and spine widths are
    not: those are style, and cnsplots sets them.

    The legend keys are also divided by the legend's marker scale, which
    cnsplots sets to a half and which the earlier drawing left at one, so that
    each key prints at MARK times the size it had rather than half of that -
    which would close the open ring of the unfilled key into a dot.

THE COHORT LEAVES THE TITLE, AND THE AXIS LABELS ARE RE-WRAPPED
    Both axes carry the same cohort name and its sample count. The cohort and
    the two counts are named in the caption and each axes keeps its own
    correlation and its own P value; on one line the pair sets about 24 mm
    against a plotting box of about 10 mm, so they take a line each.

    The x label carries the gene alone: at 7 pt the unit sets 11.0 mm and the
    two plotting boxes are about 9 mm wide, so the two axes' units meet between
    them. The y label is re-wrapped: it is rotated, so its length is vertical, and the full form sets
    24.5 mm against a 21.8 mm panel. Both are declared in RENAMES_FIGURE_3.

    The subscript of log2 is written out. Mathtext draws a subscript at 70% of
    its base size, so it would print at 4.9 pt on a 7 pt label, below the floor
    this figure set is set to.

MARKER AND RULE WIDTHS ARE MEASURED OFF THE PUBLISHED PAGE  (2026-09-11)
    Every one of these scatter panels carried a stray factor of SCALE on its
    non-type sizes - `s=30*SCALE*AREA`, `linewidth=0.8*SCALE*MARK` - on top of
    AREA and MARK, which already carry the 4x canvas across. The markers came
    out about twice as wide as the page prints them and ran together.

    Counted out of `00_GROUND_TRUTH/figures/Figure 3.pdf` geometry:

        panel A   0.520 mm across, 67 marks      panel B   0.421 mm, 445
        panel D   0.518 mm, 69                   panel E   0.424 mm, 773
        the dashed regression rule               0.595 pt

    The constants below are those numbers converted through this panel's own
    MARK. Not one coordinate, colour or statistic moves.
"""

import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCGA_BAYESPRISM_EPI, TCGA_BULK_TUMOR_ONLY, TCGA_RAW_DIR, TCGA_CLINICAL
import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402
from cnsfig import layout as cnslayout, corr_stats, rich_xlabel, rich_ylabel  # noqa: E402
from cnsfig import cache, group_key  # noqa: E402

# ==============================================================================
# Configuration
# ==============================================================================
SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 4.5                  # the earlier smallest body type, before * SCALE

MARK = style.tick_pt() / (SMALL_PT * SCALE)   # length multiplier
AREA = MARK ** 2                              # area multiplier

PANEL_LETTER = "E"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(3, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(3, PANEL_LETTER)
#: The column's margin, shared with the panel above (see 03_A / 03_B).
LEFT_MM = 12.5   # box at the same page x as the panel above (2026-09-15)
OUTPUT_DIR = Path(__file__).parent

# Data paths (from paths.py)
TCGA_EPI_EXPR = TCGA_BAYESPRISM_EPI

# Tex signature (18 genes)
TEX_GENES = [
    'PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'CTLA4', 'TOX', 'ENTPD1', 'CXCL13',
    'LAYN', 'CD38', 'BATF', 'IRF4', 'PRDM1', 'TOX2', 'ITGAE', 'NR4A1', 'NR4A2', 'NR4A3',
]

# Survival dot colors
COLOR_DEAD = '#c0392b'      # dark red, filled
COLOR_ALIVE = '#888888'     # gray edge, unfilled
COLOR_MISSING = '#cccccc'   # light gray edge, unfilled


# ==============================================================================
# ssGSEA scoring (proper weighted KS, from original validation session)
# ==============================================================================
def ssgsea_score(expr_df, gene_set):
    """Proper ssGSEA: weighted Kolmogorov-Smirnov enrichment. expr_df: samples x genes."""
    gene_set = [g for g in gene_set if g in expr_df.columns]
    print(f"  ssGSEA: {len(gene_set)}/{len(TEX_GENES)} Tex genes found")

    scores = []
    for sample in expr_df.index:
        profile = expr_df.loc[sample].dropna()
        ranked = profile.rank(ascending=True)
        sorted_idx = ranked.sort_values().index
        in_set_sorted = sorted_idx.isin(gene_set)
        ranks_sorted = ranked[sorted_idx]

        hits = in_set_sorted.astype(float)
        hw = hits * np.abs(ranks_sorted)
        hs = hw.sum()
        if hs == 0:
            scores.append(0)
            continue
        hr = np.cumsum(hw / hs)

        mw = (1 - hits)
        ms = mw.sum()
        if ms == 0:
            scores.append(0)
            continue
        mr = np.cumsum(mw / ms)
        scores.append((hr - mr).sum())

    return pd.Series(scores, index=expr_df.index)


# ==============================================================================
# Data loading
# ==============================================================================
def load_data():
    # 1. Deconvolved epithelial expression (samples × genes)
    print("[1/2] Loading TCGA BayesPrism epithelial expression...")
    epi = pd.read_csv(TCGA_EPI_EXPR, sep='\t', index_col=0)
    print(f"  Shape: {epi.shape}")

    # 2. Tumor-only FPKM (407 samples) + gene mapping from STAR raw files
    print("[2/2] Loading TCGA tumor-only FPKM for ssGSEA...")
    fpkm = pd.read_csv(TCGA_BULK_TUMOR_ONLY, sep='\t', index_col=0)
    fpkm = fpkm[fpkm.index.str.startswith('ENSG')]
    print(f"  FPKM: {fpkm.shape[0]:,} genes × {fpkm.shape[1]} samples")

    # Gene mapping from raw STAR annotation.
    #
    # This was `os.listdir(TCGA_RAW_DIR)[0]` - unfiltered and unsorted - until
    # 2026-09-10. TCGA_RAW_DIR holds 443 STAR files and a .gitkeep, and in an
    # extracted copy of the external tree the placeholder came back first, so
    # the panel died with pandas EmptyDataError: No columns to parse from file.
    # It was order-dependent, which is worse than simply broken: on the
    # filesystem the published panels were built on the same expression
    # returned a real file and the panel worked.
    #
    # Choosing differently cannot move a number, and that is measured rather
    # than assumed. The file is read for one thing - the Ensembl-ID to
    # gene-symbol map in `gene_id`/`gene_name`; the per-sample counts are in
    # columns `usecols` discards. Across all 443 STAR files that map is
    # byte-identical (md5 b41b583e25d0c49a5ba00c39702d857f of the sorted
    # id\tsymbol pairs, 1 distinct map). So a filtered, sorted, deterministic
    # choice reproduces the published panel whichever file the original run got.
    #
    # The glob is the one 04_Revision_Analyses/05_R1.6_Spatial_Confounders/
    # scripts/tcga_immune_exclusion.py already uses on this same directory.
    star = sorted(TCGA_RAW_DIR.glob('*.augmented_star_gene_counts.tsv'))
    if not star:
        raise SystemExit(
            f"no *.augmented_star_gene_counts.tsv in {TCGA_RAW_DIR}. The "
            f"Ensembl-to-symbol map is read from one of these; refusing to "
            f"guess at another file in the directory.")
    ann = pd.read_csv(star[0], sep='\t', skiprows=1,
                       usecols=['gene_id', 'gene_name']).dropna()
    ens2sym = dict(zip(ann['gene_id'], ann['gene_name']))

    fpkm.index = fpkm.index.map(lambda x: ens2sym.get(x, x))
    fpkm = fpkm[~fpkm.index.duplicated(keep='first')]

    # Samples × genes, log2(FPKM+1)
    fpkm_log = np.log2(fpkm.T + 1)
    print(f"  After mapping: {fpkm_log.shape[1]:,} genes, {fpkm_log.shape[0]} samples")

    # ssGSEA
    print("  Running ssGSEA...")
    tex_scores = ssgsea_score(fpkm_log, TEX_GENES)

    # Align: epi (410 samples) ∩ FPKM (407 tumor-only)
    common = sorted(set(epi.index) & set(tex_scores.index))
    print(f"  Common samples: {len(common)}")

    df = pd.DataFrame({
        'sample': common,
        'CEACAM5': np.log2(epi.loc[common, 'CEACAM5'].values + 1),
        'CEACAM6': np.log2(epi.loc[common, 'CEACAM6'].values + 1),
        'CD274': np.log2(epi.loc[common, 'CD274'].values + 1),
        'tex_ssgsea_raw': tex_scores.loc[common].values,
    })

    # Z-score the Tex ssGSEA (z-score)s
    mu = df['tex_ssgsea_raw'].mean()
    sd = df['tex_ssgsea_raw'].std()
    df['tex_ssgsea'] = (df['tex_ssgsea_raw'] - mu) / sd
    print(f"  Tex ssGSEA z-scored: mean={mu:.2f}, sd={sd:.2f}")

    # Merge clinical survival data
    print("  Merging clinical vital_status...")
    clin = pd.read_csv(TCGA_CLINICAL, sep='\t', usecols=['submitter_id', 'vital_status'])
    clin = clin.drop_duplicates(subset='submitter_id')
    df = df.merge(clin, left_on='sample', right_on='submitter_id', how='left')
    df['vital_status'] = df['vital_status'].fillna('Missing')
    print(f"  Vital status counts: {df['vital_status'].value_counts().to_dict()}")

    return df[['CEACAM5', 'CEACAM6', 'CD274', 'tex_ssgsea',
               'vital_status']].reset_index(drop=True)


# ==============================================================================
# Plotting
# ==============================================================================
def plot_scatter(ax, x, y, vital, title, xlabel, ylabel, fontscale, left=True):
    r_val, p_val = stats.spearmanr(x, y)

    # Plot dots by vital_status: Dead filled red, Alive open gray, Missing open light gray
    for status, fc, ec, label in [
        ('Dead',    COLOR_DEAD,  COLOR_DEAD,    'Dead'),
        ('Alive',   'none',      COLOR_ALIVE,   'Alive'),
        ('Missing', 'none',      COLOR_MISSING, None),
    ]:
        mask = vital == status
        if mask.sum() == 0:
            continue
        ax.scatter(x[mask], y[mask], facecolors=fc, edgecolors=ec,
                   s=13.0*AREA, alpha=0.85, linewidths=style.EDGE_PT, zorder=3)

    # Regression line (above dots)
    valid = np.isfinite(x) & np.isfinite(y)
    xv, yv = np.array(x)[valid], np.array(y)[valid]
    if len(xv) >= 3 and np.std(xv) > 0:
        slope, intercept = np.polyfit(xv, yv, 1)
        x_line = np.linspace(xv.min(), xv.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=style.RULE_PT, alpha=0.6, zorder=5)

    p_str = corr_stats.p_string(p_val)
    print(f'    {title}: rho = {r_val:.2f}, {p_str}')
    ax.set_title(title, fontsize=style.tick_pt())
    cnslayout.corr_annotate(ax, r_val)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(style.RULE_PT)
    ax.tick_params(axis='both', width=style.RULE_PT, length=3*fontscale*MARK)
    # The rich labels last: they measure the tick labels' reach, which the
    # tick length above changes.
    rich_xlabel(ax, xlabel)
    if left:
        rich_ylabel(ax, ylabel, y=0.40)   # 1.3 mm down: clear of the 10 pt letter cell (2026-09-15)

    # Compact legend (Dead filled red, Alive open gray). `key` undoes cnsplots'
    # legend.markerscale so each key prints at MARK times the size it had.
    key = 1.0 / plt.rcParams['legend.markerscale']
    h_dead = mlines.Line2D([], [], marker='o', color='none', markerfacecolor=COLOR_DEAD,
                           markeredgecolor=COLOR_DEAD,
                           markersize=4*fontscale/SCALE*MARK*key,
                           label='Dead')
    h_alive = mlines.Line2D([], [], marker='o', color='none', markerfacecolor='none',
                            markeredgecolor=COLOR_ALIVE,
                            markersize=4*fontscale/SCALE*MARK*key,
                            markeredgewidth=style.EDGE_PT, label='Alive')
    # RETURNED, NOT DRAWN, SINCE 2026-09-11
    #   These two keys were drawn inside each axes at loc='upper left', which
    #   put the words 'Dead' and 'Alive' on top of the scatter in both axes and
    #   twice over. The published page prints one key, outside the plotting
    #   area, at the right of the panel - which is also where panel D of this
    #   figure prints its five-group key. main() draws it once.
    return r_val, p_val, [h_dead, h_alive]


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    # From data/tcga_samples.csv (cnsfig.cache, 2026-09-15): the expression,
    # deconvolution and clinical tables are read only when it is absent.
    df = cache.table(OUTPUT_DIR, 'tcga_samples', load_data)
    n = len(df)

    # CEACAM vs Tex ssGSEA only (filter zero-CEACAM samples per subplot)
    panel = {
        'name': 'tcga_ceacam_tex',
        'combos': [
            ('CEACAM5', 'tex_ssgsea', '*CEACAM5*', 'Tex ssGSEA\n(z-score)'),
            ('CEACAM6', 'tex_ssgsea', '*CEACAM6*', 'Tex ssGSEA\n(z-score)'),
        ],
    }

    # ONE GEOMETRY FOR THE FOUR SCATTER PAIRS  (2026-09-14, evening)
    #   The author's ruling: A, B, D and E are the same size, square, the
    #   dataset name alone in the title, rho inside the box, P in the legend
    #   (cnsfig.corr_stats writes it; edits.py reads it), and the two rows
    #   2 mm apart. cnsfig.layout.scatter_pair_mm places the boxes at
    #   millimetres, so the four panels print one geometry by construction.
    results = []
    fig = style.figure_mm(PANEL_W_MM, PANEL_H_MM)
    axes = cnslayout.scatter_pair_mm(fig, left_mm=LEFT_MM)
    for col_idx, (xcol, ycol, xlabel, ylabel) in enumerate(panel['combos']):
        ax = axes[col_idx]
        # Remove samples with zero CEACAM expression
        mask = df[xcol] > 0
        sub = df[mask]
        n_sub = len(sub)
        cohort_title = 'TCGA-STAD'   # the page prints no n (2026-09-14)
        r_val, p_val, keys = plot_scatter(ax, sub[xcol].values, sub[ycol].values,
                                          sub['vital_status'].values,
                                          cohort_title, xlabel, ylabel, SCALE,
                                          left=(col_idx == 0))
        results.append({'x': xcol, 'y': ycol, 'r': r_val, 'p': p_val, 'n': n_sub})
    corr_stats.write(OUTPUT_DIR, [(r['x'], r['r'], r['p'], r['n']) for r in results])

    # The key in the strip at the right edge, as printed panel D of this
    # figure prints its five-group key.
    # The key, right of the second box, at cnsfig.legend.group_key's fixed
    # 1.3 mm circles (the handle-sized keys of 2026-09-15 were 0.42 mm).
    group_key(fig, [(keys[0].get_label(), COLOR_DEAD),
                    (keys[1].get_label(), None, COLOR_ALIVE)],
              x_mm=LEFT_MM + 2 * cnslayout.SCATTER_BOX_MM + cnslayout.SCATTER_GAP_MM + 0.8,
              y_mm=cnslayout.SCATTER_TOP_MM)

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    style.save_panel(fig, OUTPUT_DIR / panel['name'])
    print(f"Saved: {OUTPUT_DIR / panel['name']}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")

    # Summary
    print(f"\n{'='*65}")
    print(f"TCGA-STAD CORRELATION SUMMARY (Spearman)")
    print(f"{'='*65}")
    print(f"{'X':<12} {'Y':<18} {'n':>5} {'r':>6} {'P':>12} {'Sig':>5}")
    print(f"{'-'*65}")
    for r in results:
        sig = '***' if r['p'] <= 0.001 else '**' if r['p'] <= 0.01 else '*' if r['p'] <= 0.05 else 'ns'
        print(f"{r['x']:<12} {r['y']:<18} {r['n']:>5} {r['r']:>6.3f} {r['p']:>12.6f} {sig:>5}")
    print(f"{'='*65}")


if __name__ == '__main__':
    main()
