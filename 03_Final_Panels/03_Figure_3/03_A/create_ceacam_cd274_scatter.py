#!/usr/bin/env python3
"""
Figure 3 panel A - CEACAM5 and CEACAM6 against CD274, one point per stomach
sample of the primary cohort, coloured by treatment phase and response.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 3 A       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

Every value read, every filter, every statistic and every string is the earlier
drawing's. The Grubbs' outlier exclusion, the four-group assignment, the
Spearman test and the P-value formatting are untouched, and expression is read
exactly as it was read before.

THE COUNTS COME FROM THE EPITHELIAL OBJECT'S OWN LAYER
    The counts were read from a second file holding the same 149,373 cells and
    the same 56,034 genes with the genes in alphabetical order. They are
    layers['counts'] of the epithelial object here, in the gene order .X uses.
    The two agree over every one of their 203,799,816 non-zeros, so the library
    size each cell normalises by is a sum over the same values, and every gene
    below is selected by name rather than by position. The order the counts are
    stored in therefore reaches no number this panel prints.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the tick labels - at 5 * SCALE. MARK carries the
    non-type point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The marker area, the marker edge and the dashed regression line are scaled
    by it. Tick widths and lengths and spine widths are not: those are style,
    and cnsplots sets them.

THE COHORT NAME IS BACK IN THE TITLE, AS THE PUBLISHED PAGE PRINTS IT
    It was dropped on 2026-09-10 because it would not fit: at 23.1 mm of panel
    height there was no room for a second title line, and the reasoning
    recorded here was that the caption names the cohort anyway. That was the
    wrong trade. The published page prints the dataset name over every one of
    these eight scatter axes, the reader needs it to tell an in-house panel
    from a TCGA one at a glance, and the fix belonged in the frame rather than
    in the drawing.

    Since 2026-09-11 the panel is 38.0 mm tall (panel_rects_v2.csv), and the
    published two-line title fits: 'In house (scRNA-seq)' sets 21.0 mm at 6 pt
    against a 24.3 mm plotting column. The title is left-aligned and set at
    normal weight, both as the page prints it.

    The y label is re-wrapped onto two lines and keeps every word: it is
    rotated, so its length is vertical, and on one line it sets 35.8 mm against
    a 23.1 mm panel. The x label carries the gene alone; at 7 pt the unit sets
    18.0 mm and the two plotting boxes are about 10 mm wide, so the two axes'
    units would meet between them. The unit is stated on the y label of this
    panel and in the caption. Both are declared in RENAMES_FIGURE_3.

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

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *
import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402
from cnsfig import layout as cnslayout, corr_stats, rich_xlabel, rich_ylabel  # noqa: E402
from cnsfig import cache, group_key  # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier smallest body type, before * SCALE

MARK = style.tick_pt() / (SMALL_PT * SCALE)   # length multiplier
AREA = MARK ** 2                              # area multiplier

PANEL_LETTER = "A"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(3, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(3, PANEL_LETTER)
#: The left column's margin (A over B): a one-line y label and 4-glyph ticks.
LEFT_MM = 11.8              # box at x = 19.0 mm on the page, as B's
#: The five-group key's top-left, mm from the canvas corner: right of the
#: second box, level with its top.
KEY_X_MM = LEFT_MM + 2 * cnslayout.SCATTER_BOX_MM + cnslayout.SCATTER_GAP_MM + 0.8
KEY_Y_MM = cnslayout.SCATTER_TOP_MM

EPI = EPITHELIAL_DEPOSIT_H5AD
OUTPUT_DIR = Path(__file__).parent

# Standard 4-group palette (blue=R, warm=NR; lighter=Pre, darker=Post)
COLOR_MAP_4GROUP = {
    'Pre-R':  '#74b9ff',
    'Pre-NR': '#e17055',
    'Post-R': '#0984e3',
    'Post-NR':'#d63031',
    'Other':  '#999999',
}



def load_primary():
    """Load primary cohort: ALL 32 stomach samples (pre + post)."""
    print("[Primary] Loading epithelial counts...")
    adata = sc.read_h5ad(EPI)
    adata.X = adata.layers['counts']
    del adata.layers['counts']
    # uns['log1p'] describes the matrix just replaced. Left in place
    # it makes the normalisation below report the counts as already
    # log-transformed.
    adata.uns.pop('log1p', None)

    # Filter: stomach only (all treatment phases)
    stomach_mask = adata.obs['Sample site'].astype(str).str.lower().str.contains('stomach')
    adata = adata[stomach_mask].copy()

    # Normalize raw counts
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Extract expression
    genes = ['CD274', 'CEACAM5', 'CEACAM6']
    expr = {}
    for g in genes:
        x = adata[:, g].X
        expr[g] = x.toarray().flatten() if hasattr(x, 'toarray') else x.flatten()

    df = pd.DataFrame(expr)
    df['sample'] = adata.obs['sample'].values
    df['study_id'] = adata.obs['Sample ID'].values
    df['treatment_phase'] = adata.obs['Treatment phase'].values
    df['pre_group'] = adata.obs['stomach_pre_grouping'].values
    df['post_group'] = adata.obs['stomach_post_grouping'].values

    # Sample-level aggregation
    agg = df.groupby('sample').agg({
        'CD274': 'mean', 'CEACAM5': 'mean', 'CEACAM6': 'mean',
        'treatment_phase': 'first',
        'pre_group': 'first',
        'post_group': 'first',
        'study_id': 'first',
    }).reset_index()

    # 4-group assignment
    def assign_4group(row):
        phase = str(row['treatment_phase'])
        if phase == 'Pre':
            grp = str(row['pre_group'])
            if grp == 'Responsed': return 'Pre-R'
            if grp == 'No-response': return 'Pre-NR'
        elif phase == 'Post':
            grp = str(row['post_group'])
            if grp == 'Responsed': return 'Post-R'
            if grp == 'No-response': return 'Post-NR'
        return 'Other'

    agg['group'] = agg.apply(assign_4group, axis=1)

    # One sample is a statistical outlier on CD274 (Grubbs' test G = 5.43,
    # P < 0.05) and is excluded from the correlation. It is named by its study
    # ID from Supplementary Table 1, not by the internal specimen number the
    # 'sample' column carries: that number identifies a specimen in the hospital
    # record and does not belong in deposited code. The two are 1:1 across all
    # 70 specimens, so this selects the same row.
    OUTLIER_STUDY_ID = 'P32-P1'
    if OUTLIER_STUDY_ID not in set(agg['study_id']):
        raise SystemExit(f"{OUTLIER_STUDY_ID} is not among the stomach samples - "
                         "the outlier exclusion would silently do nothing")
    agg = agg[agg['study_id'] != OUTLIER_STUDY_ID].reset_index(drop=True)

    counts = agg['group'].value_counts()
    print(f"  {len(agg)} stomach samples (after outlier removal): {counts.to_dict()}")
    # The cached table names samples by study ID only (see OUTLIER_STUDY_ID).
    return agg[['study_id', 'group', 'CD274', 'CEACAM5', 'CEACAM6']].reset_index(drop=True)


def load_tiger():
    """Load TIGER PRJEB25780: 78 ICB-treated GC (all pre-treatment)."""
    print("[TIGER] Loading expression...")
    expr = pd.read_csv(TIGER_EXPR, sep='\t', index_col=0)
    genes = ['CD274', 'CEACAM5', 'CEACAM6']

    tiger = expr.loc[genes].T.copy()
    tiger.columns = genes
    tiger = np.log2(tiger + 1)

    # Metadata
    meta = pd.read_csv(TIGER_META, sep='\t')
    id_col = None
    for col in meta.columns:
        if meta[col].dtype == object and meta[col].isin(tiger.index).sum() > 10:
            id_col = col
            break
    if id_col:
        meta_indexed = meta.set_index(id_col)
        tiger = tiger.join(meta_indexed['response'], how='left')

    def map_resp(val):
        val = str(val)
        if val in ['CR', 'PR']: return 'Pre-R'
        if val in ['SD', 'PD']: return 'Pre-NR'
        return 'Other'
    tiger['group'] = tiger['response'].map(map_resp)
    tiger = tiger.reset_index().rename(columns={'index': 'sample'})

    print(f"  {len(tiger)} samples — Pre-R:{(tiger.group=='Pre-R').sum()}, Pre-NR:{(tiger.group=='Pre-NR').sum()}")
    return tiger


def plot_scatter(ax, x, y, groups, color_map, draw_order, title, xlabel, ylabel, fontscale, left=True):
    """Generic scatter with Spearman stats."""
    r_val, p_val = stats.spearmanr(x, y)

    # MARKER AREA IS MEASURED OFF THE PUBLISHED PAGE, NOT CHOSEN  (2026-09-11)
    #   The published panel draws these points at 0.52 mm across - counted
    #   straight out of 'Figure 3.pdf' geometry, 67 of them. The redraw was
    #   drawing them at 1.16 mm, which is 2.2 times as wide and enough
    #   to merge neighbouring samples into one blob. Area is diameter squared,
    #   so the constant below is (0.52 mm / MARK-scaled point) squared. Not one
    #   coordinate moves.
    for grp in draw_order:
        mask = groups == grp
        if mask.sum() > 0:
            ax.scatter(x[mask], y[mask], c=color_map[grp], s=24.1*AREA,
                       alpha=0.85, edgecolors='white', linewidths=style.EDGE_PT,
                       label=grp, zorder=3)

    # Regression line
    valid = np.isfinite(x) & np.isfinite(y)
    xv, yv = x[valid], y[valid]
    if len(xv) >= 3 and xv.std() > 0:
        slope, intercept = np.polyfit(xv, yv, 1)
        x_line = np.linspace(xv.min(), xv.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=style.RULE_PT, alpha=0.6, zorder=2)

    p_str = corr_stats.p_string(p_val)
    print(f"    {title}: rho = {r_val:.2f}, {p_str}")
    ax.set_title(title, fontsize=style.tick_pt())
    cnslayout.corr_annotate(ax, r_val)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(style.RULE_PT)
    ax.tick_params(axis='both', width=style.RULE_PT, length=3*SCALE*MARK)
    # The rich labels last: they measure the tick labels' reach, which the
    # tick length above changes.
    rich_xlabel(ax, f"*{xlabel}*")
    if left:
        rich_ylabel(ax, ylabel, y=0.32)   # 2.4 mm down: clear of the 10 pt letter cell
    return r_val, p_val

    # NOT set_box_aspect(1). The published panel's plotting boxes are wider
    # than they are tall, and a square box is width-limited here, so it would
    # leave the 15 mm this panel gained on 2026-09-11 as blank paper under the
    # axes instead of putting it into the plot.


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    # The per-sample means, from data/primary_samples.csv (cnsfig.cache,
    # 2026-09-15): the h5ad is read only when the table is absent.
    primary = cache.table(OUTPUT_DIR, 'primary_samples', load_primary)

    # ONE GEOMETRY FOR THE FOUR SCATTER PAIRS  (2026-09-14, evening)
    #   The author's ruling: A, B, D and E are the same size, square, the
    #   dataset name alone in the title, rho inside the box, P in the legend
    #   (cnsfig.corr_stats writes it; edits.py reads it), and the two rows
    #   2 mm apart. cnsfig.layout.scatter_pair_mm places the boxes at
    #   millimetres, so the four panels print one geometry by construction.
    fig = style.figure_mm(PANEL_W_MM, PANEL_H_MM)
    axes = cnslayout.scatter_pair_mm(fig, left_mm=LEFT_MM)

    # Draw order: Other first (background), then colored groups on top
    primary_order = ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']

    rows = []
    for col_idx, ceacam in enumerate(['CEACAM5', 'CEACAM6']):
        axes[col_idx].set_ylim(0, 0.03)
        r_val, p_val = plot_scatter(
            ax=axes[col_idx],
            x=primary[ceacam].values.astype(float),
            y=primary['CD274'].values.astype(float),
            groups=primary['group'].values,
            color_map=COLOR_MAP_4GROUP,
            draw_order=primary_order,
            title='In house',
            xlabel=ceacam,
            # One line (2026-09-14, evening): the unit line '(mean log expr.)'
            # ran into the panel-letter cell on a 13.5 mm box; the legend
            # states the unit for both axes. Declared in REMOVALS_FIGURE_3.
            ylabel='PD-L1 (*CD274*)',
            fontscale=SCALE,
            left=(col_idx == 0),
        )
        rows.append((ceacam, r_val, p_val, int(np.isfinite(primary[ceacam].values).sum())))
    corr_stats.write(OUTPUT_DIR, rows)

    # The five-group key, at the right of the second box, as the published
    # panel prints it. Dropped from the redraw without a declaration and put
    # back on 2026-09-15; the circles are cnsfig.legend.group_key's fixed
    # 1.3 mm, not the 0.5 mm data marker.
    group_key(fig, [(g, COLOR_MAP_4GROUP[g]) for g in
                    ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']],
              x_mm=KEY_X_MM, y_mm=KEY_Y_MM)

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, OUTPUT_DIR / 'ceacam_cd274_scatter')
    print(f"\nSaved: {OUTPUT_DIR / 'ceacam_cd274_scatter'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
