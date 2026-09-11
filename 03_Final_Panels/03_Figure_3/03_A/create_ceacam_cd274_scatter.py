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

THE COHORT LEAVES THE TITLE, AND THE AXIS LABELS ARE RE-WRAPPED
    The two axes carry the same cohort name, and so do the two axes of panel D;
    printed four times over it is four titles saying one thing. The cohort is
    named once in the caption and each axes keeps its own correlation and its
    own P value, which are the numbers that differ between them. On one line
    the pair sets 23.6 mm against a plotting box of about 10 mm, so they take a
    line each.

    The y label is re-wrapped onto two lines and keeps every word: it is
    rotated, so its length is vertical, and on one line it sets 35.8 mm against
    a 23.1 mm panel. The x label carries the gene alone; at 7 pt the unit sets
    18.0 mm and the two plotting boxes are about 10 mm wide, so the two axes'
    units would meet between them. The unit is stated on the y label of this
    panel and in the caption. Both are declared in RENAMES_FIGURE_3.
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

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier smallest body type, before * SCALE

MARK = style.tick_pt() / (SMALL_PT * SCALE)   # length multiplier
AREA = MARK ** 2                              # area multiplier

PANEL_LETTER = "A"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(3, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(3, PANEL_LETTER)

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
    return agg


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


def plot_scatter(ax, x, y, groups, color_map, draw_order, title, xlabel, ylabel, fontscale):
    """Generic scatter with Spearman stats."""
    r_val, p_val = stats.spearmanr(x, y)

    for grp in draw_order:
        mask = groups == grp
        if mask.sum() > 0:
            ax.scatter(x[mask], y[mask], c=color_map[grp], s=30*SCALE*AREA,
                       alpha=0.85, edgecolors='white', linewidths=0.3*SCALE*MARK,
                       label=grp, zorder=3)

    # Regression line
    valid = np.isfinite(x) & np.isfinite(y)
    xv, yv = x[valid], y[valid]
    if len(xv) >= 3 and xv.std() > 0:
        slope, intercept = np.polyfit(xv, yv, 1)
        x_line = np.linspace(xv.min(), xv.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=0.8*SCALE*MARK, alpha=0.6, zorder=2)

    # Stats text — 1 sig digit (floor). A very small P is written out rather
    # than set as a mathtext power of ten: mathtext draws a superscript at 70%
    # of its base, so an exponent on a 7 pt title prints at 4.9 pt, below the
    # floor this figure set is set to.
    import math
    _e = math.floor(math.log10(p_val)); _c = int(p_val / 10**_e)
    if _e >= -3:
        p_str = f'P = {_c * 10**_e:.{-_e}f}'
    else:
        p_str = f'P = {_c}e{_e}'

    # The cohort is named in the caption; see THE COHORT LEAVES THE TITLE.
    print(f"    {title}: rho = {r_val:.2f}, {p_str}")
    ax.set_title(f'ρ = {r_val:.2f}\n{p_str}', linespacing=1.4)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(0.5*SCALE*MARK)
    ax.tick_params(axis='both', width=0.5*SCALE*MARK, length=3*SCALE*MARK)

    ax.set_box_aspect(1)


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    primary = load_primary()

    # 1x2 layout (primary cohort only)
    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, 1, 2)

    # Draw order: Other first (background), then colored groups on top
    primary_order = ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']

    for col_idx, ceacam in enumerate(['CEACAM5', 'CEACAM6']):
        plot_scatter(
            ax=axes[col_idx],
            x=primary[ceacam].values.astype(float),
            y=primary['CD274'].values.astype(float),
            groups=primary['group'].values,
            color_map=COLOR_MAP_4GROUP,
            draw_order=primary_order,
            title='Primary Cohort (scRNA-seq)',
            xlabel=ceacam,
            ylabel='PD-L1 (CD274)\n(mean log expr.)',
            fontscale=SCALE,
        )
        axes[col_idx].set_ylim(0, 0.03)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
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
