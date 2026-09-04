#!/usr/bin/env python3
"""
Figure 3 panel A, RESTYLED (Version B) - CEACAM5/6 vs CD274 correlation (1x2),
primary cohort only (scRNA, ALL 32 stomach samples, 4-group colouring).

Version A is
`03_Revised_Panels/Main_Figures/03_Figure_3/03_A/create_ceacam_cd274_scatter.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every statistic and every string is Version A's.
The drawing code is the same code. In particular the Grubbs' outlier exclusion,
the 4-group assignment, the Spearman test and the P-value formatting are
untouched, and expression is read exactly as Version A reads it.

  printed panel  Figure 3 A     (PROVENANCE.csv; NOT inferred from "03_A")
  printed rect   48.7 x 23.1 mm   (panel_rects.csv)
  Version B box  107.0 x 57.0 mm

MARK
    Version A drew a 32.0 x 16.0 cm canvas (SCALE = 4) and its smallest body
    type is the tick labels at `labelsize=5 * fontscale`, fontscale = SCALE.
    So SMALL_PT = 5 and, by PANEL_SPEC.md,

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2                    = 0.1225

    Every non-type length Version A set - marker area 30*SCALE, marker edge
    0.3*SCALE, the dashed regression line 0.8*SCALE, the spines 0.5*SCALE and
    the tick width/length 0.5*SCALE / 3*SCALE - is multiplied by that, so each
    keeps its Version A size *relative to the type*.

Why the box grew from 48.7 mm to 107 mm: the two axes are square
(`set_box_aspect(1)`, Version A's) and each carries a y axis label, its tick
labels and a two-line title. At 7/8 pt that furniture is ~27 mm of the width
before a single point is drawn, and the title alone sets ~40 mm across each
axes. Nothing about what is plotted changed.
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
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *
import panel_style_cns as style  # noqa: E402

# Version A's canvas multiplier. Version B draws 1:1, so SCALE survives only to
# reproduce the exact numbers Version A set for its non-type lengths.
SCALE = 4
SMALL_PT = 5.0                       # Version A's smallest body type

# See MARK above. `style.tick_pt()` reads cnsplots' own setting, so the factor
# is derived, never a literal.
MARK = style.tick_pt() / (SMALL_PT * SCALE)   # 0.35, length multiplier
AREA = MARK ** 2                              # 0.1225, area multiplier

PRINTED_MM = (48.7, 23.1)            # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 107.0, 57.0
# Millimetres of paper. left holds the y label + its tick labels, top the
# two-line title, bottom the x label + its tick labels; wspace 1/3 of a column
# leaves 13 mm between the axes for the right subplot's own y furniture.
MARGIN = dict(left=14.0, right=2.0, top=9.0, bottom=9.0, wspace=1.0 / 3.0)

EPI_RAW = EPITHELIAL_RAW_COUNTS_H5AD
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
    print("[Primary] Loading epithelial raw counts...")
    adata = sc.read_h5ad(EPI_RAW)

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

    # Stats text — 1 sig digit (floor), scientific for very small P
    import math
    _e = math.floor(math.log10(p_val)); _c = int(p_val / 10**_e)
    if _e >= -3:
        p_str = f'P = {_c * 10**_e:.{-_e}f}'
    else:
        p_str = f'P = {_c}' + r'$\times 10^{' + str(_e) + r'}$'

    ax.set_title(f'{title}\nρ = {r_val:.2f}, {p_str}', linespacing=1.4)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(0.5*SCALE*MARK)
    ax.tick_params(axis='both', width=0.5*SCALE*MARK, length=3*SCALE*MARK)

    ax.set_box_aspect(1)


def main():
    family = style.apply()
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
            xlabel=f'{ceacam} (mean log expr.)',
            ylabel='PD-L1 (CD274) (mean log expr.)',
            fontscale=SCALE,
        )
        axes[col_idx].set_ylim(0, 0.03)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUTPUT_DIR / 'ceacam_cd274_scatter')
    print(f"\nSaved: {OUTPUT_DIR / 'ceacam_cd274_scatter'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
