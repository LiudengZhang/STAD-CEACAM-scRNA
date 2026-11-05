#!/usr/bin/env python3
"""
Figure 3A: CEACAM5 and CEACAM6 against CD274 (PD-L1), one point per stomach
sample of the in-house cohort, coloured by treatment phase and response.
Spearman rho and P are printed in each box and written to
correlation_stats.csv (cnsfig.corr_stats). One CD274 outlier (Grubbs' test) is
excluded.

Reads data/primary_samples.csv beside this script; with --recompute, computes
it from the deposited epithelial h5ad (EPITHELIAL_DEPOSIT_H5AD), normalising
layers['counts'] to log1p(CP10k).
Writes ceacam_cd274_scatter.[svg|pdf|png] beside this script.

Run: python create_ceacam_cd274_scatter.py [--recompute]
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_DEPOSIT_H5AD
from cnsfig import layout as cnslayout, corr_stats, rich_xlabel, rich_ylabel  # noqa: E402
from cnsfig import cache, group_key, Panel, style  # noqa: E402

PANEL = Panel(3, "A", __file__, small_pt=5.0)
MARK, AREA = PANEL.mark, PANEL.area
SCALE = Panel.SCALE
# Left margin, shared with the panel below so the boxes align.
LEFT_MM = 11.8
# Group key: right of the second box, level with its top (mm).
KEY_X_MM = LEFT_MM + 2 * cnslayout.SCATTER_BOX_MM + cnslayout.SCATTER_GAP_MM + 0.8
KEY_Y_MM = cnslayout.SCATTER_TOP_MM

EPI = EPITHELIAL_DEPOSIT_H5AD

# Shared four-group palette, plus grey for samples outside the four groups.
from shared.figure_config import GROUP_COLORS             # noqa: E402
COLOR_MAP_4GROUP = {**GROUP_COLORS, 'Other': '#999999'}


def load_primary():
    """Per-sample mean expression for all stomach samples (pre and post)."""
    adata = sc.read_h5ad(EPI)
    adata.X = adata.layers['counts']
    del adata.layers['counts']
    # Drop uns['log1p'], which describes the replaced matrix; otherwise scanpy
    # treats the counts as already log-transformed.
    adata.uns.pop('log1p', None)

    stomach_mask = adata.obs['Sample site'].astype(str).str.lower().str.contains('stomach')
    adata = adata[stomach_mask].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

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

    agg = df.groupby('sample').agg({
        'CD274': 'mean', 'CEACAM5': 'mean', 'CEACAM6': 'mean',
        'treatment_phase': 'first',
        'pre_group': 'first',
        'post_group': 'first',
        'study_id': 'first',
    }).reset_index()

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

    # One sample is an outlier on CD274 (Grubbs' test G = 5.43, P < 0.05) and
    # is excluded. It is named by its Supplementary Table 1 study ID.
    OUTLIER_STUDY_ID = 'P32-P1'
    if OUTLIER_STUDY_ID not in set(agg['study_id']):
        raise SystemExit(f"{OUTLIER_STUDY_ID} is not among the stomach samples - "
                         "the outlier exclusion would silently do nothing")
    agg = agg[agg['study_id'] != OUTLIER_STUDY_ID].reset_index(drop=True)

    # The cached table names samples by study ID only.
    return agg[['study_id', 'group', 'CD274', 'CEACAM5', 'CEACAM6']].reset_index(drop=True)


def plot_scatter(ax, x, y, groups, color_map, draw_order, title, xlabel, ylabel, fontscale, left=True):
    """Generic scatter with Spearman stats."""
    r_val, p_val = stats.spearmanr(x, y)

    # Marker about 0.52 mm across.
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

    ax.set_title(title, fontsize=style.tick_pt())
    cnslayout.corr_annotate(ax, r_val, p=p_val)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(style.RULE_PT)
    ax.tick_params(axis='both', width=style.RULE_PT, length=3*SCALE*MARK)
    # Axis labels last: they are placed from the tick labels' extent.
    rich_xlabel(ax, f"*{xlabel}*")
    if left:
        rich_ylabel(ax, ylabel, y=0.32)   # lowered, clear of the panel letter
    return r_val, p_val


def main():
    primary = cache.table(PANEL.dir, 'primary_samples', load_primary)

    # Two square boxes placed in mm (shared with Figures 3B, 3D and 3E).
    fig = PANEL.canvas()
    axes = cnslayout.scatter_pair_mm(fig, left_mm=LEFT_MM)

    # Draw order: Other first (background), then coloured groups on top
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
            ylabel='PD-L1 (*CD274*)',
            fontscale=SCALE,
            left=(col_idx == 0),
        )
        rows.append((ceacam, r_val, p_val, int(np.isfinite(primary[ceacam].values).sum())))
    corr_stats.write(PANEL.dir, rows)

    # Group key to the right of the second box.
    group_key(fig, [(g, COLOR_MAP_4GROUP[g]) for g in
                    ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']],
              x_mm=KEY_X_MM, y_mm=KEY_Y_MM)

    PANEL.finish(fig, 'ceacam_cd274_scatter')


if __name__ == '__main__':
    main()
