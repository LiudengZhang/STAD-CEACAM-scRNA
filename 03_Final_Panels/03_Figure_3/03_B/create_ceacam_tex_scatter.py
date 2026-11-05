#!/usr/bin/env python3
"""
Figure 3D: CEACAM5 and CEACAM6 against the exhausted CD8+ T-cell (Tex,
C6_CD8_Tex_PDCD1) fraction, one point per stomach sample, coloured by
treatment phase and response. Spearman rho and P are printed in each box
and written to correlation_stats.csv (cnsfig.corr_stats).

Reads data/samples.csv beside this script; with --recompute, computes it from
the deposited epithelial h5ad (layers['counts'] normalised to log1p(CP10k))
and the CD8 T-cell h5ad (TCD8_H5AD).
Writes ceacam_tex_scatter.[svg|pdf|png] beside this script.

Run: python create_ceacam_tex_scatter.py [--recompute]
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
from paths import TCD8_H5AD, EPITHELIAL_DEPOSIT_H5AD
from cnsfig import layout as cnslayout, corr_stats, rich_xlabel, rich_ylabel  # noqa: E402
from cnsfig import cache, group_key, Panel, style  # noqa: E402

PANEL = Panel(3, "D", __file__, small_pt=4.5)
MARK, AREA = PANEL.mark, PANEL.area
SCALE = Panel.SCALE
# Left margin, shared with Figure 3E below so the frames align.
LEFT_MM = 12.5

EPI = EPITHELIAL_DEPOSIT_H5AD

# Shared four-group palette, plus grey for samples outside the four groups.
from shared.figure_config import GROUP_COLORS             # noqa: E402
COLOR_MAP = {**GROUP_COLORS, 'Other': '#999999'}

DRAW_ORDER = ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']


def load_primary():
    """Per stomach sample: mean epithelial CEACAM5/6 and CD8 Tex fraction (%)."""
    # --- Epithelial: CEACAM expression ---
    epi = sc.read_h5ad(EPI)
    epi.X = epi.layers['counts']
    del epi.layers['counts']
    # Drop uns['log1p'], which describes the replaced matrix; otherwise scanpy
    # treats the counts as already log-transformed.
    epi.uns.pop('log1p', None)

    stomach_mask = epi.obs['Sample site'].astype(str).str.lower().str.contains('stomach')
    epi = epi[stomach_mask].copy()

    sc.pp.normalize_total(epi, target_sum=1e4)
    sc.pp.log1p(epi)

    epi_df = pd.DataFrame({
        'sample': epi.obs['sample'].values,
        'CEACAM5': epi[:, 'CEACAM5'].X.toarray().flatten() if hasattr(epi[:, 'CEACAM5'].X, 'toarray') else epi[:, 'CEACAM5'].X.flatten(),
        'CEACAM6': epi[:, 'CEACAM6'].X.toarray().flatten() if hasattr(epi[:, 'CEACAM6'].X, 'toarray') else epi[:, 'CEACAM6'].X.flatten(),
        'treatment_phase': epi.obs['Treatment phase'].values,
        'pre_group': epi.obs['stomach_pre_grouping'].values,
        'post_group': epi.obs['stomach_post_grouping'].values,
    })
    epi_sample = epi_df.groupby('sample').agg({
        'CEACAM5': 'mean', 'CEACAM6': 'mean',
        'treatment_phase': 'first',
        'pre_group': 'first',
        'post_group': 'first',
    }).reset_index()
    stomach_samples = set(epi_sample['sample'].values)
    del epi

    # --- CD8: Tex fraction ---
    cd8 = sc.read_h5ad(TCD8_H5AD)
    cd8_df = pd.DataFrame({
        'sample': cd8.obs['sample'].values,
        'minor_cell_state': cd8.obs['minor_cell_state'].values,
    })
    del cd8

    cd8_df = cd8_df[cd8_df['sample'].isin(stomach_samples)]

    cd8_list = []
    for sample, grp in cd8_df.groupby('sample'):
        n_total = len(grp)
        n_tex = (grp['minor_cell_state'] == 'C6_CD8_Tex_PDCD1').sum()
        cd8_list.append({
            'sample': sample,
            'tex_fraction': (n_tex / n_total * 100) if n_total > 0 else 0,
        })
    cd8_sample = pd.DataFrame(cd8_list)

    merged = epi_sample.merge(cd8_sample, on='sample', how='inner')

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

    merged['group'] = merged.apply(assign_4group, axis=1)

    # The cached table carries no sample identifiers.
    return merged[['group', 'CEACAM5', 'CEACAM6', 'tex_fraction']].reset_index(drop=True)


def main():
    data = cache.table(PANEL.dir, 'samples', load_primary)

    # Two square boxes placed in mm (shared with Figures 3B, 3D and 3E).
    fig = PANEL.canvas()
    axes = cnslayout.scatter_pair_mm(fig, left_mm=LEFT_MM)

    rows = []
    for col_idx, ceacam in enumerate(['CEACAM5', 'CEACAM6']):
        ax = axes[col_idx]
        x = data[ceacam].values.astype(float)
        y = data['tex_fraction'].values.astype(float)
        groups = data['group'].values

        r_val, p_val = stats.spearmanr(x, y)

        # Marker about 0.52 mm across.
        for grp in DRAW_ORDER:
            mask = groups == grp
            if mask.sum() > 0:
                ax.scatter(x[mask], y[mask], c=COLOR_MAP[grp], s=19.4*AREA,
                           alpha=0.85, edgecolors='white', linewidths=style.EDGE_PT,
                           label=grp, zorder=3)

        # Regression line
        valid = np.isfinite(x) & np.isfinite(y)
        xv, yv = x[valid], y[valid]
        if len(xv) >= 3 and xv.std() > 0:
            slope, intercept = np.polyfit(xv, yv, 1)
            x_line = np.linspace(xv.min(), xv.max(), 100)
            ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=style.RULE_PT, alpha=0.6, zorder=2)

        rows.append((ceacam, r_val, p_val, int(valid.sum())))
        ax.set_title('In house', fontsize=style.tick_pt())
        cnslayout.corr_annotate(ax, r_val, p=p_val)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for spine in ['bottom', 'left']:
            ax.spines[spine].set_linewidth(style.RULE_PT)
        ax.tick_params(axis='both', width=style.RULE_PT, length=3*SCALE*MARK)
        rich_xlabel(ax, f"*{ceacam}*")
        # One y label, on the left box. A plain "+" rather than a superscript,
        # which would print below the minimum type size.
        if col_idx == 0:
            rich_ylabel(ax, 'CD8+ Tex\nfraction (%)')
    corr_stats.write(PANEL.dir, rows)

    # Group key to the right of the second box.
    group_key(fig, [(g, COLOR_MAP[g]) for g in DRAW_ORDER],
              x_mm=LEFT_MM + 2 * cnslayout.SCATTER_BOX_MM + cnslayout.SCATTER_GAP_MM + 0.8,
              y_mm=cnslayout.SCATTER_TOP_MM)

    PANEL.finish(fig, 'ceacam_tex_scatter')


if __name__ == '__main__':
    main()
