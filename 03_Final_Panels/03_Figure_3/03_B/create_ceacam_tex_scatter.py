#!/usr/bin/env python3
"""
Figure R1-B: CEACAM5/6 vs CD8+ Tex fraction (1x2)
ALL 32 stomach samples, 4-group coloring, no external dataset.
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
from paths import TCD8_H5AD, EPITHELIAL_RAW_COUNTS_H5AD

# 4x scaling
SCALE = 4
DPI = 300
PANEL_WIDTH_CM = 8.0 * SCALE   # wide enough for 1x2
PANEL_HEIGHT_CM = 4.0 * SCALE  # half height (single row)
CM_TO_INCH = 1 / 2.54

EPI_RAW = EPITHELIAL_RAW_COUNTS_H5AD
OUTPUT_DIR = Path(__file__).parent

# Standard 4-group palette (blue=R, warm=NR; lighter=Pre, darker=Post)
COLOR_MAP = {
    'Pre-R':  '#74b9ff',
    'Pre-NR': '#e17055',
    'Post-R': '#0984e3',
    'Post-NR':'#d63031',
    'Other':  '#999999',
}

DRAW_ORDER = ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']


def load_primary():
    """Load ALL 32 stomach samples: epithelial CEACAM + CD8 Tex fraction."""
    # --- Epithelial: CEACAM expression ---
    print("[Primary] Loading epithelial raw counts...")
    epi = sc.read_h5ad(EPI_RAW)

    # Filter: stomach only (all treatment phases)
    stomach_mask = epi.obs['Sample site'].astype(str).str.lower().str.contains('stomach')
    epi = epi[stomach_mask].copy()

    # Normalize raw counts
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
    print(f"  {len(epi_sample)} stomach samples")
    del epi

    # --- CD8: Tex fraction ---
    print("[Primary] Loading CD8+ T cells...")
    cd8 = sc.read_h5ad(TCD8_H5AD)
    cd8_df = pd.DataFrame({
        'sample': cd8.obs['sample'].values,
        'minor_cell_state': cd8.obs['minor_cell_state'].values,
    })
    del cd8

    # Filter to stomach samples
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
    print(f"  CD8 data for {len(cd8_sample)} samples")

    # Merge
    merged = epi_sample.merge(cd8_sample, on='sample', how='inner')

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

    merged['group'] = merged.apply(assign_4group, axis=1)

    counts = merged['group'].value_counts()
    print(f"  Merged: {len(merged)} samples: {counts.to_dict()}")
    return merged


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    data = load_primary()

    # 1x2 layout
    fig_w = PANEL_WIDTH_CM * CM_TO_INCH
    fig_h = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h))

    for col_idx, ceacam in enumerate(['CEACAM5', 'CEACAM6']):
        ax = axes[col_idx]
        x = data[ceacam].values.astype(float)
        y = data['tex_fraction'].values.astype(float)
        groups = data['group'].values

        r_val, p_val = stats.spearmanr(x, y)

        for grp in DRAW_ORDER:
            mask = groups == grp
            if mask.sum() > 0:
                ax.scatter(x[mask], y[mask], c=COLOR_MAP[grp], s=30*SCALE,
                           alpha=0.85, edgecolors='white', linewidths=0.3*SCALE,
                           label=grp, zorder=3)

        # Regression line
        valid = np.isfinite(x) & np.isfinite(y)
        xv, yv = x[valid], y[valid]
        if len(xv) >= 3 and xv.std() > 0:
            slope, intercept = np.polyfit(xv, yv, 1)
            x_line = np.linspace(xv.min(), xv.max(), 100)
            ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=0.8*SCALE, alpha=0.6, zorder=2)

        # Stats text — asterisk notation
        p_str = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'

        ax.set_title(f'Primary Cohort (scRNA-seq)\nρ = {r_val:.2f}, {p_str}',
                     fontsize=6.5 * SCALE, fontweight='normal', linespacing=1.4)

        ax.set_xlabel(f'{ceacam} (mean log expr.)', fontsize=6 * SCALE)
        ax.set_ylabel('CD8$^+$ Tex fraction (%)', fontsize=6 * SCALE)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for spine in ['bottom', 'left']:
            ax.spines[spine].set_linewidth(0.5*SCALE)
        ax.tick_params(axis='both', labelsize=5 * SCALE, width=0.5*SCALE, length=3*SCALE)

        ax.set_box_aspect(1)

    # Single shared legend on the far right
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0, 0.5),
               fontsize=4.5 * SCALE, framealpha=0, edgecolor='none', markerscale=0.8)

    plt.tight_layout()

    out = OUTPUT_DIR / 'ceacam_tex_scatter.png'
    fig.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(str(out).replace('.png', '.pdf'), dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(str(out).replace('.png', '.svg'), format='svg', bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {out}")
    plt.close()


if __name__ == '__main__':
    main()
