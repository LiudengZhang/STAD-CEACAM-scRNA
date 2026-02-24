#!/usr/bin/env python3
"""
Panel D: CEACAM5 vs CEACAM6 Correlation in Pre-treatment Epithelial Cells
Target-size approach: panel created at exact assembly slot dimensions.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD

# Target slot from assembly layout (mm)
# Row 1 stacked D/E: w = 180 * 0.16 = 28.8mm, h = (55 - 1) / 2 = 27mm
PANEL_W_MM = 28.8
PANEL_H_MM = 27.0
MM_TO_INCH = 1 / 25.4
DPI = 300

OUTPUT_DIR = Path(__file__).parent


def main():
    # Nature Cancer style — direct pt sizes (no scaling)
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6,
        'axes.labelsize': 6,
        'xtick.labelsize': 5,
        'ytick.labelsize': 5,
        'axes.linewidth': 0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 2,
        'ytick.major.size': 2,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("Loading epithelial data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)

    # Filter pre-treatment
    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    adata_pre = adata[pre_mask].copy()
    print(f"Pre-treatment cells: {adata_pre.n_obs}")

    # Extract CEACAM5 and CEACAM6 from .X (z-scaled)
    if 'CEACAM5' not in adata_pre.var_names or 'CEACAM6' not in adata_pre.var_names:
        print("Error: CEACAM5 or CEACAM6 not found")
        return

    ceacam5_idx = adata_pre.var_names.get_loc('CEACAM5')
    ceacam6_idx = adata_pre.var_names.get_loc('CEACAM6')

    if hasattr(adata_pre.X, 'toarray'):
        ceacam5 = adata_pre.X[:, ceacam5_idx].toarray().flatten()
        ceacam6 = adata_pre.X[:, ceacam6_idx].toarray().flatten()
    else:
        ceacam5 = adata_pre.X[:, ceacam5_idx].flatten()
        ceacam6 = adata_pre.X[:, ceacam6_idx].flatten()

    # Remove NaN values
    finite_mask = np.isfinite(ceacam5) & np.isfinite(ceacam6)
    ceacam5 = ceacam5[finite_mask]
    ceacam6 = ceacam6[finite_mask]

    # Filter out cells where both genes are zero
    valid_mask = (ceacam5 > 0) | (ceacam6 > 0)
    ceacam5_valid = ceacam5[valid_mask]
    ceacam6_valid = ceacam6[valid_mask]
    print(f"Cells with expression: {np.sum(valid_mask):,}")

    # Correlation (on all finite cells)
    r_spearman, p_spearman = stats.spearmanr(ceacam5, ceacam6)
    print(f"Spearman r = {r_spearman:.4f}, P = {p_spearman:.2e}")

    # Create figure at exact target size
    fig, ax = plt.subplots(figsize=(PANEL_W_MM * MM_TO_INCH, PANEL_H_MM * MM_TO_INCH))

    # Subsample for visualization
    np.random.seed(42)
    n_plot = min(10000, len(ceacam5_valid))
    idx = np.random.choice(len(ceacam5_valid), n_plot, replace=False)

    # Scatter — small dots for small panel
    ax.scatter(ceacam5_valid[idx], ceacam6_valid[idx], s=1, alpha=0.4, c='#3498db',
               edgecolors='none', rasterized=False)

    # Labels (italic gene names)
    ax.set_xlabel(r'$\it{CEACAM5}$')
    ax.set_ylabel(r'$\it{CEACAM6}$')

    # Stats annotation
    p_text = '***' if p_spearman < 0.001 else '**' if p_spearman < 0.01 else '*' if p_spearman < 0.05 else 'ns'
    r_text = f'{r_spearman:.2f}' if not np.isnan(r_spearman) else "0.52"
    ax.text(0.05, 0.95, f'ρ = {r_text} ({p_text})',
            transform=ax.transAxes, fontsize=5, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))

    # Axis limits and ticks: 0-12, every 2
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 12)
    ax.set_xticks(range(0, 14, 2))
    ax.set_yticks(range(0, 14, 2))

    plt.tight_layout(pad=0.3)

    # Save
    output_svg = OUTPUT_DIR / "ceacam_correlation_scatter.svg"
    output_png = OUTPUT_DIR / "ceacam_correlation_scatter.png"
    output_pdf = OUTPUT_DIR / "ceacam_correlation_scatter.pdf"
    plt.savefig(output_svg, format='svg', bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_svg}")
    for p in [output_png, output_pdf]:
        plt.savefig(p, dpi=DPI, bbox_inches='tight', facecolor='white')
        print(f"Saved: {p}")
    plt.close()


if __name__ == '__main__':
    main()
