#!/usr/bin/env python3
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
"""
Correlation Analysis: Tumor vs Normal Delta vs F10-high vs F10-low Delta

Creates scatterplot with regression to assess correlation between:
1. Delta (Tumor - Normal Epi)
2. Delta (F10-high - F10-low)

Author: Generated for Round 4 Analysis
Date: 2025-11-26
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


def main():
    print("="*70)
    print("Correlation Analysis: Delta Scores")
    print("="*70)

    # Paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'summary_plots'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load GraphST results
    results_dir = Path('/path/to/home/Project_4_Gastric_Cancer/01.2_GraphST_17Types/03_Results')
    h5ad_files = sorted(results_dir.glob('*_graphst_output.h5ad'))

    print(f"Loading {len(h5ad_files)} samples...")

    # Collect all deltas
    all_delta_tumor_normal = []
    all_delta_f10 = []
    all_samples = []

    for h5ad_path in h5ad_files:
        sample_name = h5ad_path.stem.replace('_graphst_output', '')
        adata = sc.read_h5ad(h5ad_path)

        f10_high = adata.obs['Cancer F10-high'].values
        f10_low = adata.obs['Cancer F10-low'].values
        normal_epi = adata.obs['Normal Epi'].values

        tumor_avg = (f10_high + f10_low) / 2
        delta_tumor_normal = tumor_avg - normal_epi
        delta_f10 = f10_high - f10_low

        all_delta_tumor_normal.extend(delta_tumor_normal)
        all_delta_f10.extend(delta_f10)
        all_samples.extend([sample_name] * len(delta_tumor_normal))

        print(f"  {sample_name}: {len(delta_tumor_normal)} spots")

    # Convert to arrays
    x = np.array(all_delta_tumor_normal)
    y = np.array(all_delta_f10)

    print(f"\nTotal spots: {len(x)}")

    # Regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r_squared = r_value ** 2

    # Spearman correlation (more robust)
    spearman_r, spearman_p = stats.spearmanr(x, y)

    print("\n" + "="*70)
    print("Correlation Results")
    print("="*70)
    print(f"\nPearson correlation:")
    print(f"  r = {r_value:.4f}")
    print(f"  R² = {r_squared:.4f}")
    print(f"  p-value = {p_value:.2e}")

    print(f"\nSpearman correlation:")
    print(f"  rho = {spearman_r:.4f}")
    print(f"  p-value = {spearman_p:.2e}")

    print(f"\nLinear regression:")
    print(f"  y = {slope:.4f}x + {intercept:.4f}")

    # Create scatterplot
    fig, ax = plt.subplots(figsize=(10, 8))

    # Scatter with alpha for density
    ax.scatter(x, y, alpha=0.3, s=5, c='steelblue', edgecolors='none')

    # Regression line
    x_line = np.linspace(x.min(), x.max(), 100)
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, 'r-', linewidth=2, label=f'Regression line')

    # Zero lines
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)

    # Labels
    ax.set_xlabel('Delta (Tumor - Normal Epi)', fontsize=12)
    ax.set_ylabel('Delta (F10-high - F10-low)', fontsize=12)
    ax.set_title('Correlation Between Epithelial Delta Scores\n(Per Spot, All Samples)',
                fontsize=14, fontweight='bold')

    # Stats annotation
    stats_text = (f'Pearson r = {r_value:.3f}\n'
                  f'R² = {r_squared:.3f}\n'
                  f'p < {p_value:.1e}\n'
                  f'n = {len(x):,} spots')
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Correlation interpretation
    if r_value > 0:
        corr_type = "POSITIVE"
    else:
        corr_type = "NEGATIVE"

    ax.text(0.95, 0.05, f'{corr_type} correlation', transform=ax.transAxes, fontsize=12,
            verticalalignment='bottom', horizontalalignment='right',
            color='darkred' if r_value > 0 else 'darkblue', fontweight='bold')

    plt.tight_layout()

    output_path = output_dir / 'delta_correlation_scatterplot.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    print(f"\nSaved: {output_path}")

    # Interpretation
    print("\n" + "="*70)
    print("Interpretation")
    print("="*70)
    if r_value > 0.3:
        print(f"\nThe two delta scores are POSITIVELY correlated (r = {r_value:.3f}).")
        print("Spots with higher Tumor vs Normal delta tend to also have higher F10-high vs F10-low delta.")
    elif r_value < -0.3:
        print(f"\nThe two delta scores are NEGATIVELY correlated (r = {r_value:.3f}).")
        print("Spots with higher Tumor vs Normal delta tend to have lower F10-high vs F10-low delta.")
    else:
        print(f"\nThe correlation is weak (r = {r_value:.3f}).")
        print("The two delta scores show little linear relationship.")

    print("\n" + "="*70)
    print("Done!")
    print("="*70)


if __name__ == "__main__":
    main()
