#!/usr/bin/env python3
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
"""
Correlation Analysis by Sample: Delta Scores

Creates 10-panel figure (2x5 grid) showing per-sample correlation between:
- X: Delta (Tumor - Normal Epi)
- Y: Delta (F10-high - F10-low)

Author: Generated for Round 4 Analysis
Date: 2025-11-26
"""

import numpy as np
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
    print("Correlation Analysis by Sample: Delta Scores")
    print("="*70)

    # Paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'summary_plots'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load GraphST results
    results_dir = Path('/path/to/home/Project_4_Gastric_Cancer/01.2_GraphST_17Types/03_Results')
    h5ad_files = sorted(results_dir.glob('*_graphst_output.h5ad'))

    print(f"Loading {len(h5ad_files)} samples...")

    # Create figure: 2 rows x 5 columns
    fig, axes = plt.subplots(2, 5, figsize=(25, 10))
    axes = axes.flatten()

    # Store results for summary
    results = []

    for idx, h5ad_path in enumerate(h5ad_files):
        sample_name = h5ad_path.stem.replace('_graphst_output', '')
        adata = sc.read_h5ad(h5ad_path)

        # Calculate deltas
        f10_high = adata.obs['Cancer F10-high'].values
        f10_low = adata.obs['Cancer F10-low'].values
        normal_epi = adata.obs['Normal Epi'].values

        tumor_avg = (f10_high + f10_low) / 2
        x = tumor_avg - normal_epi  # Delta Tumor vs Normal
        y = f10_high - f10_low       # Delta F10-high vs F10-low

        # Regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

        results.append({
            'sample': sample_name,
            'n': len(x),
            'r': r_value,
            'p': p_value,
            'slope': slope
        })

        # Plot
        ax = axes[idx]
        ax.scatter(x, y, alpha=0.4, s=8, c='steelblue', edgecolors='none')

        # Regression line
        x_line = np.linspace(x.min(), x.max(), 100)
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, 'r-', linewidth=2)

        # Zero lines
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

        # Title with stats
        p_str = f'{p_value:.1e}' if p_value < 0.001 else f'{p_value:.3f}'
        ax.set_title(f'{sample_name}\nr = {r_value:.3f}, p = {p_str}', fontsize=11)

        ax.set_xlabel('Tumor - Normal', fontsize=9)
        ax.set_ylabel('F10-high - F10-low', fontsize=9)
        ax.tick_params(labelsize=8)

        # Color title based on correlation direction
        if r_value > 0.3:
            ax.title.set_color('darkred')
        elif r_value < -0.3:
            ax.title.set_color('darkblue')

        print(f"  {sample_name}: r = {r_value:.3f}, p = {p_str}, n = {len(x)}")

    fig.suptitle('Delta Correlation by Sample\n(Tumor-Normal vs F10high-F10low)',
                 fontsize=16, fontweight='bold', y=1.02)

    plt.tight_layout()

    output_path = output_dir / 'delta_correlation_by_sample.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    print(f"\nSaved: {output_path}")

    # Summary table
    print("\n" + "="*70)
    print("Summary: Correlation by Sample")
    print("="*70)
    print(f"\n{'Sample':<12} {'n':>6} {'r':>8} {'p-value':>12} {'Direction':<12}")
    print("-"*50)

    for res in results:
        direction = "Positive" if res['r'] > 0 else "Negative"
        strength = "strong" if abs(res['r']) > 0.3 else "weak"
        p_str = f"{res['p']:.1e}" if res['p'] < 0.001 else f"{res['p']:.4f}"
        print(f"{res['sample']:<12} {res['n']:>6} {res['r']:>8.3f} {p_str:>12} {direction} ({strength})")

    # Overall pattern
    pos_count = sum(1 for r in results if r['r'] > 0)
    neg_count = sum(1 for r in results if r['r'] < 0)
    strong_pos = sum(1 for r in results if r['r'] > 0.3)
    strong_neg = sum(1 for r in results if r['r'] < -0.3)

    print("\n" + "="*70)
    print("Overall Pattern")
    print("="*70)
    print(f"Positive correlation: {pos_count}/10 samples ({strong_pos} strong)")
    print(f"Negative correlation: {neg_count}/10 samples ({strong_neg} strong)")

    print("\n" + "="*70)
    print("Done!")
    print("="*70)


if __name__ == "__main__":
    main()
