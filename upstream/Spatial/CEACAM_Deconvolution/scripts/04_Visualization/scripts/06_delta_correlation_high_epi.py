#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
"""
Correlation Analysis by Sample - High Epithelial Spots Only

Filter to spots where Total Epithelial > 0.5 (50%), then analyze correlation.

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

# Filter threshold
EPI_THRESHOLD = 0.5


def main():
    print("="*70)
    print("Correlation Analysis - High Epithelial Spots Only (>50%)")
    print("="*70)

    # Paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'summary_plots'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load GraphST results
    results_dir = Path('/path/to/home/Project_4_Gastric_Cancer/01.2_GraphST_17Types/03_Results')
    h5ad_files = sorted(results_dir.glob('*_graphst_output.h5ad'))

    print(f"Loading {len(h5ad_files)} samples...")
    print(f"Filter: Total Epithelial > {EPI_THRESHOLD}")

    # Create figure: 2 rows x 5 columns
    fig, axes = plt.subplots(2, 5, figsize=(25, 10))
    axes = axes.flatten()

    # Store results
    results = []

    for idx, h5ad_path in enumerate(h5ad_files):
        sample_name = h5ad_path.stem.replace('_graphst_output', '')
        adata = sc.read_h5ad(h5ad_path)

        # Get proportions
        f10_high = adata.obs['Cancer F10-high'].values
        f10_low = adata.obs['Cancer F10-low'].values
        normal_epi = adata.obs['Normal Epi'].values

        # Calculate total epithelial
        total_epi = normal_epi + f10_high + f10_low

        # Filter: keep spots with high epithelial
        mask = total_epi > EPI_THRESHOLD
        n_total = len(mask)
        n_filtered = mask.sum()

        if n_filtered < 10:
            print(f"  {sample_name}: {n_filtered}/{n_total} spots pass filter - SKIPPING (too few)")
            results.append({
                'sample': sample_name,
                'n_total': n_total,
                'n_filtered': n_filtered,
                'r': np.nan,
                'p': np.nan
            })
            axes[idx].text(0.5, 0.5, f'{sample_name}\nToo few spots\n({n_filtered}/{n_total})',
                          ha='center', va='center', fontsize=12, transform=axes[idx].transAxes)
            axes[idx].set_xticks([])
            axes[idx].set_yticks([])
            continue

        # Apply filter
        f10_high_f = f10_high[mask]
        f10_low_f = f10_low[mask]
        normal_epi_f = normal_epi[mask]

        # Calculate deltas
        tumor_avg = (f10_high_f + f10_low_f) / 2
        x = tumor_avg - normal_epi_f  # Delta Tumor vs Normal
        y = f10_high_f - f10_low_f     # Delta F10-high vs F10-low

        # Regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

        results.append({
            'sample': sample_name,
            'n_total': n_total,
            'n_filtered': n_filtered,
            'r': r_value,
            'p': p_value
        })

        # Plot
        ax = axes[idx]
        ax.scatter(x, y, alpha=0.4, s=10, c='steelblue', edgecolors='none')

        # Regression line
        x_line = np.linspace(x.min(), x.max(), 100)
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, 'r-', linewidth=2)

        # Zero lines
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

        # Title with stats
        p_str = f'{p_value:.1e}' if p_value < 0.001 else f'{p_value:.3f}'
        ax.set_title(f'{sample_name}\nr = {r_value:.3f}, n = {n_filtered}/{n_total}', fontsize=11)

        ax.set_xlabel('Tumor - Normal', fontsize=9)
        ax.set_ylabel('F10-high - F10-low', fontsize=9)
        ax.tick_params(labelsize=8)

        # Color title based on correlation
        if r_value > 0.3:
            ax.title.set_color('darkred')
        elif r_value < -0.3:
            ax.title.set_color('darkblue')

        print(f"  {sample_name}: r = {r_value:.3f}, n = {n_filtered}/{n_total}")

    fig.suptitle(f'Delta Correlation by Sample (Epithelial >50% spots only)\n(Tumor-Normal vs F10high-F10low)',
                 fontsize=16, fontweight='bold', y=1.02)

    plt.tight_layout()

    output_path = output_dir / 'delta_correlation_by_sample_high_epi.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    print(f"\nSaved: {output_path}")

    # Summary
    print("\n" + "="*70)
    print("Summary: Correlation by Sample (High Epithelial Spots)")
    print("="*70)
    print(f"\n{'Sample':<12} {'n_filt':>8} {'n_total':>8} {'%kept':>8} {'r':>8}")
    print("-"*50)

    valid_results = [r for r in results if not np.isnan(r['r'])]
    for res in results:
        if np.isnan(res['r']):
            print(f"{res['sample']:<12} {res['n_filtered']:>8} {res['n_total']:>8} {res['n_filtered']/res['n_total']*100:>7.1f}% {'N/A':>8}")
        else:
            print(f"{res['sample']:<12} {res['n_filtered']:>8} {res['n_total']:>8} {res['n_filtered']/res['n_total']*100:>7.1f}% {res['r']:>8.3f}")

    # Overall pattern
    if valid_results:
        pos_count = sum(1 for r in valid_results if r['r'] > 0)
        neg_count = sum(1 for r in valid_results if r['r'] < 0)
        strong_pos = sum(1 for r in valid_results if r['r'] > 0.3)
        strong_neg = sum(1 for r in valid_results if r['r'] < -0.3)

        print("\n" + "="*70)
        print("Overall Pattern (filtered spots)")
        print("="*70)
        print(f"Positive correlation: {pos_count}/{len(valid_results)} samples ({strong_pos} strong)")
        print(f"Negative correlation: {neg_count}/{len(valid_results)} samples ({strong_neg} strong)")

    print("\n" + "="*70)
    print("Done!")
    print("="*70)


if __name__ == "__main__":
    main()
