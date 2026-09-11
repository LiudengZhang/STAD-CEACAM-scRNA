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
Tumor Score Histogram per Sample

Creates histogram of tumor_score = F10_high + F10_low for each sample.
This helps determine threshold for filtering tumor-enriched spots.

Author: Generated for Round 4 Analysis
Date: 2025-11-26
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


def load_graphst_results(results_dir):
    """Load all GraphST results."""
    print("Loading GraphST results...")
    h5ad_files = list(Path(results_dir).glob('*_graphst_output.h5ad'))

    adata_dict = {}
    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem.replace('_graphst_output', '')
        adata = sc.read_h5ad(h5ad_path)
        adata_dict[sample_name] = adata
        print(f"  {sample_name}: {adata.n_obs} spots")

    return adata_dict


def main():
    print("="*70)
    print("Tumor Score Histogram per Sample")
    print("tumor_score = F10_high + F10_low")
    print("="*70)

    # Paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'summary_plots'
    output_dir.mkdir(parents=True, exist_ok=True)

    results_dir = Path('/path/to/home/Project_4_Gastric_Cancer/01.2_GraphST_17Types/03_Results')
    adata_dict = load_graphst_results(results_dir)

    # Create 2x5 grid for 10 samples
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    axes = axes.flatten()

    samples = sorted(adata_dict.keys())

    print("\nTumor Score Statistics:")
    print(f"{'Sample':<12} {'Mean':<10} {'Median':<10} {'75th%':<10} {'90th%':<10} {'Max':<10}")
    print("-"*65)

    for idx, sample in enumerate(samples):
        adata = adata_dict[sample]
        tumor_score = adata.obs['Cancer F10-high'].values + adata.obs['Cancer F10-low'].values

        ax = axes[idx]
        ax.hist(tumor_score, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
        ax.axvline(np.median(tumor_score), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(tumor_score):.3f}')
        ax.axvline(np.percentile(tumor_score, 75), color='orange', linestyle='--', linewidth=2, label=f'75th: {np.percentile(tumor_score, 75):.3f}')

        ax.set_xlabel('Tumor Score (F10_high + F10_low)', fontsize=9)
        ax.set_ylabel('Count', fontsize=9)
        ax.set_title(f'{sample}\n(n={len(tumor_score):,})', fontsize=10, fontweight='bold')
        ax.legend(fontsize=7, loc='upper right')

        # Print statistics
        print(f"{sample:<12} {tumor_score.mean():<10.4f} {np.median(tumor_score):<10.4f} "
              f"{np.percentile(tumor_score, 75):<10.4f} {np.percentile(tumor_score, 90):<10.4f} "
              f"{tumor_score.max():<10.4f}")

    plt.suptitle('Tumor Score Distribution per Sample\n(F10_high + F10_low)', fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_path = output_dir / 'tumor_score_histogram_per_sample.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {output_path}")

    print("\n" + "="*70)
    print("Done!")
    print("="*70)


if __name__ == "__main__":
    main()
