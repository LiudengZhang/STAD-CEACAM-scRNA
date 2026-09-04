#!/usr/bin/env python3
"""
Plot Epithelial Cell Type Delta Comparisons - CEACAM Split

Creates figure:
1. CEACAM-high vs CEACAM-low Delta: Epi CEACAM-high - Epi CEACAM-low

Author: Generated for Round 4 Analysis
Date: 2025-12-18
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Set plotting parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


def main():
    print("="*60)
    print("Plotting Epithelial Cell Type Deltas (CEACAM Split)")
    print("="*60)

    # Paths
    script_dir = Path(__file__).parent
    stats_dir = script_dir.parent.parent / '03_Extract_Results' / 'results' / 'summary_statistics'
    output_dir = script_dir.parent / 'results' / 'summary_plots'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    csv_path = stats_dir / 'mean_proportions_per_sample.csv'
    df = pd.read_csv(csv_path, index_col=0)
    print(f"Loaded: {csv_path}")
    print(f"Samples: {len(df)}")

    # Calculate delta: CEACAM-high - CEACAM-low
    df['Delta_CEACAM_high_vs_low'] = df['Epi CEACAM-high'] - df['Epi CEACAM-low']

    # Figure: CEACAM-high vs CEACAM-low Delta
    print("\nPlotting Figure: CEACAM-high vs CEACAM-low Delta...")
    fig, ax = plt.subplots(figsize=(12, 6))

    samples = df.index.tolist()
    deltas = df['Delta_CEACAM_high_vs_low'].values
    colors = ['#d62728' if d > 0 else '#1f77b4' for d in deltas]

    bars = ax.bar(samples, deltas, color=colors, edgecolor='black', linewidth=0.5)

    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_ylabel('Delta (CEACAM-high - CEACAM-low)', fontsize=12)
    ax.set_title('Epi CEACAM-high vs CEACAM-low Proportion Delta',
                fontsize=14, fontweight='bold')
    ax.set_xticklabels(samples, rotation=45, ha='right')

    # Add value labels
    for bar, val in zip(bars, deltas):
        ypos = val + 0.003 if val > 0 else val - 0.008
        ax.text(bar.get_x() + bar.get_width()/2, ypos, f'{val:.3f}',
                ha='center', va='bottom' if val > 0 else 'top', fontsize=9)

    plt.tight_layout()
    output_path = output_dir / 'ceacam_high_vs_low_delta.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")

    # Print summary
    print("\n" + "="*60)
    print("Summary Statistics")
    print("="*60)
    print("\nCEACAM-high vs CEACAM-low Delta:")
    print(f"  Mean: {df['Delta_CEACAM_high_vs_low'].mean():.4f}")
    print(f"  Range: [{df['Delta_CEACAM_high_vs_low'].min():.4f}, {df['Delta_CEACAM_high_vs_low'].max():.4f}]")
    print(f"  Positive (CEACAM-high > CEACAM-low): {(df['Delta_CEACAM_high_vs_low'] > 0).sum()}/{len(df)}")

    print("\n" + "="*60)
    print("Done!")
    print("="*60)


if __name__ == "__main__":
    main()
