#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
CEACAM Tumor-Centered Analysis

Tests whether CEACAM-high tumor cells are more "tumor-centered" (in tumor core)
vs CEACAM-low cells being at the tumor-stroma interface.

Three complementary analyses:
1. Neighborhood Epithelial Density: k-NN mean epithelial proportion
2. Distance to Stroma: Distance to nearest stroma-dominant spot
3. Distance to Immune Cells: Mean distance to immune-rich spots

Author: Generated for Round 4 Analysis
Date: 2026-01-04
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree
from pathlib import Path
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'

# Parameters
K_NEIGHBORS = 15  # k for k-NN neighborhood analysis
EPI_THRESHOLD = 0.05  # Minimum epithelial to be "tumor spot"
STROMA_THRESHOLD = 0.20  # Maximum epithelial to be "stroma-dominant"
IMMUNE_THRESHOLD = 0.15  # Minimum total immune to be "immune-rich"


def load_graphst_results(results_dir):
    """Load all GraphST results with spatial coordinates."""
    print("Loading GraphST results...")
    h5ad_files = list(Path(results_dir).glob('*_graphst_output.h5ad'))

    all_data = []
    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem.replace('_graphst_output', '')
        adata = sc.read_h5ad(h5ad_path)

        # Extract data
        df = adata.obs.copy()
        df['sample'] = sample_name
        df['x'] = adata.obsm['spatial'][:, 0]
        df['y'] = adata.obsm['spatial'][:, 1]
        df['spot_id'] = df.index

        all_data.append(df)
        print(f"  {sample_name}: {len(df)} spots")

    df_all = pd.concat(all_data, ignore_index=True)
    print(f"Total spots: {len(df_all):,}")
    return df_all


def calculate_neighborhood_epithelial_density(df, k=K_NEIGHBORS):
    """
    For each spot, calculate mean epithelial proportion in k-nearest neighbors.
    """
    print(f"\n[1/3] Calculating neighborhood epithelial density (k={k})...")

    results = []

    for sample in tqdm(df['sample'].unique(), desc="Processing samples"):
        sample_df = df[df['sample'] == sample].copy()

        # Build KD-tree for this sample
        coords = sample_df[['x', 'y']].values
        tree = cKDTree(coords)

        # For each spot, find k+1 nearest neighbors (including self)
        distances, indices = tree.query(coords, k=k+1)

        # Calculate neighborhood epithelial density (excluding self)
        total_epi = sample_df['Epi CEACAM-high'].values + sample_df['Epi CEACAM-low'].values

        for i in range(len(sample_df)):
            neighbor_idx = indices[i, 1:]  # Exclude self
            neighbor_epi = total_epi[neighbor_idx].mean()

            results.append({
                'sample': sample,
                'spot_id': sample_df.iloc[i]['spot_id'],
                'neighborhood_epi_density': neighbor_epi
            })

    return pd.DataFrame(results)


def calculate_distance_to_stroma(df, stroma_threshold=STROMA_THRESHOLD):
    """
    For each tumor spot, calculate distance to nearest stroma-dominant spot.
    Stroma-dominant = Total epithelial < threshold
    """
    print(f"\n[2/3] Calculating distance to stroma (threshold={stroma_threshold})...")

    results = []

    for sample in tqdm(df['sample'].unique(), desc="Processing samples"):
        sample_df = df[df['sample'] == sample].copy()

        # Calculate total epithelial
        sample_df['Total_Epi'] = sample_df['Epi CEACAM-high'] + sample_df['Epi CEACAM-low']

        # Identify stroma-dominant spots
        stroma_mask = sample_df['Total_Epi'] < stroma_threshold
        stroma_spots = sample_df[stroma_mask]

        if len(stroma_spots) < 5:
            print(f"  Warning: {sample} has only {len(stroma_spots)} stroma spots, skipping")
            continue

        # Build KD-tree for stroma spots
        stroma_coords = stroma_spots[['x', 'y']].values
        stroma_tree = cKDTree(stroma_coords)

        # For all spots, calculate distance to nearest stroma
        all_coords = sample_df[['x', 'y']].values
        distances, _ = stroma_tree.query(all_coords, k=1)

        for i in range(len(sample_df)):
            results.append({
                'sample': sample,
                'spot_id': sample_df.iloc[i]['spot_id'],
                'distance_to_stroma': distances[i]
            })

    return pd.DataFrame(results)


def calculate_distance_to_immune(df, immune_threshold=IMMUNE_THRESHOLD):
    """
    For each tumor spot, calculate mean distance to immune-rich spots.
    Immune-rich = Total immune proportion > threshold
    """
    print(f"\n[3/3] Calculating distance to immune-rich spots (threshold={immune_threshold})...")

    immune_cols = ['B cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells',
                   'DC cells', 'Mast cells', 'Neutrophils', 'Plasma cells',
                   'Monocytes/Macrophages']

    results = []

    for sample in tqdm(df['sample'].unique(), desc="Processing samples"):
        sample_df = df[df['sample'] == sample].copy()

        # Calculate total immune
        sample_df['Total_Immune'] = sample_df[immune_cols].sum(axis=1)

        # Identify immune-rich spots
        immune_mask = sample_df['Total_Immune'] > immune_threshold
        immune_spots = sample_df[immune_mask]

        if len(immune_spots) < 5:
            print(f"  Warning: {sample} has only {len(immune_spots)} immune-rich spots, skipping")
            continue

        # Build KD-tree for immune spots
        immune_coords = immune_spots[['x', 'y']].values
        immune_tree = cKDTree(immune_coords)

        # For all spots, calculate distance to nearest 5 immune spots (mean)
        all_coords = sample_df[['x', 'y']].values
        k_immune = min(5, len(immune_spots))
        distances, _ = immune_tree.query(all_coords, k=k_immune)
        mean_distances = distances.mean(axis=1) if k_immune > 1 else distances

        for i in range(len(sample_df)):
            results.append({
                'sample': sample,
                'spot_id': sample_df.iloc[i]['spot_id'],
                'distance_to_immune': mean_distances[i]
            })

    return pd.DataFrame(results)


def merge_and_analyze(df, neigh_df, stroma_df, immune_df, epi_threshold=EPI_THRESHOLD):
    """
    Merge all metrics and compare CEACAM-high vs CEACAM-low tumor spots.
    """
    print("\n" + "="*70)
    print("Merging and analyzing results...")
    print("="*70)

    # Calculate CEACAM ratio for original data
    df['Total_Epi'] = df['Epi CEACAM-high'] + df['Epi CEACAM-low']
    df_tumor = df[df['Total_Epi'] >= epi_threshold].copy()
    df_tumor['CEACAM_ratio'] = df_tumor['Epi CEACAM-high'] / df_tumor['Total_Epi']

    print(f"Tumor spots (>={epi_threshold*100:.0f}% epithelial): {len(df_tumor):,}")

    # Merge metrics
    df_merged = df_tumor.merge(neigh_df, on=['sample', 'spot_id'], how='left')
    df_merged = df_merged.merge(stroma_df, on=['sample', 'spot_id'], how='left')
    df_merged = df_merged.merge(immune_df, on=['sample', 'spot_id'], how='left')

    # Drop rows with missing values
    df_merged = df_merged.dropna(subset=['neighborhood_epi_density', 'distance_to_stroma', 'distance_to_immune'])
    print(f"Spots with all metrics: {len(df_merged):,}")

    # Classify by CEACAM (median split)
    median_ratio = df_merged['CEACAM_ratio'].median()
    df_merged['CEACAM_group'] = np.where(
        df_merged['CEACAM_ratio'] >= median_ratio,
        'CEACAM-high',
        'CEACAM-low'
    )

    print(f"CEACAM-high spots: {(df_merged['CEACAM_group'] == 'CEACAM-high').sum():,}")
    print(f"CEACAM-low spots: {(df_merged['CEACAM_group'] == 'CEACAM-low').sum():,}")
    print(f"Median CEACAM_ratio threshold: {median_ratio:.3f}")

    return df_merged


def statistical_comparison(df_merged):
    """
    Compare metrics between CEACAM-high and CEACAM-low groups.
    """
    print("\n" + "="*70)
    print("Statistical Comparison")
    print("="*70)

    ceacam_high = df_merged[df_merged['CEACAM_group'] == 'CEACAM-high']
    ceacam_low = df_merged[df_merged['CEACAM_group'] == 'CEACAM-low']

    metrics = [
        ('neighborhood_epi_density', 'Neighborhood Epithelial Density', 'HIGHER = more tumor-centered'),
        ('distance_to_stroma', 'Distance to Stroma', 'HIGHER = more tumor-centered'),
        ('distance_to_immune', 'Distance to Immune-rich Spots', 'HIGHER = more tumor-centered')
    ]

    results = []

    for col, name, interpretation in metrics:
        high_vals = ceacam_high[col].values
        low_vals = ceacam_low[col].values

        # Mann-Whitney U test
        stat, pval = stats.mannwhitneyu(high_vals, low_vals, alternative='two-sided')

        # Effect size (Cohen's d approximation)
        pooled_std = np.sqrt((high_vals.std()**2 + low_vals.std()**2) / 2)
        cohens_d = (high_vals.mean() - low_vals.mean()) / pooled_std if pooled_std > 0 else 0

        # Direction
        diff = high_vals.mean() - low_vals.mean()
        direction = 'CEACAM-high > CEACAM-low' if diff > 0 else 'CEACAM-high < CEACAM-low'

        results.append({
            'metric': name,
            'CEACAM_high_mean': high_vals.mean(),
            'CEACAM_low_mean': low_vals.mean(),
            'difference': diff,
            'cohens_d': cohens_d,
            'pvalue': pval,
            'direction': direction,
            'interpretation': interpretation
        })

        sig = '***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else ''))
        print(f"\n{name}:")
        print(f"  CEACAM-high: {high_vals.mean():.4f} ± {high_vals.std():.4f}")
        print(f"  CEACAM-low:  {low_vals.mean():.4f} ± {low_vals.std():.4f}")
        print(f"  Difference:  {diff:+.4f}")
        print(f"  Cohen's d:   {cohens_d:.3f}")
        print(f"  p-value:     {pval:.2e} {sig}")
        print(f"  Direction:   {direction}")
        print(f"  ({interpretation})")

    return pd.DataFrame(results)


def create_visualization(df_merged, results_df, output_dir):
    """
    Create comprehensive visualization of all three analyses.
    """
    print("\nCreating visualizations...")

    ceacam_high = df_merged[df_merged['CEACAM_group'] == 'CEACAM-high']
    ceacam_low = df_merged[df_merged['CEACAM_group'] == 'CEACAM-low']

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    metrics = [
        ('neighborhood_epi_density', 'Neighborhood Epithelial Density', 'Higher = More Tumor-Centered'),
        ('distance_to_stroma', 'Distance to Stroma', 'Higher = More Tumor-Centered'),
        ('distance_to_immune', 'Distance to Immune-rich Spots', 'Higher = More Tumor-Centered')
    ]

    colors = {'CEACAM-high': '#d62728', 'CEACAM-low': '#1f77b4'}

    # Row 1: Boxplots
    for idx, (col, title, subtitle) in enumerate(metrics):
        ax = axes[0, idx]

        data = [ceacam_high[col].values, ceacam_low[col].values]
        labels = ['CEACAM-high', 'CEACAM-low']

        bp = ax.boxplot(data, patch_artist=True, labels=labels, widths=0.6)

        for patch, label in zip(bp['boxes'], labels):
            patch.set_facecolor(colors[label])
            patch.set_alpha(0.7)

        # Get p-value - match by column name instead of title
        pval = results_df.iloc[idx]['pvalue']
        sig = '***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else 'ns'))

        ax.set_ylabel(col.replace('_', ' ').title(), fontsize=10)
        ax.set_title(f'{title}\np={pval:.2e} {sig}', fontsize=11, fontweight='bold')

        # Significance bracket
        y_max = max(d.max() for d in data)
        y_range = y_max - min(d.min() for d in data)
        bracket_y = y_max + 0.05 * y_range
        ax.plot([1, 1, 2, 2], [bracket_y, bracket_y + 0.02*y_range, bracket_y + 0.02*y_range, bracket_y],
                color='black', linewidth=1)
        ax.text(1.5, bracket_y + 0.03*y_range, sig, ha='center', va='bottom', fontsize=12, fontweight='bold')

    # Row 2: Violin plots with scatter
    for idx, (col, title, subtitle) in enumerate(metrics):
        ax = axes[1, idx]

        # Sample data for visualization
        n_sample = min(2000, len(ceacam_high), len(ceacam_low))
        high_sample = ceacam_high.sample(n=n_sample, random_state=42)[col].values
        low_sample = ceacam_low.sample(n=n_sample, random_state=42)[col].values

        # Violin plot
        parts = ax.violinplot([high_sample, low_sample], positions=[1, 2], widths=0.7, showmeans=True, showmedians=True)

        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(list(colors.values())[i])
            pc.set_alpha(0.7)

        ax.set_xticks([1, 2])
        ax.set_xticklabels(['CEACAM-high', 'CEACAM-low'])
        ax.set_ylabel(col.replace('_', ' ').title(), fontsize=10)
        ax.set_title(f'{subtitle}', fontsize=10, style='italic')

    plt.suptitle('CEACAM+ Tumor Cells are More Tumor-Centered\n(3 Complementary Spatial Analyses)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = output_dir / 'ceacam_tumor_centered_analysis.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")

    # Summary bar chart
    fig, ax = plt.subplots(figsize=(10, 6))

    y_pos = np.arange(len(results_df))
    cohens_d = results_df['cohens_d'].values
    colors_bar = ['#2ca02c' if d > 0 else '#d62728' for d in cohens_d]

    bars = ax.barh(y_pos, cohens_d, color=colors_bar, edgecolor='black', alpha=0.8)

    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df['metric'])
    ax.set_xlabel("Cohen's d Effect Size\n(Positive = CEACAM-high is more tumor-centered)", fontsize=11)
    ax.set_title('CEACAM-high vs CEACAM-low: Tumor-Centered Position\n(Green = Supports Hypothesis)',
                 fontsize=12, fontweight='bold')

    # Add significance markers
    for i, (_, row) in enumerate(results_df.iterrows()):
        pval = row['pvalue']
        sig = '***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else ''))
        d = row['cohens_d']
        offset = 0.02 if d >= 0 else -0.02
        ha = 'left' if d >= 0 else 'right'
        ax.text(d + offset, i, sig, ha=ha, va='center', fontsize=12, fontweight='bold')

    plt.tight_layout()

    output_path = output_dir / 'ceacam_tumor_centered_summary.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def create_spatial_visualization(df_merged, output_dir, sample_name=None):
    """
    Create spatial maps showing CEACAM ratio vs metrics for one or all samples.

    Args:
        df_merged: Merged dataframe with all metrics
        output_dir: Output directory path
        sample_name: If provided, generate for this sample only; otherwise generate for all samples
    """
    print("\nCreating spatial visualization...")

    if sample_name is not None:
        # Generate for single sample
        samples_to_process = [sample_name]
    else:
        # Generate for all samples
        samples_to_process = sorted(df_merged['sample'].unique())

    print(f"  Processing {len(samples_to_process)} samples...")

    for sample in samples_to_process:
        sample_df = df_merged[df_merged['sample'] == sample].copy()

        if len(sample_df) < 10:
            print(f"  Skipping {sample}: only {len(sample_df)} spots")
            continue

        print(f"  Processing {sample}: {len(sample_df)} spots")

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        # Panel A: CEACAM ratio
        ax = axes[0, 0]
        scatter = ax.scatter(sample_df['x'], sample_df['y'], c=sample_df['CEACAM_ratio'],
                            cmap='RdBu_r', s=8, alpha=0.8)
        plt.colorbar(scatter, ax=ax, label='CEACAM Ratio')
        ax.set_title('A. CEACAM Ratio\n(Red=High, Blue=Low)', fontsize=11, fontweight='bold')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_aspect('equal')

        # Panel B: Neighborhood epithelial density
        ax = axes[0, 1]
        scatter = ax.scatter(sample_df['x'], sample_df['y'], c=sample_df['neighborhood_epi_density'],
                            cmap='YlOrRd', s=8, alpha=0.8)
        plt.colorbar(scatter, ax=ax, label='Neighborhood Epi Density')
        ax.set_title('B. Neighborhood Epithelial Density\n(Yellow=Low, Red=High)', fontsize=11, fontweight='bold')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_aspect('equal')

        # Panel C: Distance to stroma
        ax = axes[1, 0]
        scatter = ax.scatter(sample_df['x'], sample_df['y'], c=sample_df['distance_to_stroma'],
                            cmap='viridis', s=8, alpha=0.8)
        plt.colorbar(scatter, ax=ax, label='Distance to Stroma')
        ax.set_title('C. Distance to Stroma\n(Purple=Close, Yellow=Far)', fontsize=11, fontweight='bold')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_aspect('equal')

        # Panel D: Distance to immune
        ax = axes[1, 1]
        scatter = ax.scatter(sample_df['x'], sample_df['y'], c=sample_df['distance_to_immune'],
                            cmap='plasma', s=8, alpha=0.8)
        plt.colorbar(scatter, ax=ax, label='Distance to Immune')
        ax.set_title('D. Distance to Immune-rich Spots\n(Purple=Close, Yellow=Far)', fontsize=11, fontweight='bold')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_aspect('equal')

        plt.suptitle(f'Spatial Distribution of Tumor-Centered Metrics\n(Sample: {sample})',
                     fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()

        output_path = output_dir / f'ceacam_tumor_centered_spatial_map_{sample}.png'
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {output_path}")


def main():
    print("="*70)
    print("CEACAM Tumor-Centered Analysis")
    print("="*70)

    # Paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'ceacam_tumor_centered'
    output_dir.mkdir(parents=True, exist_ok=True)

    results_dir = script_dir.parent.parent / '02_GraphST_Analysis' / 'results'

    # Load data
    df = load_graphst_results(results_dir)

    # Calculate all three metrics
    neigh_df = calculate_neighborhood_epithelial_density(df, k=K_NEIGHBORS)
    stroma_df = calculate_distance_to_stroma(df, stroma_threshold=STROMA_THRESHOLD)
    immune_df = calculate_distance_to_immune(df, immune_threshold=IMMUNE_THRESHOLD)

    # Merge and analyze
    df_merged = merge_and_analyze(df, neigh_df, stroma_df, immune_df)

    # Statistical comparison
    results_df = statistical_comparison(df_merged)

    # Create visualizations
    create_visualization(df_merged, results_df, output_dir)
    create_spatial_visualization(df_merged, output_dir)

    # Save results
    results_csv = output_dir / 'ceacam_tumor_centered_statistics.csv'
    results_df.to_csv(results_csv, index=False)
    print(f"\nSaved: {results_csv}")

    merged_csv = output_dir / 'ceacam_tumor_centered_spot_data.csv'
    df_merged.to_csv(merged_csv, index=False)
    print(f"Saved: {merged_csv}")

    # Summary
    print("\n" + "="*70)
    print("SUMMARY - Tumor-Centered Analysis")
    print("="*70)

    print("\nHypothesis: CEACAM-high tumor cells are more tumor-centered")
    print("Expected: CEACAM-high spots have HIGHER values for all three metrics")

    cohens_label = "Cohen's d"
    print(f"\n{'Metric':<35} {'Direction':<25} {cohens_label:>10} {'p-value':>12} {'Supports?':>10}")
    print("-"*95)

    all_support = True
    for _, row in results_df.iterrows():
        supports = "YES" if row['cohens_d'] > 0 and row['pvalue'] < 0.05 else "NO"
        if supports == "NO":
            all_support = False
        sig = '***' if row['pvalue'] < 0.001 else ('**' if row['pvalue'] < 0.01 else ('*' if row['pvalue'] < 0.05 else ''))
        print(f"{row['metric']:<35} {row['direction']:<25} {row['cohens_d']:>10.3f} {row['pvalue']:>12.2e} {supports:>10} {sig}")

    print("\n" + "="*70)
    if all_support:
        print("CONCLUSION: All three analyses support the hypothesis!")
        print("CEACAM-high tumor cells ARE more tumor-centered.")
    else:
        print("CONCLUSION: Results are mixed - see details above.")
    print("="*70)

    print("\nDone!")


if __name__ == "__main__":
    main()
