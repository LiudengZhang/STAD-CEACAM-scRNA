#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Region-Level Analysis: CEACAM vs Immune Cell Infiltration

Aggregates spots by GraphST spatial_domain to compare immune infiltration
between CEACAM-high vs CEACAM-low dominant TUMOR REGIONS.

Approach:
1. Group spots by spatial_domain (GraphST clusters)
2. Filter to tumor regions (≥30% epithelial content)
3. Classify tumor regions as CEACAM-high or CEACAM-low dominant
4. Compare immune cell proportions between region types

Author: Generated for Round 4 Analysis
Date: 2026-01-04
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
plt.rcParams['font.family'] = 'sans-serif'


def load_graphst_results(results_dir):
    """Load all GraphST results with spatial domain info."""
    print("Loading GraphST results...")
    h5ad_files = list(Path(results_dir).glob('*_graphst_output.h5ad'))

    adata_dict = {}
    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem.replace('_graphst_output', '')
        adata = sc.read_h5ad(h5ad_path)
        adata_dict[sample_name] = adata
        print(f"  {sample_name}: {adata.n_obs} spots, {adata.obs['spatial_domain'].nunique()} domains")

    return adata_dict


def aggregate_to_regions(adata_dict, min_spots=10, min_epi=0.30):
    """
    Aggregate spots to region level (by spatial_domain).
    Filter to tumor regions with sufficient epithelial content.

    Args:
        min_spots: Minimum spots per region
        min_epi: Minimum total epithelial proportion to be a "tumor region"
    """
    print(f"\nAggregating to REGION level (min_spots={min_spots}, min_epi={min_epi})...")

    immune_cells = [
        'B cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells',
        'DC cells', 'Mast cells', 'Neutrophils', 'Plasma cells',
        'Monocytes/Macrophages'
    ]

    all_regions = []

    for sample_name, adata in adata_dict.items():
        domains = adata.obs['spatial_domain'].unique()

        for domain in domains:
            mask = adata.obs['spatial_domain'] == domain
            n_spots = mask.sum()

            if n_spots < min_spots:
                continue

            # Calculate region means
            region_data = {
                'sample': sample_name,
                'domain': domain,
                'n_spots': n_spots,
                'mean_CEACAM_high': adata.obs.loc[mask, 'Epi CEACAM-high'].mean(),
                'mean_CEACAM_low': adata.obs.loc[mask, 'Epi CEACAM-low'].mean(),
            }

            # Total epithelial
            region_data['total_epi'] = region_data['mean_CEACAM_high'] + region_data['mean_CEACAM_low']

            # CEACAM delta and ratio
            region_data['CEACAM_delta'] = region_data['mean_CEACAM_high'] - region_data['mean_CEACAM_low']
            if region_data['total_epi'] > 0:
                region_data['CEACAM_ratio'] = region_data['mean_CEACAM_high'] / region_data['total_epi']
            else:
                region_data['CEACAM_ratio'] = 0.5

            # Immune cell proportions
            for immune in immune_cells:
                region_data[f'mean_{immune}'] = adata.obs.loc[mask, immune].mean()

            # Total immune
            region_data['total_immune'] = sum(region_data[f'mean_{immune}'] for immune in immune_cells)

            all_regions.append(region_data)

    region_df = pd.DataFrame(all_regions)
    print(f"  Total regions: {len(region_df)}")

    # Filter to tumor regions
    tumor_regions = region_df[region_df['total_epi'] >= min_epi].copy()
    print(f"  Tumor regions (≥{min_epi*100:.0f}% epithelial): {len(tumor_regions)}")

    # Classify by CEACAM dominance (median split within tumor regions)
    median_ratio = tumor_regions['CEACAM_ratio'].median()
    tumor_regions['region_type'] = np.where(
        tumor_regions['CEACAM_ratio'] >= median_ratio,
        'CEACAM-high dominant',
        'CEACAM-low dominant'
    )

    print(f"  CEACAM-high dominant: {(tumor_regions['region_type'] == 'CEACAM-high dominant').sum()}")
    print(f"  CEACAM-low dominant: {(tumor_regions['region_type'] == 'CEACAM-low dominant').sum()}")
    print(f"  Median CEACAM_ratio threshold: {median_ratio:.3f}")

    return region_df, tumor_regions


def compare_immune_by_region_type(tumor_regions, output_dir):
    """
    Compare immune cell infiltration between CEACAM-high vs CEACAM-low tumor regions.
    """
    print("\n" + "="*70)
    print("Comparing immune infiltration by tumor region type...")
    print("="*70)

    immune_cells = [
        'B cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells',
        'DC cells', 'Mast cells', 'Neutrophils', 'Plasma cells',
        'Monocytes/Macrophages'
    ]

    ceacam_high = tumor_regions[tumor_regions['region_type'] == 'CEACAM-high dominant']
    ceacam_low = tumor_regions[tumor_regions['region_type'] == 'CEACAM-low dominant']

    results = []

    for immune in immune_cells:
        col = f'mean_{immune}'
        high_vals = ceacam_high[col].values
        low_vals = ceacam_low[col].values

        # Mann-Whitney U test
        stat, pval = stats.mannwhitneyu(high_vals, low_vals, alternative='two-sided')

        # Effect size (difference in medians)
        diff = np.median(high_vals) - np.median(low_vals)

        results.append({
            'immune_cell': immune,
            'CEACAM_high_median': np.median(high_vals),
            'CEACAM_low_median': np.median(low_vals),
            'difference': diff,
            'pvalue': pval,
            'higher_in': 'CEACAM-high' if diff > 0 else 'CEACAM-low'
        })

        print(f"  {immune}: diff={diff:.4f}, p={pval:.4f}")

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('difference', ascending=False)

    return results_df


def plot_immune_comparison_boxplots(tumor_regions, results_df, output_dir):
    """
    Create boxplots comparing immune infiltration between region types.
    """
    print("\nGenerating immune comparison boxplots...")

    immune_cells = results_df['immune_cell'].tolist()

    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    axes = axes.flatten()

    ceacam_high = tumor_regions[tumor_regions['region_type'] == 'CEACAM-high dominant']
    ceacam_low = tumor_regions[tumor_regions['region_type'] == 'CEACAM-low dominant']

    for idx, immune in enumerate(immune_cells):
        ax = axes[idx]
        col = f'mean_{immune}'

        data = [ceacam_high[col].values, ceacam_low[col].values]
        labels = ['CEACAM-high\nregions', 'CEACAM-low\nregions']
        colors = ['#d62728', '#1f77b4']

        bp = ax.boxplot(data, patch_artist=True, labels=labels, widths=0.6)

        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # Add individual points
        for i, (d, c) in enumerate(zip(data, colors), 1):
            x = np.random.normal(i, 0.08, len(d))
            ax.scatter(x, d, alpha=0.4, color=c, s=15, edgecolor='white', linewidth=0.3)

        # Get p-value
        row = results_df[results_df['immune_cell'] == immune].iloc[0]
        pval = row['pvalue']

        sig = ''
        if pval < 0.001:
            sig = '***'
        elif pval < 0.01:
            sig = '**'
        elif pval < 0.05:
            sig = '*'

        ax.set_ylabel(f'{immune}\nProportion', fontsize=10)
        ax.set_title(f'{immune}\np={pval:.3f} {sig}', fontsize=11, fontweight='bold')

        # Add significance bracket
        y_max = max(d.max() for d in data)
        y_range = y_max - min(d.min() for d in data)
        bracket_y = y_max + 0.05 * y_range
        ax.plot([1, 1, 2, 2], [bracket_y, bracket_y + 0.02*y_range, bracket_y + 0.02*y_range, bracket_y],
                color='black', linewidth=1)

    plt.suptitle('Immune Cell Infiltration: CEACAM-high vs CEACAM-low Tumor Regions\n(Region-Level Analysis)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = output_dir / 'ceacam_immune_region_boxplots.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_summary_barplot(results_df, output_dir):
    """
    Summary bar plot showing direction and significance.
    """
    print("\nGenerating summary bar plot...")

    fig, ax = plt.subplots(figsize=(10, 8))

    # Sort by difference
    df_sorted = results_df.sort_values('difference', ascending=True)

    y_pos = np.arange(len(df_sorted))
    colors = ['#d62728' if d > 0 else '#1f77b4' for d in df_sorted['difference']]

    bars = ax.barh(y_pos, df_sorted['difference'], color=colors, edgecolor='black', alpha=0.8)

    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_sorted['immune_cell'])
    ax.set_xlabel('Difference in Median Proportion\n(CEACAM-high regions - CEACAM-low regions)', fontsize=11)
    ax.set_title('Immune Recruitment by CEACAM Tumor Region Type\n'
                 '(Red = Higher in CEACAM-high, Blue = Higher in CEACAM-low)',
                 fontsize=12, fontweight='bold')

    # Add significance markers
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        pval = row['pvalue']
        sig = ''
        if pval < 0.001:
            sig = '***'
        elif pval < 0.01:
            sig = '**'
        elif pval < 0.05:
            sig = '*'

        x_pos = row['difference']
        offset = 0.002 if x_pos >= 0 else -0.002
        ha = 'left' if x_pos >= 0 else 'right'
        ax.text(x_pos + offset, i, sig, ha=ha, va='center', fontsize=12, fontweight='bold')

    plt.tight_layout()

    output_path = output_dir / 'ceacam_immune_region_summary.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_region_scatter(tumor_regions, output_dir):
    """
    Scatter plot: CEACAM_ratio vs each immune cell (region level).
    """
    print("\nGenerating region-level scatter plots...")

    immune_cells = [
        'B cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells',
        'DC cells', 'Mast cells', 'Neutrophils', 'Plasma cells',
        'Monocytes/Macrophages'
    ]

    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    axes = axes.flatten()

    for idx, immune in enumerate(immune_cells):
        ax = axes[idx]
        col = f'mean_{immune}'

        # Color by sample
        samples = tumor_regions['sample'].unique()
        cmap = plt.cm.tab10(np.linspace(0, 1, len(samples)))
        sample_colors = dict(zip(samples, cmap))

        for sample in samples:
            subset = tumor_regions[tumor_regions['sample'] == sample]
            ax.scatter(subset['CEACAM_ratio'], subset[col],
                       c=[sample_colors[sample]], s=40, alpha=0.6,
                       edgecolor='white', linewidth=0.3)

        # Correlation
        rho, pval = stats.spearmanr(tumor_regions['CEACAM_ratio'], tumor_regions[col])

        # Trend line
        z = np.polyfit(tumor_regions['CEACAM_ratio'], tumor_regions[col], 1)
        p = np.poly1d(z)
        x_line = np.linspace(tumor_regions['CEACAM_ratio'].min(), tumor_regions['CEACAM_ratio'].max(), 100)
        ax.plot(x_line, p(x_line), 'k--', alpha=0.7, linewidth=1.5)

        ax.set_xlabel('CEACAM Ratio\n(0=CEACAM-low, 1=CEACAM-high)', fontsize=9)
        ax.set_ylabel(f'{immune}', fontsize=10)

        sig = ''
        if pval < 0.001:
            sig = '***'
        elif pval < 0.01:
            sig = '**'
        elif pval < 0.05:
            sig = '*'

        ax.set_title(f'{immune}\nρ={rho:.3f}, p={pval:.3f} {sig}', fontsize=10, fontweight='bold')

    plt.suptitle('Region-Level Correlation: CEACAM Ratio vs Immune Infiltration',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = output_dir / 'ceacam_immune_region_scatter.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    print("="*70)
    print("Region-Level Analysis: CEACAM vs Immune Infiltration")
    print("="*70)

    # Paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'ceacam_immune_region_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)

    results_dir = script_dir.parent.parent / '02_GraphST_Analysis' / 'results'

    # Load GraphST results
    adata_dict = load_graphst_results(results_dir)

    # Aggregate to region level
    region_df, tumor_regions = aggregate_to_regions(adata_dict, min_spots=10, min_epi=0.30)

    # Compare immune infiltration
    results_df = compare_immune_by_region_type(tumor_regions, output_dir)

    # Generate plots
    plot_immune_comparison_boxplots(tumor_regions, results_df, output_dir)
    plot_summary_barplot(results_df, output_dir)
    plot_region_scatter(tumor_regions, output_dir)

    # Save results
    results_csv = output_dir / 'ceacam_immune_region_comparison.csv'
    results_df.to_csv(results_csv, index=False)
    print(f"\nSaved: {results_csv}")

    tumor_csv = output_dir / 'tumor_regions_data.csv'
    tumor_regions.to_csv(tumor_csv, index=False)
    print(f"Saved: {tumor_csv}")

    # Summary
    print("\n" + "="*70)
    print("SUMMARY - Region-Level Analysis")
    print("="*70)
    print(f"\nTumor regions analyzed: {len(tumor_regions)}")
    print(f"  CEACAM-high dominant: {(tumor_regions['region_type'] == 'CEACAM-high dominant').sum()}")
    print(f"  CEACAM-low dominant: {(tumor_regions['region_type'] == 'CEACAM-low dominant').sum()}")

    print(f"\n{'Immune Cell':<25} {'Higher in':>15} {'Difference':>12} {'p-value':>12}")
    print("-"*65)
    for _, row in results_df.iterrows():
        sig = '*' if row['pvalue'] < 0.05 else ''
        print(f"{row['immune_cell']:<25} {row['higher_in']:>15} {row['difference']:>12.4f} {row['pvalue']:>12.4f} {sig}")

    # Key findings
    print("\n" + "="*70)
    print("KEY FINDINGS:")
    print("="*70)

    sig_high = results_df[(results_df['difference'] > 0) & (results_df['pvalue'] < 0.05)]
    sig_low = results_df[(results_df['difference'] < 0) & (results_df['pvalue'] < 0.05)]

    if len(sig_high) > 0:
        print(f"\nSignificantly HIGHER in CEACAM-high tumor regions:")
        for _, row in sig_high.iterrows():
            print(f"  - {row['immune_cell']} (p={row['pvalue']:.4f})")
    else:
        print("\nNo immune cells significantly higher in CEACAM-high regions")

    if len(sig_low) > 0:
        print(f"\nSignificantly HIGHER in CEACAM-low tumor regions:")
        for _, row in sig_low.iterrows():
            print(f"  - {row['immune_cell']} (p={row['pvalue']:.4f})")
    else:
        print("\nNo immune cells significantly higher in CEACAM-low regions")

    print("\n" + "="*70)
    print("Done!")
    print("="*70)


if __name__ == "__main__":
    main()
