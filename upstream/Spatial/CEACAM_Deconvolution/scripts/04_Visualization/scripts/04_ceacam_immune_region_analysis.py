#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Region-Level Analysis: CEACAM vs Immune Cell Infiltration

Approach:
1. Group spots by spatial_domain (GraphST clusters)
2. Filter to tumor regions (≥30% epithelial content)
3. Classify tumor regions as CEACAM-high or CEACAM-low dominant
4. Compare immune cell proportions between region types
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
    h5ad_files = list(Path(results_dir).glob('*_graphst_output.h5ad'))

    adata_dict = {}
    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem.replace('_graphst_output', '')
        adata = sc.read_h5ad(h5ad_path)
        adata_dict[sample_name] = adata

    return adata_dict


def aggregate_to_regions(adata_dict, min_spots=10, min_epi=0.30):

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

            region_data = {
                'sample': sample_name,
                'domain': domain,
                'n_spots': n_spots,
                'mean_CEACAM_high': adata.obs.loc[mask, 'Epi CEACAM-high'].mean(),
                'mean_CEACAM_low': adata.obs.loc[mask, 'Epi CEACAM-low'].mean(),
            }

            region_data['total_epi'] = region_data['mean_CEACAM_high'] + region_data['mean_CEACAM_low']

            region_data['CEACAM_delta'] = region_data['mean_CEACAM_high'] - region_data['mean_CEACAM_low']
            if region_data['total_epi'] > 0:
                region_data['CEACAM_ratio'] = region_data['mean_CEACAM_high'] / region_data['total_epi']
            else:
                region_data['CEACAM_ratio'] = 0.5

            for immune in immune_cells:
                region_data[f'mean_{immune}'] = adata.obs.loc[mask, immune].mean()

            region_data['total_immune'] = sum(region_data[f'mean_{immune}'] for immune in immune_cells)

            all_regions.append(region_data)

    region_df = pd.DataFrame(all_regions)

    tumor_regions = region_df[region_df['total_epi'] >= min_epi].copy()

    median_ratio = tumor_regions['CEACAM_ratio'].median()
    tumor_regions['region_type'] = np.where(
        tumor_regions['CEACAM_ratio'] >= median_ratio,
        'CEACAM-high dominant',
        'CEACAM-low dominant'
    )

    return region_df, tumor_regions


def compare_immune_by_region_type(tumor_regions, output_dir):

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

        stat, pval = stats.mannwhitneyu(high_vals, low_vals, alternative='two-sided')

        diff = np.median(high_vals) - np.median(low_vals)

        results.append({
            'immune_cell': immune,
            'CEACAM_high_median': np.median(high_vals),
            'CEACAM_low_median': np.median(low_vals),
            'difference': diff,
            'pvalue': pval,
            'higher_in': 'CEACAM-high' if diff > 0 else 'CEACAM-low'
        })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('difference', ascending=False)

    return results_df




def main():

    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'ceacam_immune_region_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)

    results_dir = script_dir.parent.parent / '02_GraphST_Analysis' / 'results'

    adata_dict = load_graphst_results(results_dir)

    region_df, tumor_regions = aggregate_to_regions(adata_dict, min_spots=10, min_epi=0.30)

    results_df = compare_immune_by_region_type(tumor_regions, output_dir)


    results_csv = output_dir / 'ceacam_immune_region_comparison.csv'
    results_df.to_csv(results_csv, index=False)

    tumor_csv = output_dir / 'tumor_regions_data.csv'
    tumor_regions.to_csv(tumor_csv, index=False)
    for _, row in results_df.iterrows():
        sig = '*' if row['pvalue'] < 0.05 else ''

    sig_high = results_df[(results_df['difference'] > 0) & (results_df['pvalue'] < 0.05)]
    sig_low = results_df[(results_df['difference'] < 0) & (results_df['pvalue'] < 0.05)]

    if len(sig_high) > 0:
        for _, row in sig_high.iterrows():
            print(f"  - {row['immune_cell']} (p={row['pvalue']:.4f})")
    else:
        print("\nNo immune cells significantly higher in CEACAM-high regions")

    if len(sig_low) > 0:
        for _, row in sig_low.iterrows():
            print(f"  - {row['immune_cell']} (p={row['pvalue']:.4f})")
    else:
        print("\nNo immune cells significantly higher in CEACAM-low regions")


if __name__ == "__main__":
    main()
