#!/usr/bin/env python3
"""
Extract Cell Type Proportions and Spatial Domains from GraphST Results

This script:
1. Loads GraphST output for all samples
2. Extracts cell type proportions (14 types - CEACAM split)
3. Extracts spatial domain assignments
4. Calculates summary statistics
5. Generates comparison data

Author: Generated for Round 4 Analysis
Date: 2025-12-18
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

# Expected 14 cell types (CEACAM-based epithelial split)
CELL_TYPES_14 = [
    # Immune (8)
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "NK cells",
    "DC cells",
    "Mast cells",
    "Neutrophils",
    "Plasma cells",
    # Myeloid (1) - combined
    "Monocytes/Macrophages",
    # Stromal (3)
    "Endothelial cells",
    "Fibroblast",
    "Pericyte",
    # Epithelial (2) - CEACAM split
    "Epi CEACAM-high",
    "Epi CEACAM-low"
]


def extract_sample_results(h5ad_path, sample_name, cell_types):
    """Extract results from a single sample."""
    print(f"\nProcessing {sample_name}...")

    adata = sc.read_h5ad(h5ad_path)

    # Extract cell type proportions
    proportions = {}
    for ct in cell_types:
        if ct in adata.obs.columns:
            proportions[ct] = adata.obs[ct].values
        else:
            print(f"  Warning: {ct} not found, setting to 0")
            proportions[ct] = np.zeros(adata.n_obs)

    prop_df = pd.DataFrame(proportions, index=adata.obs_names)

    # Extract spatial domain
    if 'spatial_domain' in adata.obs.columns:
        domains = adata.obs['spatial_domain'].copy()
    else:
        print(f"  Warning: spatial_domain not found")
        domains = pd.Series(['unknown'] * adata.n_obs, index=adata.obs_names)

    # Extract spatial coordinates if available
    coords = None
    if 'spatial' in adata.obsm:
        coords = pd.DataFrame(
            adata.obsm['spatial'],
            index=adata.obs_names,
            columns=['x', 'y']
        )

    return prop_df, domains, coords


def calculate_statistics(prop_df_dict, domains_dict):
    """Calculate summary statistics."""
    print("\nCalculating summary statistics...")

    # Mean proportions per sample
    sample_means = {}
    for sample_name, prop_df in prop_df_dict.items():
        sample_means[sample_name] = prop_df.mean()

    mean_df = pd.DataFrame(sample_means).T
    print(f"  Calculated sample-level means")

    # Domain-specific proportions
    domain_stats = {}
    for sample_name in prop_df_dict.keys():
        prop_df = prop_df_dict[sample_name]
        domains = domains_dict[sample_name]

        # Calculate mean per domain
        domain_means = []
        for domain_id in domains.unique():
            mask = domains == domain_id
            domain_mean = prop_df[mask].mean()
            domain_mean.name = f"{sample_name}_domain_{domain_id}"
            domain_means.append(domain_mean)

        domain_stats[sample_name] = pd.DataFrame(domain_means)

    print(f"  Calculated domain-specific means")

    return mean_df, domain_stats


def save_results(output_dir, prop_df_dict, domains_dict, coords_dict, mean_df, domain_stats):
    """Save all extracted results."""
    print("\nSaving results...")

    # Create output directories
    deconv_dir = os.path.join(output_dir, 'deconvolution_matrices')
    domain_dir = os.path.join(output_dir, 'spatial_domains')
    stats_dir = os.path.join(output_dir, 'summary_statistics')

    for d in [deconv_dir, domain_dir, stats_dir]:
        os.makedirs(d, exist_ok=True)

    # Save proportions
    for sample_name, prop_df in prop_df_dict.items():
        path = os.path.join(deconv_dir, f'{sample_name}_proportions_14types.csv')
        prop_df.to_csv(path)

    print(f"  Saved proportions: {deconv_dir}")

    # Save domains
    for sample_name, domains in domains_dict.items():
        path = os.path.join(domain_dir, f'{sample_name}_domains.csv')
        domains.to_csv(path)

    print(f"  Saved domains: {domain_dir}")

    # Save coordinates if available
    if coords_dict:
        for sample_name, coords in coords_dict.items():
            if coords is not None:
                path = os.path.join(domain_dir, f'{sample_name}_coordinates.csv')
                coords.to_csv(path)

    # Save statistics
    mean_df.to_csv(os.path.join(stats_dir, 'mean_proportions_per_sample.csv'))
    print(f"  Saved sample means")

    # Save domain statistics
    for sample_name, df in domain_stats.items():
        path = os.path.join(stats_dir, f'{sample_name}_domain_composition.csv')
        df.to_csv(path)

    print(f"  Saved domain statistics")

    return deconv_dir, domain_dir, stats_dir


def main():
    """Main execution."""
    print("="*80)
    print("Extract GraphST Results - 14 Cell Types (CEACAM Split)")
    print("="*80)

    # Define paths
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent.parent
    results_dir = project_dir / '02_GraphST_Analysis' / 'results'
    output_dir = script_dir.parent / 'results'

    # Find all GraphST output files
    h5ad_files = list(results_dir.glob('*_graphst_output.h5ad'))

    if not h5ad_files:
        print(f"\nError: No GraphST output files found in {results_dir}")
        print("Please run GraphST analysis first.")
        sys.exit(1)

    print(f"\nFound {len(h5ad_files)} GraphST output files")

    # Extract results from all samples
    prop_df_dict = {}
    domains_dict = {}
    coords_dict = {}

    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem.replace('_graphst_output', '')

        prop_df, domains, coords = extract_sample_results(h5ad_path, sample_name, CELL_TYPES_14)

        prop_df_dict[sample_name] = prop_df
        domains_dict[sample_name] = domains
        coords_dict[sample_name] = coords

    # Calculate statistics
    mean_df, domain_stats = calculate_statistics(prop_df_dict, domains_dict)

    # Save results
    deconv_dir, domain_dir, stats_dir = save_results(
        output_dir, prop_df_dict, domains_dict, coords_dict, mean_df, domain_stats
    )

    print("\n" + "="*80)
    print("Extraction Complete!")
    print("="*80)
    print(f"\nOutput directories:")
    print(f"  Proportions: {deconv_dir}")
    print(f"  Domains: {domain_dir}")
    print(f"  Statistics: {stats_dir}")

    print(f"\nExtracted data for {len(prop_df_dict)} samples")
    print(f"  Cell types: {len(CELL_TYPES_14)}")
    print(f"  Total spots: {sum(len(df) for df in prop_df_dict.values()):,}")

    print(f"\nNext steps:")
    print(f"  - Run visualization scripts")
    print(f"  - Compare with marker-based results")


if __name__ == "__main__":
    main()
