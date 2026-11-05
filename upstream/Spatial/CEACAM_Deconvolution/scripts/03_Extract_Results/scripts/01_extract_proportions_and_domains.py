#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Extract Cell Type Proportions and Spatial Domains from GraphST Results
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

CELL_TYPES_14 = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "NK cells",
    "DC cells",
    "Mast cells",
    "Neutrophils",
    "Plasma cells",
    "Monocytes/Macrophages",
    "Endothelial cells",
    "Fibroblast",
    "Pericyte",
    "Epi CEACAM-high",
    "Epi CEACAM-low"
]


def extract_sample_results(h5ad_path, sample_name, cell_types):

    adata = sc.read_h5ad(h5ad_path)

    proportions = {}
    for ct in cell_types:
        if ct in adata.obs.columns:
            proportions[ct] = adata.obs[ct].values
        else:
            proportions[ct] = np.zeros(adata.n_obs)

    prop_df = pd.DataFrame(proportions, index=adata.obs_names)

    if 'spatial_domain' in adata.obs.columns:
        domains = adata.obs['spatial_domain'].copy()
    else:
        domains = pd.Series(['unknown'] * adata.n_obs, index=adata.obs_names)

    coords = None
    if 'spatial' in adata.obsm:
        coords = pd.DataFrame(
            adata.obsm['spatial'],
            index=adata.obs_names,
            columns=['x', 'y']
        )

    return prop_df, domains, coords


def calculate_statistics(prop_df_dict, domains_dict):

    sample_means = {}
    for sample_name, prop_df in prop_df_dict.items():
        sample_means[sample_name] = prop_df.mean()

    mean_df = pd.DataFrame(sample_means).T

    domain_stats = {}
    for sample_name in prop_df_dict.keys():
        prop_df = prop_df_dict[sample_name]
        domains = domains_dict[sample_name]

        domain_means = []
        for domain_id in domains.unique():
            mask = domains == domain_id
            domain_mean = prop_df[mask].mean()
            domain_mean.name = f"{sample_name}_domain_{domain_id}"
            domain_means.append(domain_mean)

        domain_stats[sample_name] = pd.DataFrame(domain_means)

    return mean_df, domain_stats


def save_results(output_dir, prop_df_dict, domains_dict, coords_dict, mean_df, domain_stats):

    deconv_dir = os.path.join(output_dir, 'deconvolution_matrices')
    domain_dir = os.path.join(output_dir, 'spatial_domains')
    stats_dir = os.path.join(output_dir, 'summary_statistics')

    for d in [deconv_dir, domain_dir, stats_dir]:
        os.makedirs(d, exist_ok=True)

    for sample_name, prop_df in prop_df_dict.items():
        path = os.path.join(deconv_dir, f'{sample_name}_proportions_14types.csv')
        prop_df.to_csv(path)

    for sample_name, domains in domains_dict.items():
        path = os.path.join(domain_dir, f'{sample_name}_domains.csv')
        domains.to_csv(path)

    if coords_dict:
        for sample_name, coords in coords_dict.items():
            if coords is not None:
                path = os.path.join(domain_dir, f'{sample_name}_coordinates.csv')
                coords.to_csv(path)

    mean_df.to_csv(os.path.join(stats_dir, 'mean_proportions_per_sample.csv'))

    for sample_name, df in domain_stats.items():
        path = os.path.join(stats_dir, f'{sample_name}_domain_composition.csv')
        df.to_csv(path)

    return deconv_dir, domain_dir, stats_dir


def main():

    script_dir = Path(__file__).parent
    project_dir = script_dir.parent.parent
    results_dir = project_dir / '02_GraphST_Analysis' / 'results'
    output_dir = script_dir.parent / 'results'

    h5ad_files = list(results_dir.glob('*_graphst_output.h5ad'))

    if not h5ad_files:
        sys.exit(1)

    prop_df_dict = {}
    domains_dict = {}
    coords_dict = {}

    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem.replace('_graphst_output', '')

        prop_df, domains, coords = extract_sample_results(h5ad_path, sample_name, CELL_TYPES_14)

        prop_df_dict[sample_name] = prop_df
        domains_dict[sample_name] = domains
        coords_dict[sample_name] = coords

    mean_df, domain_stats = calculate_statistics(prop_df_dict, domains_dict)

    deconv_dir, domain_dir, stats_dir = save_results(
        output_dir, prop_df_dict, domains_dict, coords_dict, mean_df, domain_stats
    )


if __name__ == "__main__":
    main()
