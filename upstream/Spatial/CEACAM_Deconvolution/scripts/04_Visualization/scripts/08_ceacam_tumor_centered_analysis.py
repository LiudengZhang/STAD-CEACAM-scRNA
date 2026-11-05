#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
CEACAM Tumor-Centered Analysis
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

K_NEIGHBORS = 15  # k for k-NN neighborhood analysis
EPI_THRESHOLD = 0.05  # Minimum epithelial to be "tumor spot"
STROMA_THRESHOLD = 0.20  # Maximum epithelial to be "stroma-dominant"
IMMUNE_THRESHOLD = 0.15  # Minimum total immune to be "immune-rich"


def load_graphst_results(results_dir):
    h5ad_files = list(Path(results_dir).glob('*_graphst_output.h5ad'))

    all_data = []
    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem.replace('_graphst_output', '')
        adata = sc.read_h5ad(h5ad_path)

        df = adata.obs.copy()
        df['sample'] = sample_name
        df['x'] = adata.obsm['spatial'][:, 0]
        df['y'] = adata.obsm['spatial'][:, 1]
        df['spot_id'] = df.index

        all_data.append(df)

    df_all = pd.concat(all_data, ignore_index=True)
    return df_all


def calculate_neighborhood_epithelial_density(df, k=K_NEIGHBORS):

    results = []

    for sample in tqdm(df['sample'].unique(), desc="Processing samples"):
        sample_df = df[df['sample'] == sample].copy()

        coords = sample_df[['x', 'y']].values
        tree = cKDTree(coords)

        distances, indices = tree.query(coords, k=k+1)

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

    results = []

    for sample in tqdm(df['sample'].unique(), desc="Processing samples"):
        sample_df = df[df['sample'] == sample].copy()

        sample_df['Total_Epi'] = sample_df['Epi CEACAM-high'] + sample_df['Epi CEACAM-low']

        stroma_mask = sample_df['Total_Epi'] < stroma_threshold
        stroma_spots = sample_df[stroma_mask]

        if len(stroma_spots) < 5:
            continue

        stroma_coords = stroma_spots[['x', 'y']].values
        stroma_tree = cKDTree(stroma_coords)

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

    immune_cols = ['B cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells',
                   'DC cells', 'Mast cells', 'Neutrophils', 'Plasma cells',
                   'Monocytes/Macrophages']

    results = []

    for sample in tqdm(df['sample'].unique(), desc="Processing samples"):
        sample_df = df[df['sample'] == sample].copy()

        sample_df['Total_Immune'] = sample_df[immune_cols].sum(axis=1)

        immune_mask = sample_df['Total_Immune'] > immune_threshold
        immune_spots = sample_df[immune_mask]

        if len(immune_spots) < 5:
            continue

        immune_coords = immune_spots[['x', 'y']].values
        immune_tree = cKDTree(immune_coords)

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

    df['Total_Epi'] = df['Epi CEACAM-high'] + df['Epi CEACAM-low']
    df_tumor = df[df['Total_Epi'] >= epi_threshold].copy()
    df_tumor['CEACAM_ratio'] = df_tumor['Epi CEACAM-high'] / df_tumor['Total_Epi']

    df_merged = df_tumor.merge(neigh_df, on=['sample', 'spot_id'], how='left')
    df_merged = df_merged.merge(stroma_df, on=['sample', 'spot_id'], how='left')
    df_merged = df_merged.merge(immune_df, on=['sample', 'spot_id'], how='left')

    df_merged = df_merged.dropna(subset=['neighborhood_epi_density', 'distance_to_stroma', 'distance_to_immune'])

    median_ratio = df_merged['CEACAM_ratio'].median()
    df_merged['CEACAM_group'] = np.where(
        df_merged['CEACAM_ratio'] >= median_ratio,
        'CEACAM-high',
        'CEACAM-low'
    )

    return df_merged


def statistical_comparison(df_merged):

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

        stat, pval = stats.mannwhitneyu(high_vals, low_vals, alternative='two-sided')

        pooled_std = np.sqrt((high_vals.std()**2 + low_vals.std()**2) / 2)
        cohens_d = (high_vals.mean() - low_vals.mean()) / pooled_std if pooled_std > 0 else 0

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

    return pd.DataFrame(results)




def main():

    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'ceacam_tumor_centered'
    output_dir.mkdir(parents=True, exist_ok=True)

    results_dir = script_dir.parent.parent / '02_GraphST_Analysis' / 'results'

    df = load_graphst_results(results_dir)

    neigh_df = calculate_neighborhood_epithelial_density(df, k=K_NEIGHBORS)
    stroma_df = calculate_distance_to_stroma(df, stroma_threshold=STROMA_THRESHOLD)
    immune_df = calculate_distance_to_immune(df, immune_threshold=IMMUNE_THRESHOLD)

    df_merged = merge_and_analyze(df, neigh_df, stroma_df, immune_df)

    results_df = statistical_comparison(df_merged)


    results_csv = output_dir / 'ceacam_tumor_centered_statistics.csv'
    results_df.to_csv(results_csv, index=False)

    merged_csv = output_dir / 'ceacam_tumor_centered_spot_data.csv'
    df_merged.to_csv(merged_csv, index=False)

    cohens_label = "Cohen's d"

    all_support = True
    for _, row in results_df.iterrows():
        supports = "YES" if row['cohens_d'] > 0 and row['pvalue'] < 0.05 else "NO"
        if supports == "NO":
            all_support = False
        sig = '***' if row['pvalue'] < 0.001 else ('**' if row['pvalue'] < 0.01 else ('*' if row['pvalue'] < 0.05 else ''))

    if all_support:
        print("CONCLUSION: All three analyses support the hypothesis!")
        print("CEACAM-high tumor cells ARE more tumor-centered.")
    else:
        print("CONCLUSION: Results are mixed - see details above.")


if __name__ == "__main__":
    main()
