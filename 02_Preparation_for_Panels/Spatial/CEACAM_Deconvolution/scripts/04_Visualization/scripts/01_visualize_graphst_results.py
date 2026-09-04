#!/usr/bin/env python3
"""
Comprehensive Visualization of GraphST Results
14 Cell Types (CEACAM Split) - Korean Gut Dataset

Creates:
1. Spatial domain maps for all samples
2. Cell type proportion maps (14 types)
3. Dominant cell type maps
4. Summary heatmaps
5. Domain composition plots

Author: Generated for Round 4 Analysis
Date: 2025-12-18
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set plotting parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# 14 Cell Types (CEACAM-based epithelial split)
CELL_TYPES_14 = [
    'B cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells', 'DC cells',
    'Mast cells', 'Neutrophils', 'Plasma cells',
    'Monocytes/Macrophages',
    'Endothelial cells', 'Fibroblast', 'Pericyte',
    'Epi CEACAM-high', 'Epi CEACAM-low'
]


def load_graphst_results(results_dir):
    """Load all GraphST results."""
    print("Loading GraphST results...")

    h5ad_files = list(Path(results_dir).glob('*_graphst_output.h5ad'))

    if not h5ad_files:
        raise FileNotFoundError(f"No GraphST results found in {results_dir}")

    adata_dict = {}
    for h5ad_path in sorted(h5ad_files):
        sample_name = h5ad_path.stem.replace('_graphst_output', '')
        adata = sc.read_h5ad(h5ad_path)
        adata_dict[sample_name] = adata
        print(f"  {sample_name}: {adata.n_obs} spots")

    return adata_dict


def plot_spatial_domains(adata_dict, output_dir):
    """Plot spatial domains for all samples."""
    print("\nGenerating spatial domain maps...")

    n_samples = len(adata_dict)
    ncols = 5
    nrows = int(np.ceil(n_samples / ncols))

    fig = plt.figure(figsize=(30, 18))
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.3, wspace=0.3)

    for idx, (sample_name, adata) in enumerate(adata_dict.items()):
        row = idx // ncols
        col = idx % ncols
        ax = fig.add_subplot(gs[row, col])

        # Get spatial coordinates and domains
        spatial_coords = adata.obsm['spatial']
        domains = adata.obs['spatial_domain'].values

        # Get unique domains and assign colors
        unique_domains = sorted(domains.unique())
        colors = plt.cm.tab20(np.linspace(0, 1, len(unique_domains)))
        domain_colors = {d: colors[i] for i, d in enumerate(unique_domains)}

        # Plot tissue image if available
        img = None
        if 'spatial' in adata.uns:
            library_id = list(adata.uns['spatial'].keys())[0]
            if 'images' in adata.uns['spatial'][library_id]:
                if 'hires' in adata.uns['spatial'][library_id]['images']:
                    img = adata.uns['spatial'][library_id]['images']['hires']
                    img = np.transpose(img, (1, 0, 2))  # Transpose y,x axes to align with spatial coords
                    scale = adata.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']
                    spatial_coords_scaled = spatial_coords * scale
                    ax.imshow(img, alpha=0.6)
                    coords_to_use = spatial_coords_scaled
                else:
                    coords_to_use = spatial_coords
            else:
                coords_to_use = spatial_coords
        else:
            coords_to_use = spatial_coords

        # Plot domains
        spot_colors = [domain_colors[d] for d in domains]
        ax.scatter(coords_to_use[:, 0], coords_to_use[:, 1],
                  c=spot_colors, s=5, alpha=0.8)

        ax.set_title(f'{sample_name}\n{len(unique_domains)} domains', fontsize=10)
        ax.axis('off')
        ax.invert_yaxis()

    fig.suptitle('Spatial Domains Across Korean Gut Samples (GraphST - 14 Types CEACAM)',
                 fontsize=16, fontweight='bold')

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'spatial_domains_all_samples.png')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def plot_cell_type_proportions(adata_dict, cell_type, output_dir):
    """Plot spatial distribution of a single cell type."""
    n_samples = len(adata_dict)
    ncols = 5
    nrows = int(np.ceil(n_samples / ncols))

    fig = plt.figure(figsize=(30, 18))
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.3, wspace=0.3)

    for idx, (sample_name, adata) in enumerate(adata_dict.items()):
        row = idx // ncols
        col = idx % ncols
        ax = fig.add_subplot(gs[row, col])

        # Get proportions
        if cell_type not in adata.obs.columns:
            continue

        proportions = adata.obs[cell_type].values
        spatial_coords = adata.obsm['spatial']

        # Get tissue image
        coords_to_use = spatial_coords
        if 'spatial' in adata.uns:
            library_id = list(adata.uns['spatial'].keys())[0]
            if 'images' in adata.uns['spatial'][library_id]:
                if 'hires' in adata.uns['spatial'][library_id]['images']:
                    img = adata.uns['spatial'][library_id]['images']['hires']
                    img = np.transpose(img, (1, 0, 2))  # Transpose y,x axes to align with spatial coords
                    scale = adata.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']
                    coords_to_use = spatial_coords * scale
                    ax.imshow(img, alpha=0.6)

        # Plot proportions
        scatter = ax.scatter(coords_to_use[:, 0], coords_to_use[:, 1],
                           c=proportions, cmap='Reds', s=5,
                           vmin=0, vmax=1, alpha=0.9)

        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Proportion', rotation=270, labelpad=15, fontsize=8)
        cbar.ax.tick_params(labelsize=7)

        # Statistics
        mean_prop = proportions.mean()
        max_prop = proportions.max()

        ax.set_title(f'{sample_name}\nMean: {mean_prop:.3f}, Max: {max_prop:.3f}',
                    fontsize=10)
        ax.axis('off')
        ax.invert_yaxis()

    fig.suptitle(f'{cell_type} Spatial Distribution (GraphST - 14 Types CEACAM)',
                 fontsize=16, fontweight='bold')

    # Sanitize filename (replace / with _)
    safe_name = cell_type.replace('/', '_')
    output_path = os.path.join(output_dir, f'{safe_name}_spatial_map.png')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()


def plot_all_cell_types(adata_dict, cell_types, output_dir):
    """Plot all cell type proportion maps."""
    print("\nGenerating cell type proportion maps...")

    prop_dir = os.path.join(output_dir, 'proportion_maps')
    os.makedirs(prop_dir, exist_ok=True)

    for cell_type in cell_types:
        print(f"  Plotting {cell_type}...")
        plot_cell_type_proportions(adata_dict, cell_type, prop_dir)

    print(f"  Saved {len(cell_types)} cell type maps")


def plot_summary_heatmap(adata_dict, cell_types, output_dir):
    """Plot summary heatmap of mean proportions."""
    print("\nGenerating summary heatmap...")

    # Collect mean proportions
    sample_names = list(adata_dict.keys())
    prop_matrix = []

    for sample_name in sample_names:
        adata = adata_dict[sample_name]
        means = []
        for ct in cell_types:
            if ct in adata.obs.columns:
                means.append(adata.obs[ct].mean())
            else:
                means.append(0)
        prop_matrix.append(means)

    prop_df = pd.DataFrame(prop_matrix, index=sample_names, columns=cell_types)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(18, 10))

    sns.heatmap(prop_df, cmap='YlOrRd', annot=True, fmt='.2f',
               cbar_kws={'label': 'Mean Proportion'}, ax=ax)

    ax.set_title('Mean Cell Type Proportions Across Samples (GraphST - 14 Types CEACAM)',
                fontsize=14, fontweight='bold')
    ax.set_xlabel('Cell Type', fontsize=12)
    ax.set_ylabel('Sample', fontsize=12)

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    output_path = os.path.join(output_dir, 'mean_proportions_heatmap.png')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")

    # Also save CSV
    csv_path = os.path.join(output_dir, 'mean_proportions_summary.csv')
    prop_df.to_csv(csv_path)
    print(f"  Saved: {csv_path}")


def main():
    """Main execution."""
    print("="*80)
    print("Visualize GraphST Results - 14 Cell Types (CEACAM Split)")
    print("="*80)

    # Define paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results'

    # Results from 14Types CEACAM analysis
    project_dir = script_dir.parent.parent
    results_dir = project_dir / '02_GraphST_Analysis' / 'results'

    # Load results
    adata_dict = load_graphst_results(results_dir)

    print(f"\nLoaded {len(adata_dict)} samples")

    # Create visualizations
    domain_dir = output_dir / 'domain_maps'
    os.makedirs(domain_dir, exist_ok=True)

    # Spatial domains
    plot_spatial_domains(adata_dict, domain_dir)

    # Cell type proportions
    plot_all_cell_types(adata_dict, CELL_TYPES_14, output_dir)

    # Summary heatmap
    summary_dir = output_dir / 'summary_plots'
    os.makedirs(summary_dir, exist_ok=True)
    plot_summary_heatmap(adata_dict, CELL_TYPES_14, summary_dir)

    print("\n" + "="*80)
    print("Visualization Complete!")
    print("="*80)
    print(f"\nOutput directories:")
    print(f"  Domain maps: {domain_dir}")
    print(f"  Proportion maps: {output_dir / 'proportion_maps'}")
    print(f"  Summary plots: {summary_dir}")


if __name__ == "__main__":
    main()
