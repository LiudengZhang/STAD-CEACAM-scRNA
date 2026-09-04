#!/usr/bin/env python3
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
"""
Plot Epithelial Cell Type Delta Comparisons - Spatial Maps

Creates 2 spatial figures showing delta at each spot:
1. Tumor vs Normal Epi Delta: (mean(F10-high, F10-low)) - Normal Epi
2. F10-high vs F10-low Delta: F10-high - F10-low

Author: Generated for Round 4 Analysis
Date: 2025-11-26
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set plotting parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


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


def plot_delta_spatial(adata_dict, delta_name, delta_values_dict, output_path, title):
    """Plot spatial distribution of delta values."""
    n_samples = len(adata_dict)
    ncols = 5
    nrows = int(np.ceil(n_samples / ncols))

    fig = plt.figure(figsize=(30, 18))
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.3, wspace=0.3)

    # Find global min/max for consistent colorbar
    all_deltas = np.concatenate([delta_values_dict[s] for s in delta_values_dict])
    vmax = max(abs(all_deltas.min()), abs(all_deltas.max()))
    vmin = -vmax  # Symmetric around 0

    for idx, (sample_name, adata) in enumerate(adata_dict.items()):
        row = idx // ncols
        col = idx % ncols
        ax = fig.add_subplot(gs[row, col])

        deltas = delta_values_dict[sample_name]
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

        # Plot deltas with diverging colormap (blue = negative, red = positive)
        scatter = ax.scatter(coords_to_use[:, 0], coords_to_use[:, 1],
                           c=deltas, cmap='RdBu_r', s=5,
                           vmin=vmin, vmax=vmax, alpha=0.9)

        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Delta', rotation=270, labelpad=15, fontsize=8)
        cbar.ax.tick_params(labelsize=7)

        # Statistics
        mean_delta = deltas.mean()
        pos_pct = (deltas > 0).sum() / len(deltas) * 100

        ax.set_title(f'{sample_name}\nMean: {mean_delta:.3f}, Pos: {pos_pct:.0f}%',
                    fontsize=10)
        ax.axis('off')
        ax.invert_yaxis()

    fig.suptitle(title, fontsize=16, fontweight='bold')

    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def main():
    print("="*70)
    print("Plotting Epithelial Cell Type Deltas - Spatial Maps")
    print("="*70)

    # Paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'summary_plots'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load GraphST results from rsrch8
    results_dir = Path('/path/to/home/Project_4_Gastric_Cancer/01.2_GraphST_17Types/03_Results')
    adata_dict = load_graphst_results(results_dir)

    print(f"\nLoaded {len(adata_dict)} samples")

    # Calculate deltas for each sample
    print("\nCalculating deltas...")
    delta_tumor_vs_normal = {}
    delta_f10high_vs_low = {}
    delta_il2high_vs_low = {}

    for sample_name, adata in adata_dict.items():
        # Get proportions
        f10_high = adata.obs['Cancer F10-high'].values
        f10_low = adata.obs['Cancer F10-low'].values
        normal_epi = adata.obs['Normal Epi'].values
        il2_high = adata.obs['IL2_high_Macrophage'].values
        il2_low = adata.obs['IL2_low_Macrophage'].values

        # Calculate deltas
        tumor_avg = (f10_high + f10_low) / 2
        delta_tumor_vs_normal[sample_name] = tumor_avg - normal_epi
        delta_f10high_vs_low[sample_name] = f10_high - f10_low
        delta_il2high_vs_low[sample_name] = il2_high - il2_low

    # Figure 1: Tumor vs Normal Epi Delta
    print("\nPlotting Figure 1: Tumor vs Normal Epi Delta (Spatial)...")
    plot_delta_spatial(
        adata_dict,
        "Tumor vs Normal",
        delta_tumor_vs_normal,
        output_dir / 'tumor_vs_normal_epi_delta_spatial.png',
        'Tumor vs Normal Epithelial Delta\n(Tumor = avg of Cancer F10-high and F10-low)'
    )

    # Figure 2: F10-high vs F10-low Delta
    print("\nPlotting Figure 2: F10-high vs F10-low Delta (Spatial)...")
    plot_delta_spatial(
        adata_dict,
        "F10-high vs F10-low",
        delta_f10high_vs_low,
        output_dir / 'f10_high_vs_low_delta_spatial.png',
        'Cancer F10-high vs F10-low Delta'
    )

    # Figure 3: IL2_high vs IL2_low Macrophage Delta
    print("\nPlotting Figure 3: IL2_high vs IL2_low Macrophage Delta (Spatial)...")
    plot_delta_spatial(
        adata_dict,
        "IL2_high vs IL2_low Mac",
        delta_il2high_vs_low,
        output_dir / 'il2_high_vs_low_macrophage_delta_spatial.png',
        'IL2_high vs IL2_low Macrophage Delta'
    )

    # Print summary
    print("\n" + "="*70)
    print("Summary Statistics")
    print("="*70)

    print("\nTumor vs Normal Epi Delta (per spot):")
    all_d1 = np.concatenate([delta_tumor_vs_normal[s] for s in delta_tumor_vs_normal])
    print(f"  Overall Mean: {all_d1.mean():.4f}")
    print(f"  Range: [{all_d1.min():.4f}, {all_d1.max():.4f}]")
    print(f"  Positive spots: {(all_d1 > 0).sum()}/{len(all_d1)} ({(all_d1 > 0).mean()*100:.1f}%)")

    print("\nF10-high vs F10-low Delta (per spot):")
    all_d2 = np.concatenate([delta_f10high_vs_low[s] for s in delta_f10high_vs_low])
    print(f"  Overall Mean: {all_d2.mean():.4f}")
    print(f"  Range: [{all_d2.min():.4f}, {all_d2.max():.4f}]")
    print(f"  Positive spots: {(all_d2 > 0).sum()}/{len(all_d2)} ({(all_d2 > 0).mean()*100:.1f}%)")

    print("\nIL2_high vs IL2_low Macrophage Delta (per spot):")
    all_d3 = np.concatenate([delta_il2high_vs_low[s] for s in delta_il2high_vs_low])
    print(f"  Overall Mean: {all_d3.mean():.4f}")
    print(f"  Range: [{all_d3.min():.4f}, {all_d3.max():.4f}]")
    print(f"  Positive spots: {(all_d3 > 0).sum()}/{len(all_d3)} ({(all_d3 > 0).mean()*100:.1f}%)")

    print("\n" + "="*70)
    print("Done!")
    print("="*70)


if __name__ == "__main__":
    main()
