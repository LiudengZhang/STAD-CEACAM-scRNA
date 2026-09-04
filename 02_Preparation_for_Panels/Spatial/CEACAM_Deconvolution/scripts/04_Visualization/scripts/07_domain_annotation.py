#!/usr/bin/env python3
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
"""
Domain Annotation - Assign Spatial Domains to Biological Categories

Classifies each GraphST spatial domain into one of 6 categories:
1. Stromal
2. Immune
3. Normal Epi
4. F10-high Tumor
5. F10-low Tumor
6. Mixed (if no category >= 40%)

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

# Cell type to category mapping (using actual column names from GraphST output)
CATEGORY_MAPPING = {
    'Stromal': ['Fibroblast', 'Endothelial cells', 'Pericyte'],
    'Immune': ['B cells', 'DC cells', 'Mast cells', 'Monocytes/Macrophages',
               'NK cells', 'Plasma cells', 'CD4+ T cells', 'CD8+ T cells', 'Neutrophils'],
    'Normal Epi': ['Normal Epi'],
    'F10-high Tumor': ['Cancer F10-high'],
    'F10-low Tumor': ['Cancer F10-low']
}

# Domain column name in GraphST output
DOMAIN_COL = 'spatial_domain'

# Color scheme
CATEGORY_COLORS = {
    'Stromal': '#FF8C00',       # Orange
    'Immune': '#1E90FF',        # Blue
    'Normal Epi': '#32CD32',    # Green
    'F10-high Tumor': '#DC143C', # Red
    'F10-low Tumor': '#FF69B4', # Pink
    'Mixed': '#808080'          # Gray
}

MIXED_THRESHOLD = 0.40  # 40%


def annotate_domains(adata):
    """Annotate domains based on cell type proportions."""

    # Get all cell type columns (the 15 types)
    cell_types = list(set(ct for cts in CATEGORY_MAPPING.values() for ct in cts))

    # Check which cell types are present
    available_cts = [ct for ct in cell_types if ct in adata.obs.columns]

    if len(available_cts) == 0:
        print("  Warning: No cell type columns found!")
        return None

    # Get unique domains
    domains = adata.obs[DOMAIN_COL].unique()

    results = []
    for domain in domains:
        mask = adata.obs[DOMAIN_COL] == domain
        n_spots = mask.sum()

        # Calculate mean proportion for each category
        category_props = {}
        for cat, cts in CATEGORY_MAPPING.items():
            cat_prop = 0
            for ct in cts:
                if ct in adata.obs.columns:
                    cat_prop += adata.obs.loc[mask, ct].mean()
            category_props[cat] = cat_prop

        # Find dominant category
        max_cat = max(category_props, key=category_props.get)
        max_prop = category_props[max_cat]

        # Apply threshold
        if max_prop >= MIXED_THRESHOLD:
            assigned_cat = max_cat
        else:
            assigned_cat = 'Mixed'

        results.append({
            'domain': domain,
            'n_spots': n_spots,
            'assigned_category': assigned_cat,
            'max_proportion': max_prop,
            **{f'prop_{cat}': prop for cat, prop in category_props.items()}
        })

    return pd.DataFrame(results)


def main():
    print("="*70)
    print("Domain Annotation - Assign Domains to Biological Categories")
    print("="*70)
    print(f"\nThreshold for Mixed: {MIXED_THRESHOLD*100:.0f}%")

    # Paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'domain_annotation'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load GraphST results
    results_dir = Path('/path/to/home/Project_4_Gastric_Cancer/01.2_GraphST_17Types/03_Results')
    h5ad_files = sorted(results_dir.glob('*_graphst_output.h5ad'))

    print(f"\nLoading {len(h5ad_files)} samples...")

    # Store all annotations
    all_annotations = []
    sample_data = []

    for h5ad_path in h5ad_files:
        sample_name = h5ad_path.stem.replace('_graphst_output', '')
        print(f"\n  Processing {sample_name}...")

        adata = sc.read_h5ad(h5ad_path)

        # Annotate domains
        df = annotate_domains(adata)
        if df is not None:
            df['sample'] = sample_name
            all_annotations.append(df)
            sample_data.append((sample_name, adata))

            # Print summary
            print(f"    Domains: {len(df)}")
            for cat in list(CATEGORY_COLORS.keys()):
                count = (df['assigned_category'] == cat).sum()
                if count > 0:
                    print(f"      {cat}: {count}")

    # Combine all annotations
    all_df = pd.concat(all_annotations, ignore_index=True)

    # Save CSV
    csv_path = output_dir / 'domain_assignments.csv'
    all_df.to_csv(csv_path, index=False)
    print(f"\n\nSaved: {csv_path}")

    # Create visualization: 10-panel figure
    fig, axes = plt.subplots(2, 5, figsize=(25, 10))
    axes = axes.flatten()

    for idx, (sample_name, adata) in enumerate(sample_data):
        ax = axes[idx]

        # Get domain annotations for this sample
        sample_df = all_df[all_df['sample'] == sample_name]

        # Create color mapping for domains
        domain_to_cat = dict(zip(sample_df['domain'], sample_df['assigned_category']))

        # Get spatial coordinates
        coords = adata.obsm['spatial']
        domains = adata.obs[DOMAIN_COL].values

        # Create color array
        colors = [CATEGORY_COLORS[domain_to_cat.get(d, 'Mixed')] for d in domains]

        # Plot
        ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=8, alpha=0.8)
        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.set_title(f'{sample_name}\n{len(sample_df)} domains', fontsize=11)
        ax.axis('off')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, label=cat)
                       for cat, color in CATEGORY_COLORS.items()]
    fig.legend(handles=legend_elements, loc='lower center', ncol=6,
               fontsize=11, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle('Spatial Domains Annotated by Biological Category\n(GraphST - 17 Types)',
                 fontsize=16, fontweight='bold', y=1.02)

    plt.tight_layout()

    plot_path = output_dir / 'annotated_domains_all_samples.png'
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()

    print(f"Saved: {plot_path}")

    # Summary statistics
    print("\n" + "="*70)
    print("Overall Summary")
    print("="*70)
    print(f"\nTotal domains: {len(all_df)}")
    print("\nCategory distribution:")
    for cat in CATEGORY_COLORS.keys():
        count = (all_df['assigned_category'] == cat).sum()
        pct = count / len(all_df) * 100
        print(f"  {cat}: {count} ({pct:.1f}%)")

    print("\n" + "="*70)
    print("Done!")
    print("="*70)


if __name__ == "__main__":
    main()
