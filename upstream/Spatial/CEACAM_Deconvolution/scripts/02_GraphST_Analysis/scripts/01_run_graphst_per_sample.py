#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
GraphST Spatial Deconvolution - Single Sample
15 Cell Types (12 Non-Epithelial + 3 Epithelial) - Korean Gut Dataset

This script runs GraphST on a single Visium sample:
1. Loads spatial data and single-cell reference
2. Preprocesses both datasets
3. Constructs spatial graphs
4. Trains GraphST model with deconvolution
5. Projects cell types to spots
6. Identifies spatial domains
7. Saves results

Author: Generated for Round 4 Analysis
Date: 2025-11-21
"""

import os
import sys
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Force unbuffered output for real-time logging
sys.stdout.reconfigure(line_buffering=True)

# Import GraphST
try:
    import GraphST
    from GraphST.GraphST import GraphST as GraphSTModel
    from GraphST.utils import project_cell_to_spot
    from GraphST.preprocess import filter_with_overlap_gene
except ImportError as e:
    print(f"Error: GraphST not found. Please ensure conda environment is activated.")
    print(f"  conda activate stad_ceacam")
    sys.exit(1)

# Set scanpy settings
sc.settings.verbosity = 1


def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def load_visium_sample(sample_path, sample_name):
    """
    Load a single Visium sample from Korean Gut format.

    Parameters:
    -----------
    sample_path : str
        Path to the sample directory
    sample_name : str
        Name of the sample

    Returns:
    --------
    adata : AnnData
        Loaded spatial data
    """
    print(f"\n{'='*80}")
    print(f"Loading Visium Sample: {sample_name}")
    print(f"{'='*80}")
    print(f"Path: {sample_path}")

    sample_dir = Path(sample_path)

    # Find and load H5 matrix
    h5_files = list(sample_dir.glob('*filtered_feature_bc_matrix.h5'))
    if not h5_files:
        raise FileNotFoundError(f"No H5 file found in {sample_path}")

    adata = sc.read_10x_h5(h5_files[0])
    print(f"  Loaded expression matrix: {adata.n_obs} spots × {adata.n_vars} genes")

    # Load spatial coordinates
    tissue_pos_files = list(sample_dir.glob('*tissue_positions_list.csv'))
    if not tissue_pos_files:
        raise FileNotFoundError(f"No tissue positions found in {sample_path}")

    tissue_positions = pd.read_csv(tissue_pos_files[0], header=None)
    tissue_positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col',
                                'pxl_row_in_fullres', 'pxl_col_in_fullres']
    tissue_positions.set_index('barcode', inplace=True)

    # Align barcodes
    common_barcodes = adata.obs_names.intersection(tissue_positions.index)
    adata = adata[common_barcodes].copy()
    tissue_positions = tissue_positions.loc[common_barcodes]

    # Add spatial coordinates
    adata.obsm['spatial'] = tissue_positions[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values

    # Load images and scale factors
    scalefactor_files = list(sample_dir.glob('*scalefactors_json.json'))
    if scalefactor_files:
        import json
        with open(scalefactor_files[0], 'r') as f:
            scalefactors = json.load(f)
    else:
        scalefactors = {}

    hires_img_files = list(sample_dir.glob('*tissue_hires_image.png'))
    lowres_img_files = list(sample_dir.glob('*tissue_lowres_image.png'))

    # Create spatial metadata
    library_id = sample_name
    adata.uns['spatial'] = {library_id: {}}

    if hires_img_files or lowres_img_files:
        adata.uns['spatial'][library_id]['images'] = {}

        if hires_img_files:
            import matplotlib.image as mpimg
            adata.uns['spatial'][library_id]['images']['hires'] = mpimg.imread(hires_img_files[0])

        if lowres_img_files:
            import matplotlib.image as mpimg
            adata.uns['spatial'][library_id]['images']['lowres'] = mpimg.imread(lowres_img_files[0])

    if scalefactors:
        adata.uns['spatial'][library_id]['scalefactors'] = scalefactors

    adata.obs['sample'] = sample_name
    adata.var_names_make_unique()

    print(f"  Final: {adata.n_obs} spots × {adata.n_vars} genes")

    return adata


def preprocess_spatial(adata):
    """
    Preprocess spatial data using GraphST.preprocess().

    This handles normalization, log1p, and HVG selection (3000 genes).
    """
    print(f"\n{'='*80}")
    print("Preprocessing Spatial Data (GraphST)")
    print(f"{'='*80}")

    # Store raw counts first
    adata.layers['counts'] = adata.X.copy()
    print(f"  ✓ Saved raw counts")

    # Ensure gene names are unique
    adata.var_names_make_unique()
    print(f"  ✓ Made gene names unique")

    # Use GraphST's preprocessing (normalize, log1p, HVG selection)
    print(f"\n  Running GraphST.preprocess()...")
    GraphST.preprocess(adata)
    print(f"  ✓ Normalized, log-transformed, identified 3000 HVGs")

    return adata


def construct_spatial_graph(adata, config):
    """Construct spatial interaction graph."""
    print(f"\n{'='*80}")
    print("Constructing Spatial Graph")
    print(f"{'='*80}")

    graph_params = config['spatial_graph']

    # Calculate neighbors
    GraphST.construct_interaction(
        adata,
        n_neighbors=graph_params['n_neighbors']
    )

    print(f"  ✓ Built spatial graph with {graph_params['n_neighbors']} neighbors")

    return adata


def add_contrastive_labels(adata):
    """Add labels for contrastive learning."""
    print(f"\n{'='*80}")
    print("Adding Contrastive Learning Labels")
    print(f"{'='*80}")

    GraphST.add_contrastive_label(adata)

    print(f"  ✓ Added contrastive labels")

    return adata


def preprocess_reference(adata_sc, cell_type_column):
    """
    Preprocess reference data manually.

    IMPORTANT: Cannot use GraphST.preprocess() because it hardcodes n_top_genes=3000
    but our reference only has ~2,900 genes. Instead, manually normalize, log-transform,
    and mark all genes as highly_variable.

    Parameters:
    -----------
    adata_sc : AnnData
        Reference single-cell data
    cell_type_column : str
        Column name for cell type annotations (from config)
    """
    print(f"\n{'='*80}")
    print("Preprocessing Reference Data (Manual)")
    print(f"{'='*80}")

    # Check if raw counts are in layers
    if 'counts' in adata_sc.layers:
        print(f"  ✓ Found raw counts in adata_sc.layers['counts']")
        # Ensure main matrix has raw counts
        adata_sc.X = adata_sc.layers['counts'].copy()
        print(f"  ✓ Restored raw counts to adata_sc.X")
    else:
        print(f"  ✓ Using existing adata_sc.X (assumed to be raw counts)")

    # Ensure gene names are unique
    adata_sc.var_names_make_unique()
    print(f"  ✓ Made gene names unique")

    # FIX: Ensure cell names are unique (GraphST construct_cell_type_matrix fails with duplicates)
    adata_sc.obs_names_make_unique()
    print(f"  ✓ Made cell names unique")

    # Manual preprocessing (can't use GraphST.preprocess due to < 3000 genes)
    print(f"\n  Manual preprocessing steps:")

    # Normalize
    sc.pp.normalize_total(adata_sc, target_sum=1e4)
    print(f"  ✓ Normalized to 10,000 counts per cell")

    # Log transform
    sc.pp.log1p(adata_sc)
    print(f"  ✓ Log-transformed (log1p)")

    # Mark all genes as highly variable (we only have ~2900 genes)
    adata_sc.var['highly_variable'] = True
    adata_sc.var['highly_variable_rank'] = range(adata_sc.n_vars)
    print(f"  ✓ Marked all {adata_sc.n_vars} genes as highly variable")

    print(f"\n  NOTE: Cannot use seurat_v3 HVG with < 3000 genes")
    print(f"  All genes marked as HVG to satisfy filter_with_overlap_gene()")

    # FIX: GraphST hardcodes 'cell_type' column - create it from config-specified column
    # Must happen BEFORE training, not after
    adata_sc.obs['cell_type'] = adata_sc.obs[cell_type_column].astype(str)
    print(f"  ✓ Created 'cell_type' column from '{cell_type_column}' ({adata_sc.obs['cell_type'].nunique()} types)")

    return adata_sc


def filter_genes(adata, adata_sc):
    """
    Filter to overlapping genes.

    IMPORTANT: This must happen BEFORE graph construction to ensure
    the spatial graph uses the final gene set.
    """
    print(f"\n{'='*80}")
    print("Filtering to Overlapping Genes")
    print(f"{'='*80}")

    print(f"  Before filtering:")
    print(f"    Spatial: {adata.n_vars} genes")
    print(f"    Reference: {adata_sc.n_vars} genes")

    adata, adata_sc = filter_with_overlap_gene(adata, adata_sc)

    print(f"\n  After filtering:")
    print(f"    Spatial: {adata.n_vars} genes")
    print(f"    Reference: {adata_sc.n_vars} genes")
    print(f"  ✓ Filtered to {adata.n_vars} overlapping genes")

    return adata, adata_sc


def train_graphst(adata, adata_sc, config):
    """Train GraphST model."""
    print(f"\n{'='*80}")
    print("Training GraphST Model")
    print(f"{'='*80}")

    params = config['graphst_parameters']

    print(f"  Configuration:")
    print(f"    Epochs: {params['epochs']}")
    print(f"    Learning rate: {params['learning_rate']}")
    print(f"    Device: {params['device']}")
    print(f"    Deconvolution: {params['deconvolution']}")
    print(f"    Random seed: {params['random_seed']}")

    # Initialize model
    model = GraphSTModel(
        adata,
        adata_sc,
        random_seed=params['random_seed'],
        epochs=params['epochs'],
        learning_rate=params['learning_rate'],
        deconvolution=params['deconvolution'],
        device=params['device']
    )

    print(f"\n  Starting training...")
    print(f"  (This may take 15-25 minutes on CPU)")

    # Train
    import time
    start_time = time.time()

    adata, adata_sc = model.train_map()

    elapsed = time.time() - start_time
    print(f"\n  ✓ Training completed in {elapsed/60:.1f} minutes")

    return adata, adata_sc


def project_cells(adata, adata_sc, config):
    """Project reference cells to spatial spots."""
    print(f"\n{'='*80}")
    print("Projecting Cells to Spots")
    print(f"{'='*80}")

    deconv_params = config['deconvolution']
    cell_type_col = config['reference']['cell_type_column']

    # Note: 'cell_type' column already created in preprocess_reference()

    project_cell_to_spot(
        adata,
        adata_sc,
        retain_percent=deconv_params['retain_percent']
    )

    # Check which cell types were mapped
    cell_type_col = config['reference']['cell_type_column']
    cell_types = adata_sc.obs[cell_type_col].unique()

    print(f"  ✓ Projected {len(cell_types)} cell types:")
    for ct in sorted(cell_types):
        if ct in adata.obs.columns:
            mean_prop = adata.obs[ct].mean()
            print(f"    {ct}: mean={mean_prop:.3f}")

    return adata


def identify_domains(adata, config):
    """Identify spatial domains using clustering."""
    print(f"\n{'='*80}")
    print("Identifying Spatial Domains")
    print(f"{'='*80}")

    clust_params = config['clustering']

    # Use GraphST embeddings for clustering
    # GraphST stores embeddings in 'emb_sp', not 'GraphST'
    sc.pp.neighbors(adata, use_rep='emb_sp')

    # Leiden clustering
    sc.tl.leiden(adata, resolution=clust_params['resolution'], key_added='spatial_domain')

    n_domains = adata.obs['spatial_domain'].nunique()
    print(f"  ✓ Identified {n_domains} spatial domains")

    domain_counts = adata.obs['spatial_domain'].value_counts()
    for domain, count in domain_counts.items():
        print(f"    Domain {domain}: {count} spots")

    return adata


def save_results(adata, sample_name, config):
    """Save processed results."""
    print(f"\n{'='*80}")
    print("Saving Results")
    print(f"{'='*80}")

    output_dir = os.path.join(
        os.path.dirname(config_path),
        config['output']['results_dir']
    )
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f'{sample_name}_graphst_output.h5ad')

    adata.write(output_path)

    file_size = os.path.getsize(output_path) / 1e6
    print(f"  ✓ Saved: {output_path}")
    print(f"  File size: {file_size:.1f} MB")

    return output_path


def detect_device(config_device):
    """
    Detect and configure compute device (GPU vs CPU).

    Returns the actual device to use.
    """
    import torch

    if config_device == "auto":
        if torch.cuda.is_available():
            device = "cuda:0"
            print(f"\n{'='*80}")
            print("GPU DETECTED!")
            print(f"{'='*80}")
            print(f"  Device: {torch.cuda.get_device_name(0)}")
            print(f"  CUDA Version: {torch.version.cuda}")
            print(f"  Using GPU for training (5-10x speedup expected)")
            print(f"{'='*80}\n")
        else:
            device = "cpu"
            print(f"\n{'='*80}")
            print("No GPU detected - using CPU")
            print(f"{'='*80}")
            print(f"  Training will be slower (~15-25 min per sample)")
            print(f"  Consider using GPU for faster processing")
            print(f"{'='*80}\n")
    else:
        device = config_device
        print(f"\n  Using configured device: {device}")

    return device


def main(config_path, sample_name):
    """Main execution function."""
    print("="*80)
    print("GraphST Spatial Deconvolution - Single Sample")
    print("15 Cell Types (12 Non-Epithelial + 3 Epithelial)")
    print("="*80)

    # Load configuration
    config = load_config(config_path)

    # Auto-detect device
    config['graphst_parameters']['device'] = detect_device(
        config['graphst_parameters']['device']
    )

    # Check if sample exists
    if sample_name not in config['samples']:
        print(f"\nError: Sample '{sample_name}' not found in config")
        print(f"Available samples: {', '.join(config['samples'].keys())}")
        sys.exit(1)

    sample_info = config['samples'][sample_name]

    # Load single-cell reference
    print(f"\n{'='*80}")
    print("Loading Single-Cell Reference")
    print(f"{'='*80}")

    ref_path = os.path.join(os.path.dirname(config_path), config['reference']['path'])
    print(f"  Path: {ref_path}")

    if not os.path.exists(ref_path):
        print(f"\n  Error: Reference not found!")
        print(f"  Please run reference preparation first:")
        print(f"    cd 01_Reference_Preparation")
        print(f"    python scripts/01_prepare_stomach_reference_13types.py config/reference_config.yaml")
        sys.exit(1)

    adata_sc = sc.read_h5ad(ref_path)
    print(f"  Loaded: {adata_sc.n_obs} cells × {adata_sc.n_vars} genes")

    cell_type_col = config['reference']['cell_type_column']
    n_types = adata_sc.obs[cell_type_col].nunique()
    print(f"  Cell types: {n_types}")

    # Load spatial sample
    adata = load_visium_sample(sample_info['path'], sample_name)

    # === WORKFLOW ORDER ===
    # 1. Preprocess spatial with GraphST (normalize, log1p, HVG)
    adata = preprocess_spatial(adata)

    # 2. Preprocess reference with GraphST (normalize, log1p, HVG)
    adata_sc = preprocess_reference(adata_sc, cell_type_col)

    # 3. Filter to overlapping HVGs (requires 'highly_variable' to be set)
    #    This ensures graph uses final gene set
    adata, adata_sc = filter_genes(adata, adata_sc)

    # 4. NOW construct spatial graph (with final genes)
    adata = construct_spatial_graph(adata, config)

    # 5. Add contrastive labels
    adata = add_contrastive_labels(adata)

    # 6. Train GraphST (will call GraphST.preprocess internally)
    adata, adata_sc = train_graphst(adata, adata_sc, config)

    # Project cells to spots
    adata = project_cells(adata, adata_sc, config)

    # Identify spatial domains
    adata = identify_domains(adata, config)

    # Save results
    output_path = save_results(adata, sample_name, config)

    print(f"\n{'='*80}")
    print("GraphST Analysis Complete!")
    print(f"{'='*80}")
    print(f"\nOutput: {output_path}")
    print(f"\nNext steps:")
    print(f"  - Run on remaining samples")
    print(f"  - Extract cell type proportions")
    print(f"  - Visualize spatial domains and cell types")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 01_run_graphst_per_sample.py <config_file> <sample_name>")
        print("\nExample:")
        print("  python 01_run_graphst_per_sample.py config/graphst_config.yaml sample_01")
        sys.exit(1)

    config_path = sys.argv[1]
    sample_name = sys.argv[2]

    if not os.path.exists(config_path):
        print(f"Error: Configuration file not found: {config_path}")
        sys.exit(1)

    main(config_path, sample_name)
