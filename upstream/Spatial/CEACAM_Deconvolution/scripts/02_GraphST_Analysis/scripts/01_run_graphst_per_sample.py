#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
GraphST Spatial Deconvolution - Single Sample
15 Cell Types (12 Non-Epithelial + 3 Epithelial) - Korean Gut Dataset
"""

import os
import sys
import yaml
import pandas as pd
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

sys.stdout.reconfigure(line_buffering=True)

try:
    import GraphST
    from GraphST.GraphST import GraphST as GraphSTModel
    from GraphST.utils import project_cell_to_spot
    from GraphST.preprocess import filter_with_overlap_gene
except ImportError as e:
    sys.exit(1)

sc.settings.verbosity = 1


def load_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def load_visium_sample(sample_path, sample_name):

    sample_dir = Path(sample_path)

    h5_files = list(sample_dir.glob('*filtered_feature_bc_matrix.h5'))
    if not h5_files:
        raise FileNotFoundError(f"No H5 file found in {sample_path}")

    adata = sc.read_10x_h5(h5_files[0])

    tissue_pos_files = list(sample_dir.glob('*tissue_positions_list.csv'))
    if not tissue_pos_files:
        raise FileNotFoundError(f"No tissue positions found in {sample_path}")

    tissue_positions = pd.read_csv(tissue_pos_files[0], header=None)
    tissue_positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col',
                                'pxl_row_in_fullres', 'pxl_col_in_fullres']
    tissue_positions.set_index('barcode', inplace=True)

    common_barcodes = adata.obs_names.intersection(tissue_positions.index)
    adata = adata[common_barcodes].copy()
    tissue_positions = tissue_positions.loc[common_barcodes]

    adata.obsm['spatial'] = tissue_positions[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values

    scalefactor_files = list(sample_dir.glob('*scalefactors_json.json'))
    if scalefactor_files:
        import json
        with open(scalefactor_files[0], 'r') as f:
            scalefactors = json.load(f)
    else:
        scalefactors = {}

    hires_img_files = list(sample_dir.glob('*tissue_hires_image.png'))
    lowres_img_files = list(sample_dir.glob('*tissue_lowres_image.png'))

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

    return adata


def preprocess_spatial(adata):
    """
    Preprocess spatial data using GraphST.preprocess().

    This handles normalization, log1p, and HVG selection (3000 genes).
    """

    # Store raw counts first
    adata.layers['counts'] = adata.X.copy()

    adata.var_names_make_unique()
    GraphST.preprocess(adata)

    return adata


def construct_spatial_graph(adata, config):

    graph_params = config['spatial_graph']

    GraphST.construct_interaction(
        adata,
        n_neighbors=graph_params['n_neighbors']
    )

    return adata


def add_contrastive_labels(adata):

    GraphST.add_contrastive_label(adata)

    return adata


def preprocess_reference(adata_sc, cell_type_column):

    # Check if raw counts are in layers
    if 'counts' in adata_sc.layers:
        # Ensure main matrix has raw counts
        adata_sc.X = adata_sc.layers['counts'].copy()
    else:
        print(f"  ✓ Using existing adata_sc.X (assumed to be raw counts)")

    adata_sc.var_names_make_unique()

    adata_sc.obs_names_make_unique()

    sc.pp.normalize_total(adata_sc, target_sum=1e4)

    sc.pp.log1p(adata_sc)

    adata_sc.var['highly_variable'] = True
    adata_sc.var['highly_variable_rank'] = range(adata_sc.n_vars)

    # FIX: GraphST hardcodes 'cell_type' column - create it from config-specified column
    # Must happen BEFORE training, not after
    adata_sc.obs['cell_type'] = adata_sc.obs[cell_type_column].astype(str)

    return adata_sc


def filter_genes(adata, adata_sc):
    """
    Filter to overlapping genes.

    IMPORTANT: This must happen BEFORE graph construction to ensure
    the spatial graph uses the final gene set.
    """

    adata, adata_sc = filter_with_overlap_gene(adata, adata_sc)

    return adata, adata_sc


def train_graphst(adata, adata_sc, config):

    params = config['graphst_parameters']

    model = GraphSTModel(
        adata,
        adata_sc,
        random_seed=params['random_seed'],
        epochs=params['epochs'],
        learning_rate=params['learning_rate'],
        deconvolution=params['deconvolution'],
        device=params['device']
    )

    import time
    start_time = time.time()

    adata, adata_sc = model.train_map()

    elapsed = time.time() - start_time

    return adata, adata_sc


def project_cells(adata, adata_sc, config):

    deconv_params = config['deconvolution']
    cell_type_col = config['reference']['cell_type_column']

    project_cell_to_spot(
        adata,
        adata_sc,
        retain_percent=deconv_params['retain_percent']
    )

    cell_type_col = config['reference']['cell_type_column']
    cell_types = adata_sc.obs[cell_type_col].unique()

    for ct in sorted(cell_types):
        if ct in adata.obs.columns:
            mean_prop = adata.obs[ct].mean()

    return adata


def identify_domains(adata, config):

    clust_params = config['clustering']

    sc.pp.neighbors(adata, use_rep='emb_sp')

    sc.tl.leiden(adata, resolution=clust_params['resolution'], key_added='spatial_domain')

    n_domains = adata.obs['spatial_domain'].nunique()

    domain_counts = adata.obs['spatial_domain'].value_counts()
    for domain, count in domain_counts.items():
        print(f"    Domain {domain}: {count} spots")

    return adata


def save_results(adata, sample_name, config):

    output_dir = os.path.join(
        os.path.dirname(config_path),
        config['output']['results_dir']
    )
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f'{sample_name}_graphst_output.h5ad')

    adata.write(output_path)

    file_size = os.path.getsize(output_path) / 1e6

    return output_path


def detect_device(config_device):
    import torch

    if config_device == "auto":
        if torch.cuda.is_available():
            device = "cuda:0"
        else:
            device = "cpu"
    else:
        device = config_device

    return device


def main(config_path, sample_name):

    config = load_config(config_path)

    config['graphst_parameters']['device'] = detect_device(
        config['graphst_parameters']['device']
    )

    if sample_name not in config['samples']:
        sys.exit(1)

    sample_info = config['samples'][sample_name]

    ref_path = os.path.join(os.path.dirname(config_path), config['reference']['path'])

    if not os.path.exists(ref_path):
        sys.exit(1)

    adata_sc = sc.read_h5ad(ref_path)

    cell_type_col = config['reference']['cell_type_column']
    n_types = adata_sc.obs[cell_type_col].nunique()

    adata = load_visium_sample(sample_info['path'], sample_name)

    # === WORKFLOW ORDER ===
    # 1. Preprocess spatial with GraphST (normalize, log1p, HVG)
    adata = preprocess_spatial(adata)

    # 2. Preprocess reference with GraphST (normalize, log1p, HVG)
    adata_sc = preprocess_reference(adata_sc, cell_type_col)

    # 3. Filter to overlapping HVGs (requires 'highly_variable' to be set)
    #    This ensures graph uses final gene set
    adata, adata_sc = filter_genes(adata, adata_sc)

    adata = construct_spatial_graph(adata, config)

    adata = add_contrastive_labels(adata)

    adata, adata_sc = train_graphst(adata, adata_sc, config)

    adata = project_cells(adata, adata_sc, config)

    adata = identify_domains(adata, config)

    output_path = save_results(adata, sample_name, config)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(1)

    config_path = sys.argv[1]
    sample_name = sys.argv[2]

    if not os.path.exists(config_path):
        sys.exit(1)

    main(config_path, sample_name)
