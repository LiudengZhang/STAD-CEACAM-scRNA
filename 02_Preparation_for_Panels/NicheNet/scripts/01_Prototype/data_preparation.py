#!/usr/bin/env python3
"""
Data Preparation Module for NicheNet Analysis

Prepares expression matrices for NicheNet R analysis.

Author: Standardized Pipeline Team
Date: October 2025
"""

import os
import logging
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite
from utils import format_cell_type_name


class NicheNetDataPreparator:
    """Prepare data for NicheNet analysis."""

    def __init__(self, config):
        self.config = config
        self.h5ad_file = config['data']['h5ad_file']
        self.sender = config['cell_types']['sender']
        self.receivers = config['cell_types']['receivers']
        self.comparisons = config['comparisons']
        self.output_dir = os.path.join(
            config['output']['base_dir'],
            config['output'].get('prepared_data_dir', '03_Results/prepared_data')
        )

        os.makedirs(self.output_dir, exist_ok=True)

    def load_data(self):
        """Load single-cell data."""
        logging.info(f"Loading data from: {self.h5ad_file}")
        adata = sc.read_h5ad(self.h5ad_file)

        if hasattr(adata, 'raw') and adata.raw is not None:
            adata = adata.raw.to_adata()

        logging.info(f"Data loaded: {adata.n_obs} cells × {adata.n_vars} genes")
        return adata

    def prepare_nichenet_format(self, adata, comparison_name):
        """
        Prepare MTX files for NicheNet analysis by extracting cell type subsets.

        Parameters:
        -----------
        adata : AnnData
            Filtered dataset
        comparison_name : str
            Name of comparison

        Returns:
        --------
        str : Path to saved data directory
        """
        logging.info(f"  Preparing NicheNet format for: {comparison_name}")

        comparison_dir = os.path.join(self.output_dir, comparison_name)
        os.makedirs(comparison_dir, exist_ok=True)

        # Get cell type column
        cell_type_col = self.config['data'].get('cell_type_column', 'cell_type')

        # Create metadata file mapping cell types to MTX directories
        import yaml
        from scipy import sparse
        mtx_paths = {}

        # For each cell type, extract and save as MTX
        for cell_type in [self.sender] + self.receivers:
            if cell_type not in adata.obs[cell_type_col].values:
                logging.warning(f"    Cell type '{cell_type}' not found in data")
                continue

            # Subset cells for this cell type
            mask = adata.obs[cell_type_col] == cell_type
            subset = adata[mask].copy()

            if subset.n_obs == 0:
                logging.warning(f"    No cells found for {cell_type}")
                continue

            # Create MTX directory for this cell type
            safe_name = cell_type.replace(' ', '_').replace('+', 'plus').replace('-', '_')
            mtx_dir = os.path.join(comparison_dir, f"mtx_{safe_name}")
            os.makedirs(mtx_dir, exist_ok=True)

            # Get expression matrix (use raw if available)
            if hasattr(subset, 'raw') and subset.raw is not None:
                expr_matrix = subset.raw.X
                gene_names = subset.raw.var_names
            else:
                expr_matrix = subset.X
                gene_names = subset.var_names

            # Convert to sparse if dense
            if not sparse.issparse(expr_matrix):
                expr_matrix = sparse.csr_matrix(expr_matrix)

            # Save as MTX format (genes x cells for Seurat)
            mmwrite(os.path.join(mtx_dir, 'matrix.mtx'), expr_matrix.T)

            # Save features (genes)
            with open(os.path.join(mtx_dir, 'features.tsv'), 'w') as f:
                for gene in gene_names:
                    f.write(f"{gene}\n")

            # Save barcodes (cell IDs)
            with open(os.path.join(mtx_dir, 'barcodes.tsv'), 'w') as f:
                for barcode in subset.obs_names:
                    f.write(f"{barcode}\n")

            mtx_paths[cell_type] = mtx_dir
            logging.info(f"    {cell_type}: {subset.n_obs} cells saved to {mtx_dir}")

        # Save metadata
        metadata_file = os.path.join(comparison_dir, 'mtx_paths.yaml')
        with open(metadata_file, 'w') as f:
            yaml.dump(mtx_paths, f, default_flow_style=False)

        logging.info(f"    Saved MTX paths metadata: {metadata_file}")

        return comparison_dir

    def prepare_all_data(self):
        """Prepare data for all comparisons."""
        logging.info("Starting data preparation for NicheNet analysis")

        adata = self.load_data()

        # Apply top-level filters (site, phase, response) FIRST
        top_level_filters = self.config.get('filters', {})
        if top_level_filters:
            logging.info(f"Applying top-level filters: {top_level_filters}")
            filtered_base = adata.copy()

            # Apply site filter
            site_col = top_level_filters.get('site_column')
            site_val = top_level_filters.get('site')
            if site_col and site_val and site_col in filtered_base.obs.columns:
                before_count = filtered_base.n_obs
                filtered_base = filtered_base[filtered_base.obs[site_col] == site_val].copy()
                logging.info(f"  Site filter ({site_col}={site_val}): {before_count} -> {filtered_base.n_obs} cells")

            # Apply phase filter (Pre, Post, or None for combined)
            phase_col = top_level_filters.get('phase_column')
            phase_val = top_level_filters.get('phase')
            if phase_col and phase_val and phase_col in filtered_base.obs.columns:
                before_count = filtered_base.n_obs
                filtered_base = filtered_base[filtered_base.obs[phase_col] == phase_val].copy()
                logging.info(f"  Phase filter ({phase_col}={phase_val}): {before_count} -> {filtered_base.n_obs} cells")
            elif phase_val is None:
                logging.info(f"  Phase filter: None (using all phases - combined analysis)")

            # Handle unified response column for combined analysis
            response_col = top_level_filters.get('response_column')
            if response_col == 'response_unified':
                # Create unified response column from pre/post groupings
                logging.info("  Creating unified response column for combined analysis")
                filtered_base.obs['response_unified'] = 'Not selected'
                pre_mask = filtered_base.obs['Treatment phase'] == 'Pre'
                post_mask = filtered_base.obs['Treatment phase'] == 'Post'
                filtered_base.obs.loc[pre_mask, 'response_unified'] = filtered_base.obs.loc[pre_mask, 'stomach_pre_grouping']
                filtered_base.obs.loc[post_mask, 'response_unified'] = filtered_base.obs.loc[post_mask, 'stomach_post_grouping']

            # Apply response filter
            response_val = top_level_filters.get('response')
            if response_col and response_val and response_col in filtered_base.obs.columns:
                before_count = filtered_base.n_obs
                filtered_base = filtered_base[filtered_base.obs[response_col] == response_val].copy()
                logging.info(f"  Response filter ({response_col}={response_val}): {before_count} -> {filtered_base.n_obs} cells")
        else:
            filtered_base = adata.copy()

        prepared_data = {}

        for comparison in self.comparisons:
            comparison_name = comparison['name']
            logging.info(f"\nPreparing comparison: {comparison_name}")

            # Apply comparison-level filters on top of base filtered data
            filters = comparison.get('filters', {})
            if filters:
                filtered_adata = filtered_base.copy()
                for col, val in filters.items():
                    if col in filtered_adata.obs.columns:
                        filtered_adata = filtered_adata[filtered_adata.obs[col] == val].copy()
            else:
                filtered_adata = filtered_base.copy()

            # Prepare data
            comparison_dir = self.prepare_nichenet_format(filtered_adata, comparison_name)

            prepared_data[comparison_name] = {
                'data_dir': comparison_dir,
                'n_cells': filtered_adata.n_obs
            }

        logging.info(f"\nData preparation completed for {len(prepared_data)} comparisons")
        return prepared_data
