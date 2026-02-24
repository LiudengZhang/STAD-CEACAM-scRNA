#!/usr/bin/env python3
"""
Step 1: Prepare scRNA-seq reference for BayesPrism deconvolution.
Uses full_dataset.h5ad .raw.X (log1p-normalized) → expm1 back-transform.
BayesPrism will use input.type = "GEP" (gene expression profile).

Exports:
  - sc_counts.tsv: cells × genes normalized expression matrix
  - sc_cell_types.tsv: cell_id, cell_type, cell_state
  - bulk_counts.tsv: samples × genes (TIGER tumor samples)
  - bulk_metadata.tsv: sample_id, response_NR
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.sparse import issparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from paths import FULL_DATASET_H5AD, TIGER_EXPR, TIGER_META

FULL_H5AD = FULL_DATASET_H5AD
OUTPUT_DIR = Path(__file__).parent

CELLS_PER_TYPE = 500
N_TOP_GENES = 5000
SEED = 42


def main():
    print("=" * 60)
    print("Step 1: Prepare scRNA-seq reference (GEP mode)")
    print("=" * 60)

    # --- Load bulk data first to get gene list ---
    print("\nLoading TIGER bulk expression...")
    bulk_df = pd.read_csv(TIGER_EXPR, sep='\t', index_col=0)
    bulk_genes = set(bulk_df.index)
    print(f"Bulk: {bulk_df.shape[0]:,} genes × {bulk_df.shape[1]} samples")

    # --- Load scRNA-seq ---
    print("\nLoading full dataset...")
    adata = sc.read_h5ad(FULL_H5AD)
    print(f"Full: {adata.shape[0]:,} cells × {adata.shape[1]} genes")

    # Filter to stomach
    stomach_mask = adata.obs['Sample site'] == 'Stomach'
    adata = adata[stomach_mask].copy()
    print(f"Stomach: {adata.shape[0]:,} cells")

    # Use .raw (57K genes, log1p-normalized)
    raw_var_names = list(adata.raw.var_names)
    print(f".raw: {adata.raw.shape[1]:,} genes")

    # Gene intersection with bulk
    common_genes = sorted(set(raw_var_names) & bulk_genes)
    print(f"Common with bulk: {len(common_genes):,}")

    # Remove MT/RP genes
    common_genes = [g for g in common_genes
                    if not g.startswith('MT-')
                    and not g.startswith('RPS')
                    and not g.startswith('RPL')]
    print(f"After removing MT/RP: {len(common_genes):,}")

    # Select top variable genes from bulk
    bulk_var = bulk_df.loc[common_genes].var(axis=1).sort_values(ascending=False)
    force_genes = {g for g in common_genes if g.startswith('CEACAM') or g == 'CD274'}
    top_genes = set(bulk_var.head(N_TOP_GENES).index) | force_genes
    final_genes = sorted(top_genes & set(common_genes))
    print(f"Final gene set: {len(final_genes):,}")
    ceacam_in = [g for g in final_genes if 'CEACAM' in g]
    print(f"  CEACAM genes: {ceacam_in}")

    # Subsample per cell type
    np.random.seed(SEED)
    keep_idx = []
    for ct in adata.obs['major_cell_type'].unique():
        ct_idx = np.where(adata.obs['major_cell_type'] == ct)[0]
        n_sample = min(CELLS_PER_TYPE, len(ct_idx))
        sampled = np.random.choice(ct_idx, n_sample, replace=False)
        keep_idx.extend(sampled)
        print(f"  {ct}: {len(ct_idx):,} → {n_sample}")

    keep_idx = sorted(keep_idx)
    adata_sub = adata[keep_idx]
    print(f"\nSubsampled: {adata_sub.shape[0]:,} cells")

    # Extract from .raw and back-transform: expm1(log1p(x)) = x
    gene_idx = [raw_var_names.index(g) for g in final_genes]
    raw_X = adata_sub.raw.X[:, gene_idx]
    if issparse(raw_X):
        raw_X = raw_X.toarray()

    # Back-transform from log1p
    expr = np.expm1(raw_X)
    # Replace NaN with 0
    expr = np.nan_to_num(expr, nan=0.0)
    print(f"Expression range: {expr.min():.2f} - {expr.max():.2f}")
    print(f"Mean: {expr.mean():.4f}")

    # Create DataFrames
    sc_counts_df = pd.DataFrame(expr, index=adata_sub.obs_names, columns=final_genes)

    # Cell labels
    state_col = 'minor_cell_state' if 'minor_cell_state' in adata_sub.obs.columns else 'major_cell_type'
    labels_df = pd.DataFrame({
        'cell_id': adata_sub.obs_names,
        'cell_type': adata_sub.obs['major_cell_type'].values,
        'cell_state': adata_sub.obs[state_col].values,
    })

    # Save
    print("\nSaving scRNA-seq reference...")
    sc_counts_df.to_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t')
    labels_df.to_csv(OUTPUT_DIR / "sc_cell_types.tsv", sep='\t', index=False)
    print(f"  sc_counts.tsv: {sc_counts_df.shape}")
    print(f"  sc_cell_types.tsv: {labels_df['cell_type'].nunique()} types")

    # --- Bulk data ---
    print("\nPreparing TIGER bulk data...")
    meta_df = pd.read_csv(TIGER_META, sep='\t')
    tumor_meta = meta_df[meta_df['Treatment'] != 'Normal'].copy()
    tumor_samples = [s for s in tumor_meta['sample_id'] if s in bulk_df.columns]

    bulk_subset = bulk_df.loc[final_genes, tumor_samples].T
    bulk_subset.to_csv(OUTPUT_DIR / "bulk_counts.tsv", sep='\t')
    print(f"  bulk_counts.tsv: {bulk_subset.shape}")

    tumor_meta_out = tumor_meta[tumor_meta['sample_id'].isin(tumor_samples)][
        ['sample_id', 'response_NR', 'response', 'TumorPurity']
    ].copy()
    tumor_meta_out.to_csv(OUTPUT_DIR / "bulk_metadata.tsv", sep='\t', index=False)
    print(f"  R: {(tumor_meta_out['response_NR'] == 'R').sum()}, N: {(tumor_meta_out['response_NR'] == 'N').sum()}")

    print("\nStep 1 complete!")


if __name__ == '__main__':
    main()
