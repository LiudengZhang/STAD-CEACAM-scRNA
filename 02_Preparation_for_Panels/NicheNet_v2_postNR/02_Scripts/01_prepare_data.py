#!/usr/bin/env python3
"""
Data Preparation for NicheNet: C3_Mac_Inflam_IL1B as Sender

Extracts expression matrices from full_dataset.h5ad (.raw layer)
for the sender (C3_Mac) and 13 receivers using final Round 4.2 annotations.

Filters to stomach post-treatment non-responder samples.
Exports MTX format (genes x cells) for NicheNet R analysis.
"""

import os
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite
from scipy.sparse import csr_matrix, issparse
from pathlib import Path


def main():
    print("=" * 70)
    print("NicheNet Data Preparation: C3_Mac_Inflam_IL1B Sender")
    print("Filter: Post-treatment Non-Responder (stomach_post_grouping == No-response)")
    print("=" * 70)

    # Paths
    script_dir = Path(__file__).parent
    app_dir = script_dir.parent
    config_path = app_dir / "01_Config" / "config.yaml"
    output_dir = app_dir / "prepared_data"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load config
    with open(config_path) as f:
        config = yaml.safe_load(f)

    print(f"\nConfig loaded from: {config_path}")
    print(f"Reference h5ad: {config['data']['reference_h5ad']}")

    # Load data
    print("\nLoading reference data...")
    adata = sc.read_h5ad(config["data"]["reference_h5ad"])
    print(f"  Full shape: {adata.n_obs} cells x {adata.n_vars} genes")

    # Filter to stomach post-treatment non-responders
    filter_col = list(config["filters"].keys())[0]
    filter_val = config["filters"][filter_col]
    mask_filter = adata.obs[filter_col] == filter_val
    adata = adata[mask_filter].copy()
    print(f"  After filter ({filter_col} == {filter_val}): {adata.n_obs} cells")

    # Use .raw layer for expression — all expressed genes (not just HVGs).
    # NicheNet needs a broad gene space for proper AUROC calculation.
    # The Zenodo ligand_target_matrix has 688 targets; with only 3,092 HVGs,
    # the overlap is too small (166 genes). Using all expressed genes from .raw
    # gives ~13K+ genes with much better NicheNet target coverage.
    print("\nExtracting .raw expression (all expressed genes)...")
    raw_adata = adata.raw.to_adata()
    # Filter to genes expressed in >1% of cells
    if issparse(raw_adata.X):
        frac_expressed = np.array((raw_adata.X > 0).mean(axis=0)).flatten()
    else:
        frac_expressed = np.mean(raw_adata.X > 0, axis=0)
    gene_mask = frac_expressed > 0.01
    raw_adata = raw_adata[:, gene_mask].copy()
    print(f"  Genes in .raw: {adata.raw.n_vars}")
    print(f"  Genes expressed in >1% cells: {gene_mask.sum()}")
    print(f"  Final shape: {raw_adata.n_obs} cells x {raw_adata.n_vars} genes")

    # Replace NaN with 0 (if any)
    X = raw_adata.X
    if issparse(X):
        X.data = np.nan_to_num(X.data, nan=0.0)
    else:
        X = np.nan_to_num(X, nan=0.0)
        raw_adata.X = X

    var_names = raw_adata.var_names
    cell_type_col = config["data"]["cell_type_column"]
    sender_col = config["data"]["sender_column"]
    sender_name = config["sender"]

    # Show cell type distribution
    print(f"\nCell type distribution ({cell_type_col}):")
    ct_counts = raw_adata.obs[cell_type_col].value_counts()
    for ct, count in ct_counts.items():
        print(f"  {ct}: {count:,}")

    # Show sender info
    sender_mask = raw_adata.obs[sender_col] == sender_name
    print(f"\nSender ({sender_name}): {sender_mask.sum()} cells")

    # Export MTX files
    print("\nExporting expression matrices per cell type...")
    mtx_paths = {}

    # --- Export sender ---
    safe_sender = sender_name.replace(" ", "_").replace("+", "plus").replace("/", "_")
    sender_dir = output_dir / f"mtx_{safe_sender}"
    sender_dir.mkdir(exist_ok=True)

    sender_X = raw_adata.X[sender_mask, :]
    if not issparse(sender_X):
        sender_X = csr_matrix(sender_X)
    sender_X_t = sender_X.T.tocsr()

    mmwrite(str(sender_dir / "matrix.mtx"), sender_X_t)
    pd.DataFrame({"gene": var_names}).to_csv(
        sender_dir / "features.tsv", sep="\t", index=False, header=False
    )
    pd.DataFrame({"barcode": raw_adata.obs_names[sender_mask]}).to_csv(
        sender_dir / "barcodes.tsv", sep="\t", index=False, header=False
    )
    mtx_paths[sender_name] = str(sender_dir)
    print(f"  Sender {sender_name}: {sender_mask.sum()} cells -> {sender_dir}")

    # --- Export each receiver ---
    for receiver in config["receivers"]:
        safe_name = (
            receiver.replace(" ", "_")
            .replace("+", "plus")
            .replace("/", "_")
            .replace("-", "_")
        )
        cell_dir = output_dir / f"mtx_{safe_name}"
        cell_dir.mkdir(exist_ok=True)

        mask = raw_adata.obs[cell_type_col] == receiver
        n_cells = mask.sum()

        if n_cells == 0:
            print(f"  {receiver}: 0 cells - SKIPPED")
            continue

        cell_X = raw_adata.X[mask, :]
        if not issparse(cell_X):
            cell_X = csr_matrix(cell_X)
        cell_X_t = cell_X.T.tocsr()

        mmwrite(str(cell_dir / "matrix.mtx"), cell_X_t)
        pd.DataFrame({"gene": var_names}).to_csv(
            cell_dir / "features.tsv", sep="\t", index=False, header=False
        )
        pd.DataFrame({"barcode": raw_adata.obs_names[mask]}).to_csv(
            cell_dir / "barcodes.tsv", sep="\t", index=False, header=False
        )
        mtx_paths[receiver] = str(cell_dir)
        print(f"  {receiver}: {n_cells:,} cells -> {cell_dir}")

    # Save MTX paths mapping
    mtx_paths_file = output_dir / "mtx_paths.yaml"
    with open(mtx_paths_file, "w") as f:
        yaml.dump(mtx_paths, f, default_flow_style=False)
    print(f"\nMTX paths saved to: {mtx_paths_file}")

    # Summary
    print("\n" + "=" * 70)
    print("Data Preparation Summary")
    print("=" * 70)
    print(f"Filter: {filter_col} == {filter_val}")
    print(f"Total cells: {adata.n_obs}")
    print(f"Genes: {raw_adata.n_vars}")
    print(f"Sender: {sender_name} ({sender_mask.sum()} cells)")
    print(f"Receivers exported: {len(mtx_paths) - 1}")
    print(f"Output directory: {output_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
