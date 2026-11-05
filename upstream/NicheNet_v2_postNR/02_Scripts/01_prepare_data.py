#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Data Preparation for NicheNet: C3_Mac_Inflam_IL1B as Sender
"""

import yaml
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite
from scipy.sparse import csr_matrix, issparse
from pathlib import Path


def main():

    script_dir = Path(__file__).parent
    app_dir = script_dir.parent
    config_path = app_dir / "01_Config" / "config.yaml"
    output_dir = app_dir / "prepared_data"
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    adata = sc.read_h5ad(config["data"]["reference_h5ad"])

    filter_col = list(config["filters"].keys())[0]
    filter_val = config["filters"][filter_col]
    mask_filter = adata.obs[filter_col] == filter_val
    adata = adata[mask_filter].copy()
    raw_adata = adata.raw.to_adata()
    # Filter to genes expressed in >1% of cells
    if issparse(raw_adata.X):
        frac_expressed = np.array((raw_adata.X > 0).mean(axis=0)).flatten()
    else:
        frac_expressed = np.mean(raw_adata.X > 0, axis=0)
    gene_mask = frac_expressed > 0.01
    raw_adata = raw_adata[:, gene_mask].copy()

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

    ct_counts = raw_adata.obs[cell_type_col].value_counts()
    for ct, count in ct_counts.items():
        print(f"  {ct}: {count:,}")

    # Show sender info
    # np.asarray is required, do not remove it. This mask indexes a scipy
    # sparse matrix below; scipy's sparse __getitem__ calls .nonzero() on a
    # boolean index, and pandas removed Series.nonzero() in 1.0, so passing
    # the Series raises AttributeError under current pandas/scipy. The mask
    # is already aligned to raw_adata.obs, so coercing it to a numpy bool
    # array is a dtype change only and selects exactly the same rows.
    sender_mask = np.asarray(raw_adata.obs[sender_col] == sender_name)
    mtx_paths = {}

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

    for receiver in config["receivers"]:
        safe_name = (
            receiver.replace(" ", "_")
            .replace("+", "plus")
            .replace("/", "_")
            .replace("-", "_")
        )
        cell_dir = output_dir / f"mtx_{safe_name}"
        cell_dir.mkdir(exist_ok=True)

        mask = np.asarray(raw_adata.obs[cell_type_col] == receiver)
        n_cells = mask.sum()

        if n_cells == 0:
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

    mtx_paths_file = output_dir / "mtx_paths.yaml"
    with open(mtx_paths_file, "w") as f:
        yaml.dump(mtx_paths, f, default_flow_style=False)


if __name__ == "__main__":
    main()
