#!/usr/bin/env python3
"""
Step 4: AUCell Regulon Activity Scoring - Full MoMac (12 Post-Stomach Samples)
=============================================================================
Purpose: Calculate per-cell activity scores for each regulon using AUCell

Input:
  - MoMac_12sample_scenic.loom (expression matrix)
  - regulons.pkl (from Step 3)

Output:
  - aucell_matrix.csv (regulon activity matrix)
  - aucell_matrix.h5ad (AnnData format)
  - aucell_binary.csv (binarized activity)

Runtime: ~10-30 minutes
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import pickle
import loompy as lp

def main():
    print("=" * 80)
    print("SCENIC Step 4: AUCell Regulon Activity Scoring - Full MoMac (12 samples)")
    print("=" * 80)
    print()

    # Central config
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
    from paths import SCENIC_RESULTS_DIR

    # Input files
    loom_file = SCENIC_RESULTS_DIR / 'MoMac_12sample_scenic.loom'
    regulons_file = SCENIC_RESULTS_DIR / 'regulons.pkl'

    # Output files
    output_matrix = SCENIC_RESULTS_DIR / 'aucell_matrix.csv'
    output_h5ad = SCENIC_RESULTS_DIR / 'aucell_matrix.h5ad'
    output_binary = SCENIC_RESULTS_DIR / 'aucell_binary.csv'

    # Check input files
    if not loom_file.exists():
        print(f"ERROR: Loom file not found: {loom_file}")
        sys.exit(1)

    if not regulons_file.exists():
        print(f"ERROR: Regulons file not found: {regulons_file}")
        print("Please run 03_run_cisTarget.py first")
        sys.exit(1)

    # Step 1: Load regulons
    print("[1/4] Loading regulons...")
    with open(regulons_file, 'rb') as f:
        regulons = pickle.load(f)
    print(f"   Loaded {len(regulons)} regulons")
    print()

    # Step 2: Load expression matrix from loom
    print("[2/4] Loading expression matrix from loom...")
    print(f"   Loom file: {loom_file}")

    with lp.connect(str(loom_file), mode='r', validate=False) as ds:
        ex_matrix = pd.DataFrame(
            data=ds[:, :],
            index=ds.ra.Gene,
            columns=ds.ca.CellID
        ).T  # Transpose so cells are rows
        print(f"   Loaded: {ex_matrix.shape[0]:,} cells x {ex_matrix.shape[1]:,} genes")
    print()

    # Step 3: Run AUCell
    print("[3/4] Running AUCell scoring...")
    print("   This calculates regulon activity scores per cell")
    print("   Runtime: ~10-30 minutes")
    print()

    from pyscenic.aucell import aucell

    # Run AUCell
    auc_mtx = aucell(ex_matrix, regulons, num_workers=8)

    print(f"   AUCell complete!")
    print(f"   Activity matrix shape: {auc_mtx.shape}")
    print()

    # Step 4: Save results
    print("[4/4] Saving results...")

    # Save as CSV
    auc_mtx.to_csv(output_matrix)
    print(f"   Saved: {output_matrix}")

    # Save as h5ad
    import anndata as ad
    adata = ad.AnnData(X=auc_mtx.values, obs=pd.DataFrame(index=auc_mtx.index), var=pd.DataFrame(index=auc_mtx.columns))
    adata.write(output_h5ad)
    print(f"   Saved: {output_h5ad}")

    # Binarize activity scores using threshold
    print()
    print("   Binarizing activity scores...")

    from pyscenic.binarization import binarize

    # Binarize with default AUCell thresholds
    binary_mtx, thresholds = binarize(auc_mtx)
    binary_mtx.to_csv(output_binary)
    print(f"   Saved: {output_binary}")
    print()

    # Print statistics
    print("   Top regulons by mean activity:")
    mean_activity = auc_mtx.mean().sort_values(ascending=False)
    for reg, activity in mean_activity.head(15).items():
        print(f"     {reg}: {activity:.4f}")
    print()

    # Key TF regulons
    key_tfs = ['FOSB', 'FOS', 'JUN', 'ATF3', 'NFKB1', 'NFKB2', 'RELA', 'ETS2',
               'NR4A1', 'NR4A2', 'EGR1', 'STAT1', 'STAT3', 'IRF1', 'CEBPB']

    print("   Key inflammatory TF regulon activities:")
    for tf in key_tfs:
        matching = [col for col in auc_mtx.columns if col.startswith(tf + '(')]
        if matching:
            for reg in matching:
                mean_act = auc_mtx[reg].mean()
                print(f"     {reg}: {mean_act:.4f}")
    print()

    # Summary
    print("=" * 80)
    print("AUCell Scoring Complete!")
    print("=" * 80)
    print(f"Results:")
    print(f"  - Cells: {auc_mtx.shape[0]:,}")
    print(f"  - Regulons: {auc_mtx.shape[1]}")
    print(f"  - Output files:")
    print(f"      {output_matrix}")
    print(f"      {output_h5ad}")
    print(f"      {output_binary}")
    print()
    print("SCENIC pipeline complete!")
    print("=" * 80)

if __name__ == "__main__":
    main()
