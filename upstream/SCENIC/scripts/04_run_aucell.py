#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
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
"""

import pandas as pd
from pathlib import Path
import sys
import pickle
import loompy as lp

def main():

    base_dir = Path('/path/to/upstream-pipeline')
    project_dir = base_dir / '02_Playground/02_Cell-Cell_Interaction/01_Analysis/09_SCENIC/17_Full_MoMac_SCENIC'

    loom_file = project_dir / 'Results/MoMac_12sample_scenic.loom'
    regulons_file = project_dir / 'Results/regulons.pkl'

    output_matrix = project_dir / 'Results/aucell_matrix.csv'
    output_h5ad = project_dir / 'Results/aucell_matrix.h5ad'
    output_binary = project_dir / 'Results/aucell_binary.csv'

    if not loom_file.exists():
        sys.exit(1)

    if not regulons_file.exists():
        sys.exit(1)

    with open(regulons_file, 'rb') as f:
        regulons = pickle.load(f)

    with lp.connect(str(loom_file), mode='r', validate=False) as ds:
        ex_matrix = pd.DataFrame(
            data=ds[:, :],
            index=ds.ra.Gene,
            columns=ds.ca.CellID
        ).T  # Transpose so cells are rows

    from pyscenic.aucell import aucell

    auc_mtx = aucell(ex_matrix, regulons, num_workers=8)

    auc_mtx.to_csv(output_matrix)

    import anndata as ad
    adata = ad.AnnData(X=auc_mtx.values, obs=pd.DataFrame(index=auc_mtx.index), var=pd.DataFrame(index=auc_mtx.columns))
    adata.write(output_h5ad)

    from pyscenic.binarization import binarize

    # Binarize with default AUCell thresholds
    binary_mtx, thresholds = binarize(auc_mtx)
    binary_mtx.to_csv(output_binary)
    mean_activity = auc_mtx.mean().sort_values(ascending=False)
    for reg, activity in mean_activity.head(15).items():
        print(f"     {reg}: {activity:.4f}")

    key_tfs = ['FOSB', 'FOS', 'JUN', 'ATF3', 'NFKB1', 'NFKB2', 'RELA', 'ETS2',
               'NR4A1', 'NR4A2', 'EGR1', 'STAT1', 'STAT3', 'IRF1', 'CEBPB']

    for tf in key_tfs:
        matching = [col for col in auc_mtx.columns if col.startswith(tf + '(')]
        if matching:
            for reg in matching:
                mean_act = auc_mtx[reg].mean()

if __name__ == "__main__":
    main()
