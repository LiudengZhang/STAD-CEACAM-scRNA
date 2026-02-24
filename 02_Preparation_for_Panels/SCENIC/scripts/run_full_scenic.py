#!/usr/bin/env python3
"""
Full MoMac SCENIC Pipeline - All 68,289 cells
Run all steps: prepare -> GRN -> cisTarget -> AUCell
"""

import scanpy as sc
import pandas as pd
import numpy as np
import loompy as lp
from pathlib import Path
import sys
import time
from datetime import datetime

# pySCENIC imports
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD, SCENIC_RESULTS_DIR, SCENIC_DB_DIR

# Paths
DATA_PATH = MOMAC_H5AD
RESULTS_DIR = SCENIC_RESULTS_DIR

# Database paths (download from https://resources.aertslab.org/cistarget/)
DB_DIR = SCENIC_DB_DIR
RANKING_DBS = list(DB_DIR.glob('*.feather'))
MOTIF_ANNOTATIONS = DB_DIR / 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
TF_LIST = DB_DIR / 'allTFs_hg38.txt'


def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")
    sys.stdout.flush()


def step1_prepare():
    """Prepare loom file from h5ad"""
    log("=" * 60)
    log("STEP 1: Preparing data")
    log("=" * 60)

    adata = sc.read_h5ad(DATA_PATH)
    log(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # Get expression matrix (use raw if available)
    if adata.raw is not None:
        log("Using .raw.X")
        expr_matrix = adata.raw.X
        gene_names = adata.raw.var_names.tolist()
    else:
        log("Using .X")
        expr_matrix = adata.X
        gene_names = adata.var_names.tolist()

    if hasattr(expr_matrix, 'toarray'):
        expr_matrix = expr_matrix.toarray()

    # Create expression dataframe
    expr_df = pd.DataFrame(expr_matrix, index=adata.obs_names, columns=gene_names)
    log(f"Expression matrix: {expr_df.shape}")

    # Save metadata
    metadata_cols = ['minor_cell_state', 'major_cell_type', 'stomach_pre_grouping',
                     'stomach_post_grouping', 'stomach_treatment', 'liver_pre_grouping']
    metadata_cols = [c for c in metadata_cols if c in adata.obs.columns]
    metadata_df = adata.obs[metadata_cols].copy()

    if 'X_umap' in adata.obsm:
        metadata_df['UMAP_1'] = adata.obsm['X_umap'][:, 0]
        metadata_df['UMAP_2'] = adata.obsm['X_umap'][:, 1]

    metadata_df['cell_id'] = adata.obs_names
    metadata_df.to_csv(RESULTS_DIR / 'cell_metadata.csv', index=False)
    log(f"Saved metadata: {RESULTS_DIR / 'cell_metadata.csv'}")

    # Save gene names
    with open(RESULTS_DIR / 'gene_names.txt', 'w') as f:
        f.write('\n'.join(gene_names))

    return expr_df, gene_names


def step2_grn_inference(expr_df):
    """Run GRNBoost2 for GRN inference"""
    log("=" * 60)
    log("STEP 2: GRN Inference (GRNBoost2)")
    log("=" * 60)

    # Load TF list
    tf_names = pd.read_csv(TF_LIST, header=None)[0].tolist()
    tf_names = [tf for tf in tf_names if tf in expr_df.columns]
    log(f"TFs in expression data: {len(tf_names)}")

    start = time.time()
    adjacencies = grnboost2(expr_df, tf_names=tf_names, verbose=True)
    elapsed = time.time() - start
    log(f"GRNBoost2 completed in {elapsed/3600:.2f} hours")

    adjacencies.to_csv(RESULTS_DIR / 'adjacencies.tsv', sep='\t', index=False)
    log(f"Saved: {RESULTS_DIR / 'adjacencies.tsv'}")

    return adjacencies, tf_names


def step3_cistarget(expr_df, adjacencies, tf_names):
    """Run cisTarget for motif enrichment"""
    log("=" * 60)
    log("STEP 3: cisTarget Motif Enrichment")
    log("=" * 60)

    # Load ranking databases
    dbs = [RankingDatabase(fname=str(db), name=db.stem) for db in RANKING_DBS]
    log(f"Loaded {len(dbs)} ranking databases")

    # Create modules from adjacencies
    modules = list(modules_from_adjacencies(adjacencies, expr_df))
    log(f"Created {len(modules)} modules")

    start = time.time()
    df_motifs = prune2df(dbs, modules, str(MOTIF_ANNOTATIONS))
    elapsed = time.time() - start
    log(f"cisTarget completed in {elapsed/3600:.2f} hours")

    df_motifs.to_csv(RESULTS_DIR / 'motif_enrichment.csv', index=False)

    # Convert to regulons
    regulons = df2regulons(df_motifs)
    log(f"Created {len(regulons)} regulons")

    # Save regulons
    regulon_list = []
    for reg in regulons:
        for target, weight in reg.gene2weight.items():
            regulon_list.append({'TF': reg.name, 'target': target, 'importance': weight})

    regulon_df = pd.DataFrame(regulon_list)
    regulon_df.to_csv(RESULTS_DIR / 'regulons.csv', index=False)
    log(f"Saved: {RESULTS_DIR / 'regulons.csv'}")

    return regulons


def step4_aucell(expr_df, regulons):
    """Run AUCell for regulon activity scoring"""
    log("=" * 60)
    log("STEP 4: AUCell Scoring")
    log("=" * 60)

    start = time.time()
    auc_matrix = aucell(expr_df, regulons, num_workers=8)
    elapsed = time.time() - start
    log(f"AUCell completed in {elapsed/60:.2f} minutes")

    auc_matrix.to_csv(RESULTS_DIR / 'aucell_matrix.csv')
    log(f"Saved: {RESULTS_DIR / 'aucell_matrix.csv'}")

    # Also save as h5ad for easy loading
    import anndata
    adata_auc = anndata.AnnData(X=auc_matrix.values,
                                 obs=pd.DataFrame(index=auc_matrix.index),
                                 var=pd.DataFrame(index=auc_matrix.columns))
    adata_auc.write_h5ad(RESULTS_DIR / 'aucell_matrix.h5ad')
    log(f"Saved: {RESULTS_DIR / 'aucell_matrix.h5ad'}")

    return auc_matrix


def main():
    log("=" * 60)
    log("FULL MoMac SCENIC PIPELINE")
    log(f"Start time: {datetime.now()}")
    log("=" * 60)

    total_start = time.time()

    # Step 1: Prepare
    expr_df, gene_names = step1_prepare()

    # Step 2: GRN
    adjacencies, tf_names = step2_grn_inference(expr_df)

    # Step 3: cisTarget
    regulons = step3_cistarget(expr_df, adjacencies, tf_names)

    # Step 4: AUCell
    auc_matrix = step4_aucell(expr_df, regulons)

    total_elapsed = time.time() - total_start
    log("=" * 60)
    log(f"PIPELINE COMPLETE")
    log(f"Total runtime: {total_elapsed/3600:.2f} hours")
    log(f"End time: {datetime.now()}")
    log("=" * 60)


if __name__ == "__main__":
    main()
