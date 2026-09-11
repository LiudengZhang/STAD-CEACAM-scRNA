#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
01_prepare_per_sample_data.py
Prepare per-sample expression matrices for 3CA NMF analysis.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from pathlib import Path
from datetime import datetime

def log(msg):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
CONFIG_PATH = BASE_DIR / "00_Config" / "config.yaml"
INPUT_PATH = BASE_DIR / "01_Input_Data" / "epithelial_with_raw_counts_full.h5ad"
OUTPUT_DIR = BASE_DIR / "03_Output" / "01_Per_Sample_NMF"

with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

log("=" * 80)
log("PHASE 1.1: PREPARE PER-SAMPLE DATA FOR 3CA NMF")
log("=" * 80)

log(f"\n[1/4] Loading data from: {INPUT_PATH}")
adata = sc.read_h5ad(INPUT_PATH)
log(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

sample_col = None
for col in ['sample_id', 'Sample', 'sample', 'SampleID']:
    if col in adata.obs.columns:
        sample_col = col
        break

if sample_col is None:
    raise ValueError("Sample column not found")

log(f"Using sample column: {sample_col}")
samples = adata.obs[sample_col].unique()
log(f"Found {len(samples)} unique samples")

log(f"\n[2/4] Computing highly variable genes...")
n_hvg = config['data']['n_hvg']
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor='seurat_v3', subset=False)
hvg_genes = adata.var_names[adata.var['highly_variable']].tolist()
log(f"Selected {len(hvg_genes)} highly variable genes")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
hvg_file = OUTPUT_DIR / "hvg_genes.txt"
with open(hvg_file, 'w') as f:
    for gene in hvg_genes:
        f.write(f"{gene}\n")

log(f"\n[3/4] Splitting data by sample...")
min_cells = config['data']['min_cells_per_sample']
sample_info = []

for sample in samples:
    sample_mask = adata.obs[sample_col] == sample
    n_cells = sample_mask.sum()
    if n_cells < min_cells:
        log(f"  Skipping {sample}: only {n_cells} cells")
        continue

    adata_sample = adata[sample_mask, hvg_genes].copy()
    if hasattr(adata_sample.X, 'toarray'):
        counts = adata_sample.X.toarray()
    else:
        counts = adata_sample.X.copy()

    df = pd.DataFrame(counts, index=adata_sample.obs_names, columns=hvg_genes)
    sample_dir = OUTPUT_DIR / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    counts_file = sample_dir / f"{sample}_counts.csv"
    df.to_csv(counts_file)
    meta_file = sample_dir / f"{sample}_metadata.csv"
    adata_sample.obs.to_csv(meta_file)
    log(f"  Saved {sample}: {n_cells} cells x {len(hvg_genes)} genes")
    sample_info.append({'sample_id': sample, 'n_cells': n_cells, 'n_genes': len(hvg_genes)})

summary_df = pd.DataFrame(sample_info)
summary_file = OUTPUT_DIR / "sample_summary.csv"
summary_df.to_csv(summary_file, index=False)

log("\n" + "=" * 80)
log("PHASE 1.1 COMPLETE")
log(f"Prepared {len(sample_info)} samples")
log("=" * 80)
