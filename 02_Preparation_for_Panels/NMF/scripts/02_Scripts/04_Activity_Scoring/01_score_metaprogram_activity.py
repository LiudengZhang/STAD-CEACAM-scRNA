#!/usr/bin/env python3
"""01_score_metaprogram_activity.py - Score meta-program activity per cell"""

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from pathlib import Path
from datetime import datetime

def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
CONFIG_PATH = BASE_DIR / "00_Config" / "config.yaml"
INPUT_H5AD = BASE_DIR / "01_Input_Data" / "epithelial_with_raw_counts_full.h5ad"
SIG_DIR = BASE_DIR / "03_Output" / "03_MetaPrograms" / "metaprogram_signatures"
OUTPUT_DIR = BASE_DIR / "03_Output" / "04_Activity_Scores"

with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

ctrl_size = config['scoring']['ctrl_size']
n_bins = config['scoring']['n_bins']

log("PHASE 4.1: SCORE META-PROGRAM ACTIVITY")

adata = sc.read_h5ad(INPUT_H5AD)
adata_score = adata.copy()

if hasattr(adata_score.X, 'toarray'):
    max_val = adata_score.X.toarray().max()
else:
    max_val = adata_score.X.max()

if max_val > 100:
    sc.pp.normalize_total(adata_score, target_sum=1e4)
    sc.pp.log1p(adata_score)

sig_files = list(SIG_DIR.glob("MP*_genes.txt"))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

activity_scores = {}
scored_mps = []

for sig_file in sorted(sig_files):
    mp_id = sig_file.stem.replace("_genes", "")
    with open(sig_file, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    genes_present = [g for g in genes if g in adata_score.var_names]
    if len(genes_present) < 10:
        continue

    score_name = f"{mp_id}_score"
    sc.tl.score_genes(adata_score, gene_list=genes_present, score_name=score_name, ctrl_size=ctrl_size, n_bins=n_bins)
    activity_scores[mp_id] = adata_score.obs[score_name].values
    scored_mps.append({'metaprogram_id': mp_id, 'n_genes_present': len(genes_present)})

scores_df = pd.DataFrame(activity_scores, index=adata_score.obs_names)
for col in ['sample_id', 'Sample', 'sample']:
    if col in adata_score.obs.columns:
        scores_df['sample_id'] = adata_score.obs[col].values
        break

scores_df.to_csv(OUTPUT_DIR / "cell_activity_scores.csv")

adata_out = adata.copy()
for mp_id in activity_scores:
    adata_out.obs[f"{mp_id}_score"] = activity_scores[mp_id]
adata_out.write_h5ad(OUTPUT_DIR / "epithelial_with_mp_scores.h5ad")

log(f"Scored {len(scored_mps)} meta-programs")
log("PHASE 4.1 COMPLETE")
