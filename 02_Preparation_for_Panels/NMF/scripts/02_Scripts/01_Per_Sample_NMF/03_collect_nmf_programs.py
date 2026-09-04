#!/usr/bin/env python3
"""03_collect_nmf_programs.py - Collect all NMF programs"""

import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime

def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
CONFIG_PATH = BASE_DIR / "00_Config" / "config.yaml"
INPUT_DIR = BASE_DIR / "03_Output" / "01_Per_Sample_NMF"
OUTPUT_DIR = BASE_DIR / "03_Output" / "02_Robust_Programs"

with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

k_values = config['nmf']['k_values']
top_genes = config['nmf']['top_genes_per_program']

log("PHASE 1.3: COLLECT ALL NMF PROGRAMS")

sample_summary = pd.read_csv(INPUT_DIR / "sample_summary.csv")
samples = sample_summary['sample_id'].tolist()

all_programs = []
program_genes = {}

for sample in samples:
    sample_dir = INPUT_DIR / sample
    if not sample_dir.exists():
        continue
    for k in k_values:
        genes_file = sample_dir / f"{sample}_nmf_k{k}_genes.csv"
        W_file = sample_dir / f"{sample}_nmf_k{k}_W.csv"
        if not genes_file.exists():
            continue
        genes_df = pd.read_csv(genes_file)
        W_df = pd.read_csv(W_file, index_col=0)
        for i, col in enumerate(genes_df.columns):
            program_id = f"{sample}_k{k}_P{i+1}"
            genes = genes_df[col].dropna().tolist()
            program_col = f"Program_{i+1}"
            weights = W_df[program_col].to_dict() if program_col in W_df.columns else {g: 1.0 for g in genes}
            all_programs.append({'program_id': program_id, 'sample_id': sample, 'k': k, 'program_idx': i+1, 'n_genes': len(genes), 'genes': ','.join(genes)})
            program_genes[program_id] = {'genes': genes, 'weights': {g: weights.get(g, 0) for g in genes}}

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
programs_df = pd.DataFrame(all_programs)
programs_df.to_csv(OUTPUT_DIR / "all_nmf_programs.csv", index=False)

gene_lists_dir = OUTPUT_DIR / "gene_lists"
gene_lists_dir.mkdir(exist_ok=True)
for program_id, data in program_genes.items():
    with open(gene_lists_dir / f"{program_id}_genes.txt", 'w') as f:
        for gene in data['genes']:
            f.write(f"{gene}\n")
    pd.DataFrame({'gene': data['genes'], 'weight': [data['weights'].get(g, 0) for g in data['genes']]}).to_csv(gene_lists_dir / f"{program_id}_weights.csv", index=False)

log(f"Collected {len(all_programs)} NMF programs")
log("PHASE 1.3 COMPLETE")
