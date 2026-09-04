#!/usr/bin/env python3
"""03_remove_redundancy.py - Remove redundant programs"""

import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime

def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

def compute_overlap(genes1, genes2):
    return len(set(genes1).intersection(set(genes2)))

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
CONFIG_PATH = BASE_DIR / "00_Config" / "config.yaml"
INPUT_DIR = BASE_DIR / "03_Output" / "02_Robust_Programs"

with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

redundancy_threshold = config['robust_filtering']['redundancy_threshold']
redundancy_genes = int(redundancy_threshold * 50)
log("PHASE 2.3: REDUNDANCY REMOVAL")

programs_df = pd.read_csv(INPUT_DIR / "cross_sample_filtered.csv")
programs_df['gene_list'] = programs_df['genes'].apply(lambda x: x.split(','))
samples = programs_df['sample_id'].unique()

robust_programs = []
for sample in samples:
    sample_programs = programs_df[programs_df['sample_id'] == sample].sort_values('cross_sample_count', ascending=False)
    selected = []
    for idx, row in sample_programs.iterrows():
        genes = row['gene_list']
        is_redundant = False
        for selected_genes in selected:
            if compute_overlap(genes, selected_genes) > redundancy_genes:
                is_redundant = True
                break
        if not is_redundant:
            selected.append(genes)
            robust_programs.append(row.to_dict())

robust_df = pd.DataFrame(robust_programs)
if 'gene_list' in robust_df.columns:
    robust_df = robust_df.drop(columns=['gene_list'])
robust_df.to_csv(INPUT_DIR / "robust_programs.csv", index=False)

gene_lists_dir = INPUT_DIR / "robust_gene_lists"
gene_lists_dir.mkdir(exist_ok=True)
for _, row in robust_df.iterrows():
    with open(gene_lists_dir / f"{row['program_id']}.txt", 'w') as f:
        for gene in row['genes'].split(','):
            f.write(f"{gene}\n")

all_programs = pd.read_csv(INPUT_DIR / "all_nmf_programs.csv")
within_filtered = pd.read_csv(INPUT_DIR / "within_sample_filtered.csv")
cross_filtered = pd.read_csv(INPUT_DIR / "cross_sample_filtered.csv")

summary = {
    'filtering_pipeline': [
        {'stage': 'Initial NMF programs', 'count': len(all_programs)},
        {'stage': 'After within-sample filter', 'count': len(within_filtered)},
        {'stage': 'After cross-sample filter', 'count': len(cross_filtered)},
        {'stage': 'Final robust programs', 'count': len(robust_df)}
    ]
}
with open(INPUT_DIR / "filtering_summary.yaml", 'w') as f:
    yaml.dump(summary, f)

log(f"Final robust programs: {len(robust_df)}")
log("PHASE 2.3 COMPLETE")
