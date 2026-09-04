#!/usr/bin/env python3
"""03_extract_gene_signatures.py - Extract 50-gene signatures per meta-program"""

import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime
from collections import Counter

def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
CONFIG_PATH = BASE_DIR / "00_Config" / "config.yaml"
ROBUST_DIR = BASE_DIR / "03_Output" / "02_Robust_Programs"
INPUT_DIR = BASE_DIR / "03_Output" / "03_MetaPrograms"

with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

top_genes = config['nmf']['top_genes_per_program']
log("PHASE 3.3: EXTRACT META-PROGRAM GENE SIGNATURES")

assignments_df = pd.read_csv(INPUT_DIR / "metaprogram_assignments.csv")
robust_df = pd.read_csv(ROBUST_DIR / "robust_programs.csv")
robust_df['gene_list'] = robust_df['genes'].apply(lambda x: x.split(','))

gene_weights = {}
gene_lists_dir = ROBUST_DIR / "gene_lists"
for _, row in robust_df.iterrows():
    weight_file = gene_lists_dir / f"{row['program_id']}_weights.csv"
    if weight_file.exists():
        weight_df = pd.read_csv(weight_file)
        gene_weights[row['program_id']] = dict(zip(weight_df['gene'], weight_df['weight']))

meta_programs = assignments_df['metaprogram_id'].unique()
sig_dir = INPUT_DIR / "metaprogram_signatures"
sig_dir.mkdir(exist_ok=True)

mp_signatures = []
for mp_id in meta_programs:
    mp_programs = assignments_df[assignments_df['metaprogram_id'] == mp_id]['program_id'].tolist()
    gene_freq = Counter()
    gene_avg_weight = {}

    for program_id in mp_programs:
        program_row = robust_df[robust_df['program_id'] == program_id]
        if program_row.empty:
            continue
        genes = program_row.iloc[0]['gene_list']
        weights = gene_weights.get(program_id, {})
        for gene in genes:
            gene_freq[gene] += 1
            gene_avg_weight.setdefault(gene, []).append(weights.get(gene, 1.0))

    for gene in gene_avg_weight:
        gene_avg_weight[gene] = np.mean(gene_avg_weight[gene])

    gene_scores = [(gene, freq, gene_avg_weight.get(gene, 0)) for gene, freq in gene_freq.items()]
    gene_scores.sort(key=lambda x: (-x[1], -x[2]))
    top_gene_list = [g[0] for g in gene_scores[:top_genes]]

    with open(sig_dir / f"{mp_id}_genes.txt", 'w') as f:
        for gene in top_gene_list:
            f.write(f"{gene}\n")

    pd.DataFrame(gene_scores[:top_genes], columns=['gene', 'frequency', 'avg_weight']).to_csv(sig_dir / f"{mp_id}_gene_scores.csv", index=False)
    mp_signatures.append({'metaprogram_id': mp_id, 'n_programs': len(mp_programs), 'genes': ','.join(top_gene_list)})

pd.DataFrame(mp_signatures).to_csv(INPUT_DIR / "metaprogram_signatures_summary.csv", index=False)

all_genes = set()
for mp in mp_signatures:
    all_genes.update(mp['genes'].split(','))
sig_matrix = pd.DataFrame(0, index=sorted(all_genes), columns=[mp['metaprogram_id'] for mp in mp_signatures])
for mp in mp_signatures:
    for gene in mp['genes'].split(','):
        if gene in sig_matrix.index:
            sig_matrix.loc[gene, mp['metaprogram_id']] = 1
sig_matrix.to_csv(INPUT_DIR / "signature_gene_matrix.csv")

log(f"Created signatures for {len(meta_programs)} meta-programs")
log("PHASE 3.3 COMPLETE")
