#!/usr/bin/env python3
"""01_compute_jaccard_similarity.py - Compute Jaccard similarity matrix"""

import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime

def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

def jaccard_similarity(genes1, genes2):
    set1, set2 = set(genes1), set(genes2)
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
INPUT_DIR = BASE_DIR / "03_Output" / "02_Robust_Programs"
OUTPUT_DIR = BASE_DIR / "03_Output" / "03_MetaPrograms"

log("PHASE 3.1: COMPUTE JACCARD SIMILARITY MATRIX")

robust_df = pd.read_csv(INPUT_DIR / "robust_programs.csv")
robust_df['gene_list'] = robust_df['genes'].apply(lambda x: x.split(','))
program_ids = robust_df['program_id'].tolist()
gene_lists = robust_df['gene_list'].tolist()
n_programs = len(program_ids)

similarity_matrix = np.zeros((n_programs, n_programs))
for i in range(n_programs):
    for j in range(i, n_programs):
        sim = 1.0 if i == j else jaccard_similarity(gene_lists[i], gene_lists[j])
        similarity_matrix[i, j] = sim
        similarity_matrix[j, i] = sim

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
pd.DataFrame(similarity_matrix, index=program_ids, columns=program_ids).to_csv(OUTPUT_DIR / "jaccard_similarity_matrix.csv")
np.save(OUTPUT_DIR / "jaccard_similarity_matrix.npy", similarity_matrix)
with open(OUTPUT_DIR / "program_order.txt", 'w') as f:
    for pid in program_ids:
        f.write(f"{pid}\n")

log(f"Computed {n_programs}x{n_programs} similarity matrix")
log("PHASE 3.1 COMPLETE")
