#!/usr/bin/env python3
"""02_cluster_metaprograms.py - Cluster programs into meta-programs"""

import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
CONFIG_PATH = BASE_DIR / "00_Config" / "config.yaml"
INPUT_DIR = BASE_DIR / "03_Output" / "03_MetaPrograms"
ROBUST_DIR = BASE_DIR / "03_Output" / "02_Robust_Programs"

with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

jaccard_threshold = config['clustering']['jaccard_threshold']
min_programs = config['clustering']['min_programs_per_mp']
linkage_method = config['clustering']['linkage_method']

log("PHASE 3.2: CLUSTER META-PROGRAMS")

sim_matrix = np.load(INPUT_DIR / "jaccard_similarity_matrix.npy")
with open(INPUT_DIR / "program_order.txt", 'r') as f:
    program_ids = [line.strip() for line in f]

n_programs = len(program_ids)
remaining = set(range(n_programs))
meta_programs = []
mp_id = 0

while remaining:
    remaining_list = list(remaining)
    if len(remaining_list) == 1:
        idx = remaining_list[0]
        meta_programs.append({'mp_id': f"MP{mp_id + 1}", 'program_indices': [idx], 'program_ids': [program_ids[idx]]})
        remaining.remove(idx)
        mp_id += 1
        continue

    avg_sim = [np.mean([sim_matrix[i, j] for j in remaining_list if j != i]) for i in remaining_list]
    seed_idx = remaining_list[np.argmax(avg_sim)]
    cluster = [seed_idx] + [i for i in remaining_list if i != seed_idx and sim_matrix[seed_idx, i] >= jaccard_threshold]

    for idx in cluster:
        remaining.remove(idx)

    meta_programs.append({'mp_id': f"MP{mp_id + 1}", 'program_indices': cluster, 'program_ids': [program_ids[i] for i in cluster]})
    mp_id += 1

valid_mps = [mp for mp in meta_programs if len(mp['program_indices']) >= min_programs]

distance_matrix = 1 - sim_matrix
np.fill_diagonal(distance_matrix, 0)
Z = linkage(squareform(distance_matrix), method=linkage_method)
np.save(INPUT_DIR / "hierarchical_linkage.npy", Z)

assignments = [{'program_id': pid, 'metaprogram_id': mp['mp_id'], 'mp_size': len(mp['program_ids'])} for mp in valid_mps for pid in mp['program_ids']]
pd.DataFrame(assignments).to_csv(INPUT_DIR / "metaprogram_assignments.csv", index=False)

mp_summary = [{'metaprogram_id': mp['mp_id'], 'n_programs': len(mp['program_ids']), 'program_ids': ','.join(mp['program_ids'])} for mp in valid_mps]
pd.DataFrame(mp_summary).to_csv(INPUT_DIR / "metaprogram_summary.csv", index=False)

log(f"Created {len(valid_mps)} meta-programs")
log("PHASE 3.2 COMPLETE")
