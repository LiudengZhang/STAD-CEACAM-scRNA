#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""01_within_sample_recurrence.py - Filter by within-sample recurrence"""

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

within_threshold = config['robust_filtering']['within_sample_threshold']
log("PHASE 2.1: WITHIN-SAMPLE RECURRENCE FILTERING")

programs_df = pd.read_csv(INPUT_DIR / "all_nmf_programs.csv")
programs_df['gene_list'] = programs_df['genes'].apply(lambda x: x.split(','))
samples = programs_df['sample_id'].unique()

passing_programs = []
for sample in samples:
    sample_programs = programs_df[programs_df['sample_id'] == sample].copy()
    for idx, row in sample_programs.iterrows():
        genes = row['gene_list']
        k = row['k']
        other_programs = sample_programs[sample_programs['k'] != k]
        max_overlap = 0
        recurrence_count = 0
        for _, other_row in other_programs.iterrows():
            overlap = compute_overlap(genes, other_row['gene_list'])
            max_overlap = max(max_overlap, overlap)
            if overlap >= within_threshold:
                recurrence_count += 1
        if recurrence_count > 0:
            passing_programs.append({**row.to_dict(), 'max_overlap': max_overlap, 'recurrence_count': recurrence_count})

filtered_df = pd.DataFrame(passing_programs)
if 'gene_list' in filtered_df.columns:
    filtered_df = filtered_df.drop(columns=['gene_list'])
filtered_df.to_csv(INPUT_DIR / "within_sample_filtered.csv", index=False)
log(f"Passing: {len(filtered_df)}/{len(programs_df)}")
log("PHASE 2.1 COMPLETE")
