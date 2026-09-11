#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""02_cross_sample_recurrence.py - Filter by cross-sample recurrence"""

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

cross_threshold = config['robust_filtering']['cross_sample_threshold']
log("PHASE 2.2: CROSS-SAMPLE RECURRENCE FILTERING")

programs_df = pd.read_csv(INPUT_DIR / "within_sample_filtered.csv")
programs_df['gene_list'] = programs_df['genes'].apply(lambda x: x.split(','))

cross_sample_counts = []
cross_sample_max_overlap = []

for idx, row in programs_df.iterrows():
    sample = row['sample_id']
    genes = row['gene_list']
    other_programs = programs_df[programs_df['sample_id'] != sample]
    max_overlap = 0
    overlap_count = 0
    for _, other_row in other_programs.iterrows():
        overlap = compute_overlap(genes, other_row['gene_list'])
        max_overlap = max(max_overlap, overlap)
        if overlap >= cross_threshold:
            overlap_count += 1
    cross_sample_counts.append(overlap_count)
    cross_sample_max_overlap.append(max_overlap)

programs_df['cross_sample_count'] = cross_sample_counts
programs_df['cross_sample_max_overlap'] = cross_sample_max_overlap

filtered_df = programs_df[programs_df['cross_sample_count'] > 0].copy()
filtered_df = filtered_df.drop(columns=['gene_list'])
filtered_df.to_csv(INPUT_DIR / "cross_sample_filtered.csv", index=False)
log(f"Passing: {len(filtered_df)}/{len(programs_df)}")
log("PHASE 2.2 COMPLETE")
