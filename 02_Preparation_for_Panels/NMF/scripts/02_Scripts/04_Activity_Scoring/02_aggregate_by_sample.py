#!/usr/bin/env python3
"""02_aggregate_by_sample.py - Aggregate cell-level scores to sample level"""

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
INPUT_DIR = BASE_DIR / "03_Output" / "04_Activity_Scores"

with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

log("PHASE 4.2: AGGREGATE SCORES BY SAMPLE")

scores_df = pd.read_csv(INPUT_DIR / "cell_activity_scores.csv", index_col=0)
metadata_cols = ['sample_id', 'response', 'treatment_phase']
mp_cols = [col for col in scores_df.columns if col not in metadata_cols]

sample_df = scores_df.groupby('sample_id')[mp_cols].mean()
sample_df['n_cells'] = scores_df.groupby('sample_id').size()

pre_responders = config['sample_groups']['pre_responders']
pre_non_responders = config['sample_groups']['pre_non_responders']
post_responders = config['sample_groups']['post_responders']
post_non_responders = config['sample_groups']['post_non_responders']

def assign_group(sample_id):
    if sample_id in pre_responders:
        return 'Pre_Responder'
    elif sample_id in pre_non_responders:
        return 'Pre_NonResponder'
    elif sample_id in post_responders:
        return 'Post_Responder'
    elif sample_id in post_non_responders:
        return 'Post_NonResponder'
    return 'Unknown'

sample_df['clinical_group'] = sample_df.index.map(assign_group)
sample_df['timepoint'] = sample_df['clinical_group'].apply(lambda x: 'Pre' if 'Pre' in x else ('Post' if 'Post' in x else 'Unknown'))
sample_df['response_status'] = sample_df['clinical_group'].apply(lambda x: 'Responder' if 'Responder' in x and 'Non' not in x else ('NonResponder' if 'NonResponder' in x else 'Unknown'))

sample_df.to_csv(INPUT_DIR / "sample_activity_scores.csv")
log(f"Aggregated to {len(sample_df)} samples")
log("PHASE 4.2 COMPLETE")
