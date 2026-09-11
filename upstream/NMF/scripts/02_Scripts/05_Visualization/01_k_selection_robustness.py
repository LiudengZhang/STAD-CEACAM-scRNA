#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""01_k_selection_robustness.py - Diagnostic plot for program filtering"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
CONFIG_PATH = BASE_DIR / "00_Config" / "config.yaml"
ROBUST_DIR = BASE_DIR / "03_Output" / "02_Robust_Programs"
OUTPUT_DIR = BASE_DIR / "03_Output" / "05_Figures"

with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

with open(ROBUST_DIR / "filtering_summary.yaml", 'r') as f:
    filtering_summary = yaml.safe_load(f)

all_programs = pd.read_csv(ROBUST_DIR / "all_nmf_programs.csv")
robust_programs = pd.read_csv(ROBUST_DIR / "robust_programs.csv")

fig, axes = plt.subplots(1, 3, figsize=(14, 5))

# Panel A: Program counts at each stage
stages = ['Initial\nNMF', 'Within-\nSample', 'Cross-\nSample', 'Final\nRobust']
counts = [d['count'] for d in filtering_summary['filtering_pipeline']]
colors = plt.cm.Blues(np.linspace(0.3, 0.9, 4))
bars = axes[0].bar(stages, counts, color=colors, edgecolor='black')
for bar, count in zip(bars, counts):
    axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, str(count), ha='center', fontsize=11, fontweight='bold')
axes[0].set_ylabel('Number of Programs', fontsize=12)
axes[0].set_title('A. Program Filtering Pipeline', fontsize=12, fontweight='bold')

# Panel B: Robust programs per sample
sample_counts = robust_programs.groupby('sample_id').size().sort_values(ascending=False)
axes[1].bar(range(len(sample_counts)), sample_counts.values, color='steelblue', edgecolor='black')
axes[1].set_xticks(range(len(sample_counts)))
axes[1].set_xticklabels(sample_counts.index, rotation=45, ha='right', fontsize=9)
axes[1].set_ylabel('Number of Robust Programs', fontsize=12)
axes[1].set_title('B. Robust Programs per Sample', fontsize=12, fontweight='bold')

# Panel C: Retention rate by K
k_values = config['nmf']['k_values']
initial_by_k = all_programs.groupby('k').size()
robust_by_k = robust_programs.groupby('k').size()
retention_rates = [(robust_by_k.get(k, 0) / initial_by_k.get(k, 1) * 100) for k in k_values]
axes[2].bar(k_values, retention_rates, color='coral', edgecolor='black')
axes[2].set_xlabel('K (NMF rank)', fontsize=12)
axes[2].set_ylabel('Retention Rate (%)', fontsize=12)
axes[2].set_title('C. Program Retention by K', fontsize=12, fontweight='bold')
axes[2].set_xticks(k_values)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "01_k_robustness_diagnostic.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print("Figure 1 complete")
