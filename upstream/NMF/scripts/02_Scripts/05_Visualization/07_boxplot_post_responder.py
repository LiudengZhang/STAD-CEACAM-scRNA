#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""07_boxplot_post_responder.py - Post-treatment R vs NR boxplot"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
SCORES_DIR = BASE_DIR / "03_Output" / "04_Activity_Scores"
OUTPUT_DIR = BASE_DIR / "03_Output" / "05_Figures"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

scores_df = pd.read_csv(SCORES_DIR / "sample_activity_scores.csv", index_col=0)
metadata_cols = ['clinical_group', 'timepoint', 'response_status', 'n_cells', 'response', 'treatment_phase']
mp_cols = [col for col in scores_df.columns if col not in metadata_cols]

post_df = scores_df[scores_df['timepoint'] == 'Post'].copy()
responders = post_df[post_df['response_status'] == 'Responder']
non_responders = post_df[post_df['response_status'] == 'NonResponder']

n_mps = len(mp_cols)
n_cols = min(4, n_mps)
n_rows = (n_mps + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
axes = np.atleast_2d(axes)

for idx, mp in enumerate(sorted(mp_cols)):
    row, col = idx // n_cols, idx % n_cols
    ax = axes[row, col] if n_rows > 1 else axes[0, col]

    plot_data = [{'Response': r['response_status'], 'Activity': r[mp]} for _, r in post_df.iterrows()]
    plot_df = pd.DataFrame(plot_data)

    colors = {'Responder': '#2196F3', 'NonResponder': '#FF9800'}
    sns.boxplot(data=plot_df, x='Response', y='Activity', palette=colors, ax=ax, width=0.6)
    sns.stripplot(data=plot_df, x='Response', y='Activity', color='black', size=8, ax=ax, alpha=0.7)

    if len(responders) >= 2 and len(non_responders) >= 2:
        _, pval = mannwhitneyu(responders[mp].values, non_responders[mp].values, alternative='two-sided')
        sig_text = '***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else 'ns'))
        y_max = plot_df['Activity'].max()
        y_range = plot_df['Activity'].max() - plot_df['Activity'].min()
        ax.text(0.5, y_max + 0.1 * y_range, sig_text, ha='center', fontsize=14, fontweight='bold')
        ax.text(0.5, y_max + 0.2 * y_range, f'p={pval:.3f}', ha='center', fontsize=9)

    ax.set_title(mp.replace('_score', ''), fontsize=11, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Activity Score', fontsize=10)

for idx in range(n_mps, n_rows * n_cols):
    row, col = idx // n_cols, idx % n_cols
    (axes[row, col] if n_rows > 1 else axes[0, col]).axis('off')

plt.suptitle('Post-Treatment: Responder vs Non-Responder', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "07_boxplot_post_RvsNR.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print("Figure 7 complete")
