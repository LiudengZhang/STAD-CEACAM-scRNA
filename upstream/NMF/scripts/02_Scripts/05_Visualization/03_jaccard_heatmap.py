#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""03_jaccard_heatmap.py - Jaccard similarity heatmap"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
MP_DIR = BASE_DIR / "03_Output" / "03_MetaPrograms"
OUTPUT_DIR = BASE_DIR / "03_Output" / "05_Figures"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

sim_matrix = np.load(MP_DIR / "jaccard_similarity_matrix.npy")
with open(MP_DIR / "program_order.txt", 'r') as f:
    program_ids = [line.strip() for line in f]

assignments_file = MP_DIR / "metaprogram_assignments.csv"
program_to_mp = {}
if assignments_file.exists():
    assignments_df = pd.read_csv(assignments_file)
    program_to_mp = dict(zip(assignments_df['program_id'], assignments_df['metaprogram_id']))

distance_matrix = 1 - sim_matrix
np.fill_diagonal(distance_matrix, 0)
Z = linkage(squareform(distance_matrix), method='average')
order = dendrogram(Z, no_plot=True)['leaves']
sim_ordered = sim_matrix[np.ix_(order, order)]
programs_ordered = [program_ids[i] for i in order]

mp_ids = list(set(program_to_mp.values()))
mp_colors = plt.cm.tab20(np.linspace(0, 1, len(mp_ids)))
mp_color_map = dict(zip(mp_ids, mp_colors))

row_colors = []
for pid in programs_ordered:
    mp = program_to_mp.get(pid, 'Unassigned')
    row_colors.append(mp_color_map.get(mp, [0.8, 0.8, 0.8, 1.0]))

fig = plt.figure(figsize=(12, 10))
ax_main = fig.add_axes([0.15, 0.15, 0.7, 0.7])
im = ax_main.imshow(sim_ordered, cmap=sns.color_palette("RdYlBu_r", as_cmap=True), aspect='auto', vmin=0, vmax=1)
ax_main.set_xticks([])
ax_main.set_yticks([])
ax_main.set_xlabel('NMF Programs (clustered)', fontsize=12)
ax_main.set_ylabel('NMF Programs (clustered)', fontsize=12)

cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax).set_label('Jaccard Similarity', fontsize=11)

ax_top = fig.add_axes([0.15, 0.86, 0.7, 0.02])
for i, color in enumerate(row_colors):
    ax_top.axvspan(i, i+1, color=color)
ax_top.set_xlim(0, len(row_colors))
ax_top.set_xticks([])
ax_top.set_yticks([])

fig.suptitle('Jaccard Similarity of Robust NMF Programs', fontsize=14, fontweight='bold', y=0.95)
plt.savefig(OUTPUT_DIR / "03_jaccard_similarity_heatmap.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print("Figure 3 complete")
