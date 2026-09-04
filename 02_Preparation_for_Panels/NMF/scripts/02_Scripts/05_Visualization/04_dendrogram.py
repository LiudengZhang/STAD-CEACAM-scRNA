#!/usr/bin/env python3
"""04_dendrogram.py - Meta-program dendrogram"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
MP_DIR = BASE_DIR / "03_Output" / "03_MetaPrograms"
OUTPUT_DIR = BASE_DIR / "03_Output" / "05_Figures"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

Z = np.load(MP_DIR / "hierarchical_linkage.npy")
with open(MP_DIR / "program_order.txt", 'r') as f:
    program_ids = [line.strip() for line in f]

assignments_file = MP_DIR / "metaprogram_assignments.csv"
program_to_mp = {}
if assignments_file.exists():
    assignments_df = pd.read_csv(assignments_file)
    program_to_mp = dict(zip(assignments_df['program_id'], assignments_df['metaprogram_id']))

mp_ids = sorted(set(program_to_mp.values()))
colors = plt.cm.tab20(np.linspace(0, 1, max(len(mp_ids), 1)))
mp_color_map = {mp: colors[i] for i, mp in enumerate(mp_ids)}

fig, ax = plt.subplots(figsize=(14, 8))
dend = dendrogram(Z, labels=program_ids, leaf_rotation=90, leaf_font_size=6, ax=ax, color_threshold=0, above_threshold_color='gray')

for lbl in ax.get_xmajorticklabels():
    program = lbl.get_text()
    mp = program_to_mp.get(program, 'Unassigned')
    if mp in mp_color_map:
        c = mp_color_map[mp]
        lbl.set_color('#%02x%02x%02x' % (int(c[0]*255), int(c[1]*255), int(c[2]*255)))

ax.set_xlabel('NMF Programs', fontsize=12)
ax.set_ylabel('Distance (1 - Jaccard)', fontsize=12)
ax.set_title('Hierarchical Clustering of Robust NMF Programs', fontsize=14, fontweight='bold')

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=mp_color_map[mp], label=mp) for mp in mp_ids]
ax.legend(handles=legend_elements, loc='upper right', fontsize=9, title='Meta-Programs', ncol=2)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "04_metaprogram_dendrogram.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print("Figure 4 complete")
