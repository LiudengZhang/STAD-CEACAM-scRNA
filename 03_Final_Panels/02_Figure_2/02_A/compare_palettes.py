#!/usr/bin/env python3
"""
Compare 8 color palettes for epithelial UMAP (9 clusters).
Generates a 2x4 grid for side-by-side comparison.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

BASE_DIR = Path(__file__).parent

# ── 8 palettes for 9 clusters ──────────────────────────────────────────

PALETTES = {
    "1. Set1": [
        '#E41A1C', '#377EB8', '#4DAF4A', '#FF7F00', '#984EA3',
        '#A65628', '#F781BF', '#66C2A5', '#FC8D62',
    ],
    "2. Set2": [
        '#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854',
        '#FFD92F', '#E5C494', '#B3B3B3', '#8DD3C7',
    ],
    "3. Paired": [
        '#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99',
        '#E31A1C', '#FDBF6F', '#FF7F00', '#CAB2D6',
    ],
    "4. Dark2": [
        '#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E',
        '#E6AB02', '#A6761D', '#666666', '#8DD3C7',
    ],
    "5. NPG (Nature)": [
        '#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F',
        '#8491B4', '#91D1C2', '#DC0000', '#7E6148',
    ],
    "6. Paul Tol": [
        '#332288', '#88CCEE', '#44AA99', '#117733', '#999933',
        '#DDCC77', '#CC6677', '#882255', '#AA4499',
    ],
    "7. Okabe-Ito": [
        '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
        '#D55E00', '#CC79A7', '#000000', '#999999',
    ],
    "8. dittoSeq": [
        '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
        '#D55E00', '#CC79A7', '#7A6A65', '#90AD1C',
    ],
}

SHORT_LABELS = [
    'C0_PTMA', 'C1_KRT19', 'C2_CEACAM5/6', 'C3_Chief_Like',
    'C4_MUC5AC', 'C5_Stem_TPX2', 'C6_Stem_SPINK4', 'C7_MT1E', 'C8_CD74',
]

ORIG_CATS = [
    'C0_Epi_PTMA', 'C1_Epi_KRT19', 'C2_Epi_CEACAM6',
    'C3_Epi_Chief_Like_PGC', 'C4_Epi_MUC5AC', 'C5_Epi_Stem_Like_TPX2',
    'C6_Epi_Stem_Like_SPINK4', 'C7_Epi_MT1E', 'C8_Epi_CD74',
]


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("Loading data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    print(f"Loaded {adata.n_obs} cells")

    umap = adata.obsm['X_umap']
    states = adata.obs['minor_cell_state'].values

    # Shuffle once for consistent ordering across panels
    idx = np.random.RandomState(42).permutation(len(states))
    umap_shuffled = umap[idx]
    states_shuffled = states[idx]

    # 2 rows x 4 cols
    fig, axes = plt.subplots(2, 4, figsize=(32, 18))
    axes = axes.flatten()

    for i, (name, colors) in enumerate(PALETTES.items()):
        ax = axes[i]

        for j, cat in enumerate(ORIG_CATS):
            mask = states_shuffled == cat
            ax.scatter(
                umap_shuffled[mask, 0],
                umap_shuffled[mask, 1],
                c=colors[j],
                s=2,
                alpha=0.6,
                edgecolors='none',
                rasterized=True,
                label=SHORT_LABELS[j],
            )

        ax.set_title(name, fontsize=16, fontweight='bold', pad=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Legend
        handles = [mpatches.Patch(color=colors[j], label=SHORT_LABELS[j])
                   for j in range(len(ORIG_CATS))]
        ax.legend(
            handles=handles,
            loc='center left',
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=8,
            handlelength=1,
            handletextpad=0.4,
            labelspacing=0.5,
        )

    plt.tight_layout()
    output = BASE_DIR / 'palette_comparison.png'
    plt.savefig(output, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output}")


if __name__ == '__main__':
    main()
