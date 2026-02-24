#!/usr/bin/env python3
"""
S2_G: 5-Metaprogram Gene Signature Heatmap (flipped — genes on x, MPs on y)
Seaborn heatmap style matching Figure 4A.
Data source: NMF intermediate files.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import NMF_PER_SAMPLE, NMF_INTERMEDIATE

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

SCALE = 4
DPI = 300
CM = 1 / 2.54

# 5-color palette for metaprograms
MP_COLORS = {
    'MP1': '#E64B35',
    'MP2': '#4DBBD5',
    'MP3': '#00A087',
    'MP4': '#3C5488',
    'MP5': '#F39B7F',
}

TOP_N_GENES = 10
ALL_MPS = ['MP1', 'MP2', 'MP3', 'MP4', 'MP5']

OUT_DIR = Path(__file__).parent
STOMACH_NMF_DIR = NMF_PER_SAMPLE
INTERMEDIATE = NMF_INTERMEDIATE

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})


def get_program_genes(program_id, nmf_dir):
    """Get top genes for a given NMF program."""
    parts = program_id.rsplit('_', 2)
    sample, k, p_num = parts[0], int(parts[1][1:]), int(parts[2][1:])
    genes_file = nmf_dir / sample / f'{sample}_nmf_k{k}_genes.csv'
    if genes_file.exists():
        df = pd.read_csv(genes_file)
        return list(df[df.columns[p_num - 1]].dropna().tolist()[:50])
    return []


def get_mp_consensus_genes(assignments, mp_id, nmf_dir, top_n=50):
    """Get consensus genes for a metaprogram using weighted voting."""
    programs = assignments[assignments['metaprogram_id'] == mp_id]['program_id'].tolist()
    gene_counts = Counter()
    for pid in programs:
        for i, g in enumerate(get_program_genes(pid, nmf_dir)):
            gene_counts[g] += (50 - i)
    return [g for g, _ in gene_counts.most_common(top_n)], gene_counts


def main():
    print("S2_G: 5-Metaprogram Gene Signature Heatmap (flipped)")
    print("=" * 60)

    # Load MP assignments
    assign_file = INTERMEDIATE / 'panel_A2_stomach_all_mp_assignments.csv'
    print(f"Loading assignments: {assign_file}")
    assignments = pd.read_csv(assign_file)

    # Get consensus genes and scores for each MP
    mp_top_genes = {}
    mp_all_counts = {}
    for mp in ALL_MPS:
        genes, counts = get_mp_consensus_genes(assignments, mp, STOMACH_NMF_DIR)
        mp_top_genes[mp] = genes[:TOP_N_GENES]
        mp_all_counts[mp] = counts
        print(f"  {mp}: top {TOP_N_GENES} genes = {genes[:5]}...")

    # Build unique gene list (in MP order, no duplicates)
    all_genes = []
    gene_to_mp = {}
    for mp in ALL_MPS:
        for g in mp_top_genes[mp]:
            if g not in gene_to_mp:
                all_genes.append(g)
                gene_to_mp[g] = mp

    print(f"  Total unique genes: {len(all_genes)}")

    # Build score matrix: rows = MPs, cols = genes (FLIPPED from original)
    score_matrix = np.zeros((len(ALL_MPS), len(all_genes)))
    for i, mp in enumerate(ALL_MPS):
        counts = mp_all_counts[mp]
        max_count = max(counts.values()) if counts else 1
        for j, gene in enumerate(all_genes):
            score_matrix[i, j] = counts.get(gene, 0) / max_count

    score_df = pd.DataFrame(score_matrix, index=[f'S-{mp}' for mp in ALL_MPS], columns=all_genes)

    # Create heatmap — seaborn style matching Fig 4A
    fig_w = (0.3 * len(all_genes) + 3.0) * SCALE * CM
    fig_h = (0.5 * len(ALL_MPS) + 2.0) * SCALE * CM
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sns.heatmap(
        score_df,
        cmap='YlOrRd',
        vmin=0,
        vmax=1,
        ax=ax,
        xticklabels=True,
        yticklabels=True,
        linewidths=0.5,
        linecolor='lightgray',
        cbar_kws={'label': 'Normalized Loading', 'shrink': 0.6},
    )

    # De-rasterize (match Fig 4A)
    for coll in ax.collections:
        coll.set_rasterized(False)

    # Style x-axis gene labels
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=4.5 * SCALE, rotation=90,
                       ha='center', style='italic')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6 * SCALE, rotation=0)

    # Add vertical lines between MP gene blocks
    cumulative = 0
    for mp in ALL_MPS:
        block_size = sum(1 for g in all_genes if gene_to_mp[g] == mp)
        cumulative += block_size
        if cumulative < len(all_genes):
            ax.axvline(x=cumulative, color='black', linewidth=1.5)

    # Add MP color strips on the left
    for i, mp in enumerate(ALL_MPS):
        ax.add_patch(plt.Rectangle((-0.7, i), 0.2, 1,
                                    facecolor=MP_COLORS[mp], edgecolor='none',
                                    clip_on=False, transform=ax.transData))

    ax.set_title('Metaprogram Gene Signatures', fontsize=6 * SCALE, fontweight='normal')

    # Unified structural elements
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['left', 'bottom']:
        ax.spines[spine].set_linewidth(0.5 * SCALE)
    ax.tick_params(width=0.5 * SCALE, length=3 * SCALE)

    plt.tight_layout()

    stem = 'panel_S2_G'
    for ext in ['svg', 'png']:
        kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if ext == 'png':
            kw['dpi'] = DPI
        fig.savefig(OUT_DIR / f'{stem}.{ext}', **kw)
    print(f"Saved: {stem}")
    plt.close()


if __name__ == '__main__':
    main()
