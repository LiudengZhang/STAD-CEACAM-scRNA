#!/usr/bin/env python3
"""
Panel G: MP4/5 Top Gene Loading Dotplot (rotated — all 5 MPs)
- Columns: Top 8 genes from MP4 and MP5
- Rows: All 5 MPs (MP1–MP5) to show specificity
- Dot size: loading score
- Color: red intensity by loading score
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter

# Central config
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

BASE_DIR = Path(__file__).parent
INTERMEDIATE = NMF_INTERMEDIATE
STOMACH_NMF_DIR = NMF_PER_SAMPLE

DOT_COLOR = '#C62828'

# Nature Cancer 4× scaling
DPI = 300
SCALE = 4
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 9.0 * SCALE    # wide horizontal layout
PANEL_HEIGHT_CM = 3.0 * SCALE   # short — 5 rows of MPs


def get_program_genes(program_id, nmf_dir):
    parts = program_id.rsplit('_', 2)
    sample, k, p_num = parts[0], int(parts[1][1:]), int(parts[2][1:])
    genes_file = nmf_dir / sample / f'{sample}_nmf_k{k}_genes.csv'
    if genes_file.exists():
        df = pd.read_csv(genes_file)
        return list(df[df.columns[p_num - 1]].dropna().tolist()[:50])
    return []


def get_mp_consensus_genes(assignments, mp_id, nmf_dir, top_n=50):
    programs = assignments[assignments['metaprogram_id'] == mp_id]['program_id'].tolist()
    gene_counts = Counter()
    for pid in programs:
        for i, g in enumerate(get_program_genes(pid, nmf_dir)):
            gene_counts[g] += (50 - i)
    return [g for g, c in gene_counts.most_common(top_n)], gene_counts


def create_panel_G():
    print("Creating Panel G: Gene Loading Dotplot (rotated, all 5 MPs)...")

    stomach_assign = pd.read_csv(INTERMEDIATE / 'panel_A2_stomach_all_mp_assignments.csv')

    all_mps = ['MP1', 'MP2', 'MP3', 'MP4', 'MP5']
    mp_genes = {}
    mp_counts = {}
    for mp in all_mps:
        genes, counts = get_mp_consensus_genes(stomach_assign, mp, STOMACH_NMF_DIR)
        mp_genes[f'S-{mp}'] = genes
        mp_counts[f'S-{mp}'] = counts
        print(f"  S-{mp}: {len(genes)} genes")

    # Select top 8 genes from MP4 and MP5
    top_n = 8
    mp4_top = mp_genes['S-MP4'][:top_n]
    mp5_top = [g for g in mp_genes['S-MP5'] if g not in mp4_top][:top_n]

    # Combine: MP4 genes first, then MP5 genes (left to right)
    all_genes = mp4_top + mp5_top
    n_mp4 = len(mp4_top)
    print(f"  Total genes: {len(all_genes)} (MP4: {n_mp4}, MP5: {len(mp5_top)})")

    # Build score matrix: rows=MPs, cols=genes (rotated from before)
    mp_names = [f'S-{mp}' for mp in all_mps]
    score_matrix = np.zeros((len(mp_names), len(all_genes)))

    for i, mp_name in enumerate(mp_names):
        counts = mp_counts[mp_name]
        max_count = max(counts.values()) if counts else 1
        for j, gene in enumerate(all_genes):
            score_matrix[i, j] = counts.get(gene, 0) / max_count

    # Create figure
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    fig, ax = plt.subplots(figsize=(PANEL_WIDTH_CM * CM_TO_INCH, PANEL_HEIGHT_CM * CM_TO_INCH))

    # Dotplot: x=genes, y=MPs
    for i, mp_name in enumerate(mp_names):
        for j, gene in enumerate(all_genes):
            score = score_matrix[i, j]
            if score > 0:
                size = (score * 200 + 20) * SCALE
                alpha = 0.3 + 0.7 * score
                ax.scatter(j, i, s=size, c=DOT_COLOR, alpha=alpha,
                          edgecolors='black', linewidth=0.5)

    # Styling
    ax.set_xticks(range(len(all_genes)))
    ax.set_xticklabels(all_genes, fontsize=6 * SCALE, style='italic', rotation=45, ha='right')
    ax.set_yticks(range(len(mp_names)))
    ax.set_yticklabels(all_mps, fontsize=6 * SCALE, fontweight='normal')

    ax.set_xlim(-0.5, len(all_genes) - 0.5)
    ax.set_ylim(-0.5, len(mp_names) - 0.5)

    # Vertical dashed line separating MP4 genes from MP5 genes
    ax.axvline(x=n_mp4 - 0.5, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.tick_params(axis='both', width=0.5, length=4)

    # Legend for size — upper right corner
    legend_sizes = [0.25, 0.5, 0.75, 1.0]
    legend_labels = ['25%', '50%', '75%', '100%']
    legend_elements = []
    for sz, lab in zip(legend_sizes, legend_labels):
        legend_elements.append(plt.scatter([], [], s=(sz*200+20)*SCALE, c=DOT_COLOR, alpha=0.6,
                                          edgecolors='black', linewidth=0.5, label=lab))

    ax.legend(handles=legend_elements, title='Loading', loc='upper left',
              bbox_to_anchor=(1.02, 1.0), fontsize=5 * SCALE, title_fontsize=6 * SCALE,
              frameon=True, ncol=1)

    # Extra left margin so y-axis labels don't sit at SVG edge
    # (prevents panel label "J" from colliding with "MP5" in assembly)
    plt.tight_layout(rect=[0.03, 0, 1, 1])

    output = BASE_DIR / 'mp45_gene_dotplot.png'
    plt.savefig(output, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.svg'), dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(BASE_DIR / 'mp45_gene_dotplot.pdf', format='pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output}")


if __name__ == '__main__':
    create_panel_G()
