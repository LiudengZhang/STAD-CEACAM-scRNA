#!/usr/bin/env python3
"""
Figure 2G: top consensus genes of the stomach metaprograms MP4 and MP5,
scored against all five metaprograms.

Columns are the top 8 genes of MP4 and of MP5; rows are MP1-MP5, so
specificity is visible. Dot area and opacity show the consensus loading,
scaled to each metaprogram's maximum. The consensus-gene function is shared
with 02_G/create_mp45_horizontal_boxplot.py.

Reads NMF_INTERMEDIATE/panel_A2_stomach_all_mp_assignments.csv and the
per-sample NMF gene tables under NMF_PER_SAMPLE (00_Config/paths.py).
Writes mp45_gene_dotplot.[svg|pdf|png] beside this script.

Run: python create_mp45_gene_dotplot.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "02_G"))
from paths import NMF_INTERMEDIATE, NMF_PER_SAMPLE         # noqa: E402
from cnsfig import Panel, style  # noqa: E402
from cnsfig import legend as cnslegend                    # noqa: E402
from create_mp45_horizontal_boxplot import get_mp_consensus_genes  # noqa: E402

INTERMEDIATE = NMF_INTERMEDIATE
STOMACH_NMF_DIR = NMF_PER_SAMPLE

DOT_COLOR = '#C62828'
# Largest dot: 3.27 mm across, as scatter's s (diameter^2 in pt).
DOT_S_MAX = (3.27 * style.PT_PER_MM) ** 2
# Width of the size-key column at the right of the matrix.
KEY_COLUMN_MM = 12.0

SCALE = 4                           # sets MARK relative to the tick size
SMALL_PT = 5.0

PANEL = Panel(2, "G", __file__, small_pt=SMALL_PT)


def create_panel_G():
    stomach_assign = pd.read_csv(INTERMEDIATE / 'panel_A2_stomach_all_mp_assignments.csv')

    all_mps = ['MP1', 'MP2', 'MP3', 'MP4', 'MP5']
    mp_genes = {}
    mp_counts = {}
    for mp in all_mps:
        genes, counts = get_mp_consensus_genes(stomach_assign, mp, STOMACH_NMF_DIR)
        mp_genes[f'S-{mp}'] = genes
        mp_counts[f'S-{mp}'] = counts

    top_n = 8
    mp4_top = mp_genes['S-MP4'][:top_n]
    mp5_top = [g for g in mp_genes['S-MP5'] if g not in mp4_top][:top_n]

    # Combine: MP4 genes first, then MP5 genes (left to right)
    all_genes = mp4_top + mp5_top
    n_mp4 = len(mp4_top)

    mp_names = [f'S-{mp}' for mp in all_mps]
    score_matrix = np.zeros((len(mp_names), len(all_genes)))

    for i, mp_name in enumerate(mp_names):
        counts = mp_counts[mp_name]
        max_count = max(counts.values()) if counts else 1
        for j, gene in enumerate(all_genes):
            score_matrix[i, j] = counts.get(gene, 0) / max_count

    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)

    # Dot area grows linearly with loading; a loading of 1.0 is DOT_S_MAX.
    def dot_s(score):
        return DOT_S_MAX * (score * 200 + 20) / 220.0

    # Dotplot: x=genes, y=MPs
    drawn = []
    for i, mp_name in enumerate(mp_names):
        for j, gene in enumerate(all_genes):
            score = score_matrix[i, j]
            if score > 0:
                size = dot_s(score)
                drawn.append(size)
                alpha = 0.3 + 0.7 * score
                ax.scatter(j, i, s=size, c=DOT_COLOR, alpha=alpha,
                           edgecolors='black', linewidth=style.EDGE_PT)

    ax.set_xticks(range(len(all_genes)))
    ax.set_xticklabels(all_genes, style='italic', rotation=45, ha='right')
    ax.set_yticks(range(len(mp_names)))
    ax.set_yticklabels(all_mps)

    ax.set_xlim(-0.5, len(all_genes) - 0.5)
    ax.set_ylim(-0.5, len(mp_names) - 0.5)

    # Vertical dashed line separating MP4 genes from MP5 genes
    ax.axvline(x=n_mp4 - 0.5, color='gray', linestyle='--',
               linewidth=style.RULE_PT, alpha=0.8)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Size key in a column to the right, drawn with the dots' own scale.
    PANEL.fit(fig)
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0,
                     pos.width - (KEY_COLUMN_MM + 1.0) / PANEL.w_mm,
                     pos.height])
    key_areas = [dot_s(f) for f in (0.25, 0.5, 0.75, 1.0)]
    key = cnslegend.compact_key(
        fig, size_areas=key_areas, size_labels=['25%', '50%', '75%', '100%'],
        column_mm=KEY_COLUMN_MM, right_mm=0.3,
        top_mm=(1.0 - pos.y1) * PANEL.h_mm, size_title='Loading',
        label_pt=style.tick_pt(), title_pt=style.tick_pt(),
        colour=DOT_COLOR)
    cnslegend.require_size_key(key, dot_areas=drawn, panel='Figure 2 panel G')

    PANEL.finish(fig, 'mp45_gene_dotplot')


if __name__ == '__main__':
    create_panel_G()
