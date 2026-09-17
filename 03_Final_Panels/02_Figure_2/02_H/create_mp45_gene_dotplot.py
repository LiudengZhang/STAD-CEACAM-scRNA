#!/usr/bin/env python3
"""
Figure 2 panel G - top gene loadings of the stomach metaprograms MP4 and MP5,
against all five metaprograms.

- Columns: the top 8 genes of MP4 and the top 8 of MP5
- Rows: all five metaprograms, so that specificity is visible
- Dot area and colour intensity: the consensus loading score

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 2 G       (PROVENANCE.csv; NOT inferred from "02_H")

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the legend - at 5 * SCALE. MARK carries the non-type
    point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The dot areas take AREA; the dot edge width and the MP4/MP5 divider rule
    take MARK. Tick widths and lengths and spine widths do not: those are
    style, and cnsplots sets them. The legend frame is kept, because it is not
    type.

LEGEND KEYS
    cnsplots sets `legend.markerscale = 0.5`, which a scatter handler applies
    as an *area* factor of 0.25. The earlier drawing had no such reduction, so
    carrying the handle areas over unchanged would draw each size key at a
    quarter of the dot it stands for and the key would stop meaning what it
    says. Each key is therefore divided by `markerscale ** 2`, so that after
    the legend applies its factor the key prints at exactly the area of the dot
    it labels.

Every value read, every gene, every consensus score and every string is the
earlier drawing's. The drawing code is the same code.
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
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402
import slots                                             # noqa: E402
from cnsfig import legend as cnslegend                    # noqa: E402

BASE_DIR = Path(__file__).parent
INTERMEDIATE = NMF_INTERMEDIATE                           # noqa: F405
STOMACH_NMF_DIR = NMF_PER_SAMPLE                          # noqa: F405

DOT_COLOR = '#C62828'
#: The largest dot on the published page, 3.27 mm across (counted off
#: 'Figure 2.pdf' geometry, 2026-09-14), as scatter's s = diameter^2 in pt.
DOT_S_MAX = (3.27 * style.PT_PER_MM) ** 2
#: The width of the size-key column at the right of the matrix.
KEY_COLUMN_MM = 12.0

SCALE = 4                           # the earlier canvas multiplier
SMALL_PT = 5.0                      # the earlier smallest body type

PANEL_LETTER = "G"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(2, PANEL_LETTER)


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

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    # THE DOTS ARE THE PUBLISHED PAGE'S SIZE  (2026-09-14)
    #   Counted off 'Figure 2.pdf' geometry: 70 circles from 1.00 to 3.27 mm
    #   across, the largest being a loading of 1.0. scatter's s is the
    #   diameter squared in points, so the largest is DOT_S_MAX and every
    #   other dot is the same (score * 200 + 20) / 220 fraction of it the
    #   earlier drawing used. Loading and alpha are untouched.
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

    # Styling
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

    # THE SIZE KEY IS THE SHARED COMPACT ONE  (2026-09-14)
    #   The framed ax.legend that was here took a column at the right with
    #   one circle per line and paper around each; the author read it as
    #   waste, and asked that this panel and Figure 4E decode their dots the
    #   same way. cnsfig.legend.compact_key lays the four circles out on a
    #   pitch taken from their own radii in a KEY_COLUMN_MM column flush with
    #   the matrix, under the heading the page prints. The circles are drawn
    #   from dot_s, so the key is the dots' own scale.
    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0,
                     pos.width - (KEY_COLUMN_MM + 1.0) / PANEL_W_MM,
                     pos.height])
    key_areas = [dot_s(f) for f in (0.25, 0.5, 0.75, 1.0)]
    key = cnslegend.compact_key(
        fig, size_areas=key_areas, size_labels=['25%', '50%', '75%', '100%'],
        column_mm=KEY_COLUMN_MM, right_mm=0.3,
        top_mm=(1.0 - pos.y1) * PANEL_H_MM, size_title='Loading',
        label_pt=style.tick_pt(), title_pt=style.tick_pt(),
        colour=DOT_COLOR)
    cnslegend.require_size_key(key, dot_areas=drawn, panel='Figure 2 panel G')

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, BASE_DIR / 'mp45_gene_dotplot')
    print(f"  Saved: {BASE_DIR / 'mp45_gene_dotplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    create_panel_G()
