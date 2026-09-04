#!/usr/bin/env python3
"""
Figure 2, printed panel G, RESTYLED (Version B) - MP4/MP5 top gene loading
dotplot (all 5 metaprograms).

- Columns: Top 8 genes from MP4 and MP5
- Rows: All 5 MPs (MP1-MP5) to show specificity
- Dot size: loading score
- Color: red intensity by loading score

Version A is
`03_Final_Panels/02_Figure_2/02_H/create_mp45_gene_dotplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every gene, every consensus score and every string is
Version A's. The drawing code is the same code.

  printed panel  Figure 2 G     (PROVENANCE.csv; NOT inferred from "02_H")
  printed rect   90.5 x 32.5 mm    (panel_rects.csv)
  Version B box  120.0 x 42.0 mm

MARK
    Version A drew at SCALE = 4 (9.0 x 3.0 cm x 4 = 360 x 120 mm) and its
    smallest body type is the legend at `5 * SCALE`. So

        SCALE = 4, SMALL_PT = 5
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA  = MARK ** 2 = 0.1225

    The dot areas `(score * 200 + 20) * SCALE` take AREA; the dot edge width
    and the MP4/MP5 divider rule take MARK.

    Version A's spine linewidths and `tick_params(width=..., length=...)` are
    NOT carried over: those are axes furniture, cnsplots has its own settings
    for them, and following the library rather than rescaling the old numbers
    is the standard-methods rule. `frameon=True` on the legend IS carried over -
    it is Version A's explicit choice and it is not type.

LEGEND KEYS
    cnsplots sets `legend.markerscale = 0.5`, which a scatter handler applies as
    an *area* factor of 0.25. Version A had no such reduction, so carrying the
    handle areas over unchanged would draw the size key at a quarter of the dot
    it stands for and the key would stop meaning what it says. Each key is
    therefore divided by `markerscale ** 2` so that, after the legend applies
    it, it prints at exactly the area of the dot it labels. Same correction, and
    same reason, as the stage-2 exemplar S8_F.

    The panel grew from 90.5 x 32.5 mm to 120 x 42 mm: at 7 pt the sixteen
    45-degree gene labels need about 5.9 mm of column pitch and 11 mm of depth
    below the axis, and the size legend needs 17 mm to the right of it.
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
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402

BASE_DIR = Path(__file__).parent
INTERMEDIATE = NMF_INTERMEDIATE                           # noqa: F405
STOMACH_NMF_DIR = NMF_PER_SAMPLE                          # noqa: F405

DOT_COLOR = '#C62828'

PRINTED_MM = (90.5, 32.5)           # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 120.0, 42.0
MARGIN = dict(left=9.0, right=17.0, top=2.5, bottom=11.5)

SCALE = 4                           # Version A's canvas multiplier
SMALL_PT = 5.0                      # Version A's smallest body type


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

    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    # Dotplot: x=genes, y=MPs
    for i, mp_name in enumerate(mp_names):
        for j, gene in enumerate(all_genes):
            score = score_matrix[i, j]
            if score > 0:
                size = (score * 200 + 20) * SCALE * AREA
                alpha = 0.3 + 0.7 * score
                ax.scatter(j, i, s=size, c=DOT_COLOR, alpha=alpha,
                           edgecolors='black', linewidth=0.5 * MARK)

    # Styling
    ax.set_xticks(range(len(all_genes)))
    ax.set_xticklabels(all_genes, style='italic', rotation=45, ha='right')
    ax.set_yticks(range(len(mp_names)))
    ax.set_yticklabels(all_mps)

    ax.set_xlim(-0.5, len(all_genes) - 0.5)
    ax.set_ylim(-0.5, len(mp_names) - 0.5)

    # Vertical dashed line separating MP4 genes from MP5 genes
    ax.axvline(x=n_mp4 - 0.5, color='gray', linestyle='--',
               linewidth=0.5 * MARK, alpha=0.5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend for size - upper right corner. See LEGEND KEYS in the header for
    # the markerscale correction.
    key = 1.0 / plt.rcParams["legend.markerscale"] ** 2
    legend_sizes = [0.25, 0.5, 0.75, 1.0]
    legend_labels = ['25%', '50%', '75%', '100%']
    legend_elements = []
    for sz, lab in zip(legend_sizes, legend_labels):
        legend_elements.append(
            plt.scatter([], [], s=(sz * 200 + 20) * SCALE * AREA * key,
                        c=DOT_COLOR, alpha=0.6, edgecolors='black',
                        linewidth=0.5 * MARK, label=lab))

    ax.legend(handles=legend_elements, title='Loading', loc='upper left',
              bbox_to_anchor=(1.02, 1.0), frameon=True, ncol=1)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    style.save_panel(fig, BASE_DIR / 'mp45_gene_dotplot')
    print(f"  Saved: {BASE_DIR / 'mp45_gene_dotplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    create_panel_G()
