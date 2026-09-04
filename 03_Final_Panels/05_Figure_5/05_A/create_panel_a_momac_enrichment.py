#!/usr/bin/env python3
"""
Panel A: MoMac GSEA horizontal barplot - top 9 Hallmark gene sets by |NES|.
4x scaling method.

Reads the MAST prerank GSEA table, GSEA/MoMac_mast_prerank_gsea.csv.

That input is the whole history of this panel. The twelve per-cell-type
*_mast_prerank_gsea.csv were moved under GSEA/_archived/ after the figures were
made, and this script was repointed at GSEA/post/MoMac_gsea_hallmark.csv - a
different run, built on the doubly normalised .X (00_Data_Audit/FINDINGS.md
sections 1 and 7). That repoint, not the figure, is why Panel A stopped
reproducing, and twice led an agent to conclude the published panel was wrong.
Both conclusions were retracted. On 2026-09-01 the original table was restored
here from _archived/ (byte-identical, md5 0e005bfe87764018f8edbb8e994bc637) and
the script pointed back at it.

Against the archived table the nine bars reproduce the printed panel to within
0.0001 NES, measured off the vector rectangles in
00_GROUND_TRUTH/figures/Figure 5.pdf against its -2/0/2 tick centres:

    TNF-alpha Signaling via NF-kB  +2.140    Coagulation      -1.488
    Inflammatory Response          +1.860    Spermatogenesis  -1.803
    Interferon Gamma Response      +1.750    Mitotic Spindle  -1.875
    IL-6/JAK/STAT3 Signaling       +1.748    G2-M Checkpoint  -2.223
                                             E2F Targets      -2.450

Do not repoint this at GSEA/post/ or at 05_GSEA_Summary/gsea_data/gsea_momac.csv.
Neither reproduces the paper, and gsea_momac.csv swaps out four of the nine gene
sets that Figure 5B and the Results sentence both rest on.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *
from shared.figure_config import use_panel_style

BASE_DIR = Path(__file__).parent

# 4× scaling
DPI = 300
SCALE = 4
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 6.0 * SCALE
PANEL_HEIGHT_CM = 4.2 * SCALE

# Standard R/NR colors
COLOR_POS = '#B2182B'   # NR-upregulated (positive NES) — red
COLOR_NEG = '#2166AC'   # R-upregulated (negative NES) — blue

GSEA_CSV = GSEA_DIR / "MoMac_mast_prerank_gsea.csv"


def main():
    use_panel_style(font_pt=7)

    print("Loading MAST Prerank GSEA results...")
    results = pd.read_csv(GSEA_CSV)
    results['NES'] = pd.to_numeric(results['NES'], errors='coerce')
    results['abs_nes'] = results['NES'].abs()
    results = results.sort_values('abs_nes', ascending=False)
    print(f"  Loaded {len(results)} pathways")

    # Top 9, sorted for barh (ascending so positive ends up on top)
    top9 = results.head(9).copy()
    top9 = top9.sort_values('NES', ascending=True)
    top9['clean_name'] = top9['Term'].str.replace('HALLMARK_', '').str.replace('_', ' ')

    colors = [COLOR_POS if nes > 0 else COLOR_NEG for nes in top9['NES']]

    fig, ax = plt.subplots(figsize=(PANEL_WIDTH_CM * CM_TO_INCH, PANEL_HEIGHT_CM * CM_TO_INCH))

    y_pos = np.arange(len(top9))
    ax.barh(y_pos, top9['NES'], color=colors, edgecolor='black', linewidth=0.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top9['clean_name'], fontsize=6 * SCALE)
    # The printed panel names the contrast here - 'NES (Post-NR/Post-R)'. No script
    # on disk emits that string; it was added when the figure was assembled, and
    # writing it here instead would be inventing provenance for it.
    ax.set_xlabel('NES', fontsize=6 * SCALE)
    ax.tick_params(axis='x', labelsize=6 * SCALE, width=1.0, length=4)
    ax.tick_params(axis='y', width=1.0, length=4)

    ax.axvline(x=0, color='black', linewidth=0.8)

    max_abs = max(abs(top9['NES'].min()), abs(top9['NES'].max()))
    ax.set_xlim(-max_abs - 0.3, max_abs + 0.3)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)

    plt.tight_layout()

    output = BASE_DIR / 'momac_pathway_enrichment.png'
    plt.savefig(output, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.pdf'), format='pdf', dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output}")


if __name__ == "__main__":
    main()
