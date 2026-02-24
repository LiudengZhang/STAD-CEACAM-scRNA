#!/usr/bin/env python3
"""
Panel A: MoMac GSEA Horizontal Barplot — Top 9 Hallmark Pathways by |NES|
4× scaling method for Nature Cancer.
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

GSEA_CSV = GSEA_POST_DIR / "MoMac_gsea_hallmark.csv"


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

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
