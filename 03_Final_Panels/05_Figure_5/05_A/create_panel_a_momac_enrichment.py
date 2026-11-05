#!/usr/bin/env python3
"""Figure 5A: Hallmark GSEA in MoMac cells (Post-NR vs Post-R, MAST prerank),
the nine gene sets with the largest absolute NES as horizontal bars.

Input:  GSEA_DIR/MoMac_mast_prerank_gsea.csv (00_Config/paths.py)
Output: momac_pathway_enrichment.{svg,pdf,png} in this directory
Run:    python create_panel_a_momac_enrichment.py
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *
from cnsfig import Panel, style  # noqa: E402

PANEL = Panel(5, "A", __file__, small_pt=6.0)

COLOR_POS = '#2166AC'   # NR-upregulated (positive NES) — blue
COLOR_NEG = '#B2182B'   # R-upregulated (negative NES) — red

GSEA_CSV = GSEA_DIR / "MoMac_mast_prerank_gsea.csv"

#: Display names with literal Greek letters.
PRINTED_NAMES = {
    "TNF-alpha Signaling via NF-kB": "TNF-α Signaling via NF-κB",
}


def main():
    PANEL.apply()

    results = pd.read_csv(GSEA_CSV)
    results['NES'] = pd.to_numeric(results['NES'], errors='coerce')
    results['abs_nes'] = results['NES'].abs()
    results = results.sort_values('abs_nes', ascending=False)

    # Top 9 by |NES|; ascending so positive NES is on top
    top9 = results.head(9).copy()
    top9 = top9.sort_values('NES', ascending=True)
    top9['clean_name'] = (top9['Term'].str.replace('HALLMARK_', '')
                          .str.replace('_', ' ').replace(PRINTED_NAMES))

    colors = [COLOR_POS if nes > 0 else COLOR_NEG for nes in top9['NES']]

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)

    y_pos = np.arange(len(top9))
    ax.barh(y_pos, top9['NES'], color=colors, edgecolor='black',
            linewidth=style.EDGE_PT)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top9['clean_name'])
    ax.set_xlabel('NES (Post-NR/Post-R)')

    ax.axvline(x=0, color='black', linewidth=style.RULE_PT)

    max_abs = max(abs(top9['NES'].min()), abs(top9['NES'].max()))
    ax.set_xlim(-max_abs - 0.3, max_abs + 0.3)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    PANEL.fit(fig)
    PANEL.finish(fig, 'momac_pathway_enrichment')


if __name__ == "__main__":
    main()
