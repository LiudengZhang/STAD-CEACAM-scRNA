#!/usr/bin/env python3
"""
Figure 5 panel A - MoMac GSEA horizontal barplot, top 9 Hallmark gene sets by
absolute NES.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 5 A       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

THE INPUT PATH IS NOT TOUCHED, AND MUST NOT BE.
-----------------------------------------------
This reads the MAST prerank GSEA table GSEA/MoMac_mast_prerank_gsea.csv. The
twelve per-cell-type *_mast_prerank_gsea.csv were moved under GSEA/_archived/
after the figures were made and this script was once repointed at
GSEA/post/MoMac_gsea_hallmark.csv - a different run, built on the doubly
normalised .X. That repoint, not the figure, is why the panel stopped
reproducing. Against the restored table the nine bars reproduce the printed
panel to within 0.0001 NES.

Do not repoint this at GSEA/post/ or at
05_GSEA_Summary/gsea_data/gsea_momac.csv. Neither reproduces the paper, and
gsea_momac.csv swaps out four of the nine gene sets that the panel and the
Results sentence both rest on.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the y tick labels, the x tick labels and the x axis
    label - at 6 * SCALE. MARK carries the non-type point sizes across to the
    1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The bar edge and the zero rule are scaled by it. Tick widths and lengths
    and spine widths are not: those are style, and cnsplots sets them.

The x axis label reads 'NES' and not the contrast string the printed panel
carries. No script on disk emits that string; it was added at assembly, and
writing it here would be inventing provenance for it.
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
import panel_style_cns as style
import slots

BASE_DIR = Path(__file__).parent

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 6.0                  # the earlier smallest body type, before * SCALE

PANEL_LETTER = "A"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(5, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(5, PANEL_LETTER)

# Standard R/NR colors
COLOR_POS = '#B2182B'   # NR-upregulated (positive NES) — red
COLOR_NEG = '#2166AC'   # R-upregulated (negative NES) — blue

GSEA_CSV = GSEA_DIR / "MoMac_mast_prerank_gsea.csv"

# How a gene-set name is printed where the table spells a Greek letter out.
# The table is not touched: the label is a display form, the row it names, its
# NES and its place in the ranking are unchanged. Every other NF-κB on this
# page is set in Greek - the four dotplots of panel N, the axis titles of K and
# M - so the page names one pathway one way. The two glyphs are set as literal
# characters and not as mathtext, which draws a symbol at 70% of the base size
# and would put them below the figure's 6 pt floor.
PRINTED_NAMES = {
    "TNF-alpha Signaling via NF-kB": "TNF-α Signaling via NF-κB",
}


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

    print("Loading MAST Prerank GSEA results...")
    results = pd.read_csv(GSEA_CSV)
    results['NES'] = pd.to_numeric(results['NES'], errors='coerce')
    results['abs_nes'] = results['NES'].abs()
    results = results.sort_values('abs_nes', ascending=False)
    print(f"  Loaded {len(results)} pathways")

    # Top 9, sorted for barh (ascending so positive ends up on top)
    top9 = results.head(9).copy()
    top9 = top9.sort_values('NES', ascending=True)
    top9['clean_name'] = (top9['Term'].str.replace('HALLMARK_', '')
                          .str.replace('_', ' ').replace(PRINTED_NAMES))

    colors = [COLOR_POS if nes > 0 else COLOR_NEG for nes in top9['NES']]

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    y_pos = np.arange(len(top9))
    ax.barh(y_pos, top9['NES'], color=colors, edgecolor='black',
            linewidth=style.EDGE_PT)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top9['clean_name'])
    ax.set_xlabel('NES (Post-NR/Post-R)')   # the page's label (2026-09-14)

    ax.axvline(x=0, color='black', linewidth=style.RULE_PT)

    max_abs = max(abs(top9['NES'].min()), abs(top9['NES'].max()))
    ax.set_xlim(-max_abs - 0.3, max_abs + 0.3)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, BASE_DIR / 'momac_pathway_enrichment')
    print(f"  Saved: momac_pathway_enrichment.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
