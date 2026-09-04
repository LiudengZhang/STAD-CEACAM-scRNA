#!/usr/bin/env python3
"""
Figure 5 panel A, RESTYLED (Version B) - MoMac GSEA horizontal barplot, top 9
Hallmark gene sets by |NES|.

Version A is
`03_Final_Panels/05_Figure_5/05_A/create_panel_a_momac_enrichment.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

THE INPUT PATH IS NOT TOUCHED, AND MUST NOT BE.
-----------------------------------------------
Version A's own header records why, and it is repeated here because this is the
panel this project retracted twice. It reads the MAST prerank GSEA table
`GSEA/MoMac_mast_prerank_gsea.csv`. The twelve per-cell-type
`*_mast_prerank_gsea.csv` were moved under `GSEA/_archived/` after the figures
were made and this script was once repointed at `GSEA/post/MoMac_gsea_hallmark.csv`
- a different run, built on the doubly normalised .X. That repoint, not the
figure, is why panel A stopped reproducing, and twice led an agent to conclude
the published panel was wrong. Both conclusions were retracted. Against the
restored table the nine bars reproduce the printed panel to within 0.0001 NES.

Do not repoint this at `GSEA/post/` or at
`05_GSEA_Summary/gsea_data/gsea_momac.csv`. A restyle changes how a panel is
drawn and nothing else; if the redraw ever appeared to disagree with the
published panel, the bug would be in this copy, not in the figure.

Every value read, sorted, selected and plotted is Version A's. The drawing code
is the same code.

  printed panel  Figure 5 A        (PROVENANCE.csv; "05_A" happens to agree,
                 but it was looked up, not inferred)
  printed rect   58.7 x 40.6 mm    (panel_rects.csv)
  Version B box  100.0 x 50.0 mm

    The nine Hallmark names are the panel's widest element. The longest,
    "TNF-alpha Signaling via NF-kB", sets about 35 mm at 7 pt, so the label
    column alone is more than half of the published 58.7 mm width. Nine
    categories also need about 4 mm of height each to be read at 7 pt.

MARK
    Version A drew 24 x 16.8 cm at SCALE = 4 and set its smallest body type -
    the y tick labels, the x tick labels and the x axis label, all three - at
    `6 * SCALE`, so SMALL_PT = 6 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 24 = 0.2917

    The bar edge and the zero rule are scaled by it. Tick widths/lengths and
    spine widths are not: those are style, and cnsplots sets them
    (axes.linewidth 0.5, ticks size 2 width 0.6).

Note on the x axis label: Version A deliberately writes 'NES' and not the
contrast string the printed panel carries, because no script on disk emits that
string - it was added at assembly. That comment and that behaviour are kept.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *
import panel_style_cns as style

BASE_DIR = Path(__file__).parent

SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 6.0                  # Version A's smallest body type, before * SCALE

PRINTED_MM = (58.7, 40.6)       # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 100.0, 50.0
MARGIN = dict(left=38.0, right=2.0, top=1.5, bottom=9.0)

# Standard R/NR colors
COLOR_POS = '#B2182B'   # NR-upregulated (positive NES) — red
COLOR_NEG = '#2166AC'   # R-upregulated (negative NES) — blue

GSEA_CSV = GSEA_DIR / "MoMac_mast_prerank_gsea.csv"


def main():
    family = style.apply()
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
    top9['clean_name'] = top9['Term'].str.replace('HALLMARK_', '').str.replace('_', ' ')

    colors = [COLOR_POS if nes > 0 else COLOR_NEG for nes in top9['NES']]

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    y_pos = np.arange(len(top9))
    ax.barh(y_pos, top9['NES'], color=colors, edgecolor='black',
            linewidth=0.5 * MARK)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top9['clean_name'])
    # The printed panel names the contrast here - 'NES (Post-NR/Post-R)'. No script
    # on disk emits that string; it was added when the figure was assembled, and
    # writing it here instead would be inventing provenance for it.
    ax.set_xlabel('NES')

    ax.axvline(x=0, color='black', linewidth=0.8 * MARK)

    max_abs = max(abs(top9['NES'].min()), abs(top9['NES'].max()))
    ax.set_xlim(-max_abs - 0.3, max_abs + 0.3)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, BASE_DIR / 'momac_pathway_enrichment')
    print(f"  Saved: momac_pathway_enrichment.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
