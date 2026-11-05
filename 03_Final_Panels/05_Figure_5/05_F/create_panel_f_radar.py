#!/usr/bin/env python3
"""Figure 5H: Hallmark TNF-alpha/NF-kB NES (method "ttest") in thirteen cell
types, pre- and post-treatment, as a radar plot.

Input:  SOUND13_NFKB_CSV (00_Config/paths.py; a deposited table)
Output: nfkb_radar_celltype_enrichment.{svg,pdf,png} and nfkb_rankings_used.csv
        in this directory
Run:    python create_panel_f_radar.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *  # noqa: E402,F403
from cnsfig import Panel, style  # noqa: E402

SRC = SOUND13_NFKB_CSV
METHOD = "ttest"

SCALE = 4                       # tick-pad scaling constants
SMALL_PT = 4.0

PANEL = Panel(5, "H", __file__, small_pt=SMALL_PT)

#: Ring radius and centre on the canvas, in mm. The outer ring is NES 2.4.
RING_R_MM = 9.8
RING_CX_MM = 23.0
RING_CY_MM = 17.3
#: Radial distance, in NES units, from the outer ring to the start of a name.
SPOKE_PAD = 0.4
COLOR_PRE = '#2166AC'
COLOR_POST = '#B2182B'

LABELS = {"B_cells": "B cells", "DC_cells": "DC",
          "Endothelial_cells": "Endothelial", "Epithelial": "Epithelial",
          "Fibroblast": "Fibroblast", "Mast_cells": "Mast", "MoMac": "MoMac",
          "Neutrophils": "Neutrophils", "NK_cells": "NK",
          "Pericyte": "Pericyte", "Plasma_cells": "Plasma",
          "TCD4_cells": "CD4+ T", "TCD8_cells": "CD8+ T"}
ORDER = ["B_cells", "DC_cells", "Endothelial_cells", "Epithelial", "Fibroblast",
         "Mast_cells", "MoMac", "Neutrophils", "NK_cells", "Pericyte",
         "Plasma_cells", "TCD4_cells", "TCD8_cells"]


def load():
    if not SRC.exists():
        raise SystemExit(f"missing {SRC} - run 07_R1.8_NFkB_Specificity first")
    df = pd.read_csv(SRC)
    df = df[df["method"] == METHOD]
    wide = df.pivot(index="cell_type", columns="phase",
                    values=["nes", "fdr_q", "rank", "n_sets"])
    missing = [c for c in ORDER if c not in wide.index]
    if missing:
        raise SystemExit(f"cell types absent from {SRC.name}: {missing}")
    out = pd.DataFrame({
        "cell_type": ORDER,
        "label": [LABELS[c] for c in ORDER],
        "pre_nes": [wide.loc[c, ("nes", "pre")] for c in ORDER],
        "post_nes": [wide.loc[c, ("nes", "post")] for c in ORDER],
        "pre_fdr": [wide.loc[c, ("fdr_q", "pre")] for c in ORDER],
        "post_fdr": [wide.loc[c, ("fdr_q", "post")] for c in ORDER],
        "pre_rank": [wide.loc[c, ("rank", "pre")] for c in ORDER],
        "post_rank": [wide.loc[c, ("rank", "post")] for c in ORDER],
        "pre_n_sets": [wide.loc[c, ("n_sets", "pre")] for c in ORDER],
        "post_n_sets": [wide.loc[c, ("n_sets", "post")] for c in ORDER],
    })
    if out[["pre_nes", "post_nes"]].isna().any().any():
        raise SystemExit("NES missing for at least one cell type or timepoint")
    return out


def main():
    PANEL.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)

    df = load()

    labels = df['label'].tolist()
    n_vars = len(labels)
    angles = np.linspace(0, 2 * np.pi, n_vars, endpoint=False).tolist()
    angles_closed = angles + [angles[0]]

    pre_nes = df['pre_nes'].tolist()
    post_nes = df['post_nes'].tolist()
    pre_nes_closed = pre_nes + [pre_nes[0]]
    post_nes_closed = post_nes + [post_nes[0]]

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm,
                                subplot_kw=dict(polar=True))

    ax.fill(angles_closed, pre_nes_closed, color=COLOR_PRE, alpha=0.10)
    ax.fill(angles_closed, post_nes_closed, color=COLOR_POST, alpha=0.10)
    ax.plot(angles_closed, pre_nes_closed, 'o--', linewidth=style.RULE_PT,
            color=COLOR_PRE, label='Pre-treatment',
            markersize=0.41 * style.PT_PER_MM)

    ax.plot(angles_closed, post_nes_closed, 'o-', linewidth=style.RULE_PT,
            color=COLOR_POST, label='Post-treatment',
            markersize=0.41 * style.PT_PER_MM)

    # Spoke names drawn as text just outside the ring (a polar axis re-lays
    # its tick labels at draw time).
    ax.set_xticks(angles)
    ax.set_xticklabels([])
    r_label = 2.4 + SPOKE_PAD
    for name, ang in zip(labels, angles):
        c, s_ = np.cos(ang), np.sin(ang)
        ha = 'center' if abs(c) < 0.2 else ('left' if c > 0 else 'right')
        va = 'center' if abs(s_) < 0.2 else ('bottom' if s_ > 0 else 'top')
        ax.text(ang, r_label, name, ha=ha, va=va, fontsize=style.tick_pt())
    ax.tick_params(axis='x', pad=2 * SCALE * MARK)
    ax.set_ylim(-2, 2.4)
    # Label rings 0 and 2, in the gap between spokes that crosses the fewest
    # curves, markers or names.
    ax.set_yticks([0, 2])
    ax.set_yticklabels(['0', '2'], color='gray')
    mids = [(a + b) / 2 for a, b in zip(angles, angles[1:] + [angles[0] + 2 * np.pi])]

    def crossings(deg):
        ax.set_rlabel_position(deg)
        fig.canvas.draw()
        rr = fig.canvas.get_renderer()
        # Padded: PDF glyph boxes are taller than the rendered ink.
        boxes_ = [t.get_window_extent(rr).expanded(1.6, 2.2)
                  for t in ax.yaxis.get_ticklabels() if t.get_text()]
        hits = 0
        for bb in boxes_:
            for ln in ax.lines:
                path = ln.get_path().transformed(ln.get_transform())
                if path.intersects_bbox(bb, filled=False):
                    hits += 1
            for t in ax.texts:
                if t.get_window_extent(rr).overlaps(bb):
                    hits += 1
        return hits
    scores = [(crossings(np.degrees(m) % 360), k) for k, m in enumerate(mids)]
    best = min(scores)[1]
    ax.set_rlabel_position(np.degrees(mids[best]) % 360)
    ax.legend(loc='lower right', bbox_to_anchor=(1.0 - 0.4 / PANEL.w_mm, 0.4 / PANEL.h_mm),
              bbox_transform=fig.transFigure,
              handlelength=1.4, borderpad=0.2, labelspacing=0.3,
              edgecolor='0.6', framealpha=1.0, frameon=True)
    ax.grid(True, linestyle='-', alpha=0.3)
    ax.set_title('')
    fig.text(RING_CX_MM / PANEL.w_mm, 1.0 - 0.4 / PANEL.h_mm, 'NES',
             ha='center', va='top', fontsize=style.body_pt())

    theta_circle = np.linspace(0, 2 * np.pi, 100)
    ax.plot(theta_circle, [0] * 100, 'k--', linewidth=style.RULE_PT, alpha=0.7)

    ax.set_position([(RING_CX_MM - RING_R_MM) / PANEL.w_mm,
                     1.0 - (RING_CY_MM + RING_R_MM) / PANEL.h_mm,
                     2 * RING_R_MM / PANEL.w_mm, 2 * RING_R_MM / PANEL.h_mm])
    PANEL.finish(fig, 'nfkb_radar_celltype_enrichment')

    df.to_csv(PANEL.dir / 'nfkb_rankings_used.csv', index=False)


if __name__ == "__main__":
    main()
