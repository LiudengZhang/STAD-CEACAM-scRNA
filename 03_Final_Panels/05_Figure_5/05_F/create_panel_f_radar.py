#!/usr/bin/env python3
"""
Panel H of Figure 5: NF-kB NES radar, pre versus post treatment.

REVISED FOR CIR-26-0753-ET. The submitted panel read
GSEA/nfkb_rankings_13types_mast.csv, a prepared table from the differential
expression run that the data audit found to be built on a doubly normalised
.X (00_Data_Audit/FINDINGS.md, sections 1 and 7). Its values disagree with the
clean recompute in two places that the manuscript now makes claims about:
B cells after treatment (+1.01 there, -0.99 in the recompute) and
monocytes/macrophages before treatment (-0.98 there, +1.06 in the recompute).

The panel is therefore rebuilt from
04_Revision_Analyses/12_R1.8_DEG_Recompute, by way of
07_R1.8_NFkB_Specificity/outputs/nfkb_per_celltype.csv, which is the same
table Fig. S9E and Fig. S10C read. The Welch t-test branch is used, matching
the Methods; the MAST branch is the deposited sensitivity analysis.

4x scaling method for Nature Cancer.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *  # noqa: E402,F403
from shared.figure_config import use_panel_style

BASE_DIR = Path(__file__).parent
SRC = (Path(__file__).resolve().parents[4] / "04_Revision_Analyses"
       / "07_R1.8_NFkB_Specificity" / "outputs" / "nfkb_per_celltype.csv")
METHOD = "ttest"

# 4x scaling
DPI = 300
SCALE = 4
CM_TO_INCH = 1 / 2.54
PANEL_SIZE_CM = 5.0 * SCALE  # square

COLOR_PRE = '#2166AC'
COLOR_POST = '#B2182B'

# The names the submitted panel used, so the axis labels do not move.
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
    use_panel_style(font_pt=7)

    df = load()

    labels = df['label'].tolist()
    n_vars = len(labels)
    angles = np.linspace(0, 2 * np.pi, n_vars, endpoint=False).tolist()
    angles_closed = angles + [angles[0]]

    pre_nes = df['pre_nes'].tolist()
    post_nes = df['post_nes'].tolist()
    pre_nes_closed = pre_nes + [pre_nes[0]]
    post_nes_closed = post_nes + [post_nes[0]]

    fig_size = PANEL_SIZE_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_size, fig_size), subplot_kw=dict(polar=True))

    ax.fill(angles_closed, pre_nes_closed, color=COLOR_PRE, alpha=0.10)
    ax.fill(angles_closed, post_nes_closed, color=COLOR_POST, alpha=0.10)
    ax.plot(angles_closed, pre_nes_closed, 'o--', linewidth=2,
            color=COLOR_PRE, label='Pre-treatment', markersize=6)

    ax.plot(angles_closed, post_nes_closed, 'o-', linewidth=2,
            color=COLOR_POST, label='Post-treatment', markersize=6)

    ax.set_xticks(angles)
    ax.set_xticklabels(labels, size=4.1 * SCALE)
    ax.tick_params(axis='x', pad=2 * SCALE)
    ax.set_ylim(-2, 2.4)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_yticklabels(['-2', '-1', '0', '1', '2'], size=4 * SCALE, color='gray')
    ax.set_rlabel_position(90)
    ax.legend(loc='upper left', bbox_to_anchor=(0.99, 0.04), fontsize=4 * SCALE,
              handlelength=1.6, borderpad=0.35, labelspacing=0.35,
              edgecolor='0.6', framealpha=1.0)
    ax.grid(True, linestyle='-', alpha=0.3)
    ax.set_title('NES', fontsize=7 * SCALE, pad=15)

    # Dashed zero circle
    theta_circle = np.linspace(0, 2 * np.pi, 100)
    ax.plot(theta_circle, [0] * 100, 'k--', linewidth=1, alpha=0.7)

    plt.tight_layout()

    output = BASE_DIR / 'nfkb_radar_celltype_enrichment.png'
    plt.savefig(output, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.pdf'), format='pdf', dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output}")

    df.to_csv(BASE_DIR / 'nfkb_rankings_used.csv', index=False)
    print(df[['label', 'pre_nes', 'post_nes']].to_string(index=False))


if __name__ == "__main__":
    main()
