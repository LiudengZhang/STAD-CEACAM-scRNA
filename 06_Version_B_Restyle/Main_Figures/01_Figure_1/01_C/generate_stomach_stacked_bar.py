#!/usr/bin/env python3
# Paths below refer to the upstream Round_4 processing pipeline, which is
# not part of this release. This script is included as a record of how the
# input was produced; it is not called by _run_all_panels.sh.
"""
Figure 1 panel C, RESTYLED (Version B) - cell type composition per stomach
sample.

Version A is
`03_Final_Panels/01_Figure_1/01_C/generate_stomach_stacked_bar.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every ordering and every string is Version A's.
The drawing code is the same code.

  printed panel  Figure 1 C     (PROVENANCE.csv; NOT inferred from "01_C")
  printed rect   160.9 x 56.8 mm   (panel_rects.csv)
  Version B box  118.0 x 56.0 mm

MARK
    Version A drew a 10 x 4 inch canvas (254.0 x 101.6 mm) with no SCALE
    constant, and set three different type sizes on it: 12 pt tick labels,
    7 pt axis label, 5 pt legend. The assembler then fitted the saved SVG
    (251.8 mm wide) into 160.9 mm, a fit of 0.639 - so those three printed at
    7.67, 4.47 and 3.20 pt. The panel's smallest body type is the 5 pt legend,
    so SMALL_PT = 5 and

        MARK = tick_pt / SMALL_PT = 7 / 5 = 1.4

    The one non-type length here is the 0.3 pt white rule between the stacked
    segments; it printed at 0.19 pt and now prints at 0.42, the same 2.19x the
    type grew by.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import FULL_DATASET_H5AD                      # noqa: E402
import panel_style_cns as style                          # noqa: E402

# Sample ID mapping from Supplementary Table 1
ST1_CSV = Path('/path/to/Project_4_05232025/Round_5/04_Manuscript/04_Tables/ST1_patient_sample_characteristics.csv')

# Set2 + Set3 palette (12 major cell types, consistent across Fig 1B & 1C)
CELL_TYPE_COLORS = {
    'Epithelial cells': '#e5c494',
    'T/NK cells': '#66c2a5',
    'Monocytes/Macrophages': '#fc8d62',
    'Plasma cells': '#8da0cb',
    'B cells': '#ffd92f',
    'Endothelial cells': '#a6d854',
    'Fibroblasts': '#e78ac3',
    'Neutrophils': '#b3b3b3',
    'Mast cells': '#fb8072',
    'Pericytes': '#bebada',
    'Dendritic cells': '#80b1d3',
    'Hepatocytes': '#bc80bd',
}

DATA_FILE = FULL_DATASET_H5AD
OUTPUT_DIR = Path(__file__).parent

PRINTED_MM = (160.9, 56.8)          # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 118.0, 56.0
MARGIN = dict(left=11.0, right=32.0, top=1.5, bottom=13.5)

SMALL_PT = 5.0                      # Version A's smallest body type


def main():
    family = style.apply()
    MARK = style.tick_pt() / SMALL_PT
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("=" * 80)
    print("GENERATING CELL TYPE STACKED BAR - STOMACH SAMPLES ONLY")
    print("=" * 80)

    print(f"\nLoading data from: {DATA_FILE}")
    adata = sc.read_h5ad(DATA_FILE)
    print(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    adata_stomach = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    print(f"Stomach samples: {adata_stomach.shape[0]:,} cells")
    print(f"Number of samples: {adata_stomach.obs['sample'].nunique()}")

    adata_stomach.obs['major_cell_type_merged'] = adata_stomach.obs['major_cell_type'].replace({
        'CD4+ T cells': 'T/NK cells',
        'CD8+ T cells': 'T/NK cells',
        'NK cells': 'T/NK cells',
        'DC cells': 'Dendritic cells',
        'Fibroblast': 'Fibroblasts',
        'Pericyte': 'Pericytes',
        'Hepatocyte': 'Hepatocytes'
    })

    proportions = adata_stomach.obs.groupby(['sample', 'major_cell_type_merged'], observed=True).size().unstack(fill_value=0)
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    st1 = pd.read_csv(ST1_CSV)
    st1_stomach = st1[st1['Anatomical site'] == 'Stomach']
    orig_to_sample = dict(zip(st1_stomach['Original Sample ID'], st1_stomach['Sample']))

    sample_order = sorted(proportions.index.tolist())
    proportions = proportions.reindex(sample_order)

    x_labels = [orig_to_sample.get(s, s) for s in sample_order]

    cell_type_order = proportions.mean().sort_values(ascending=False).index.tolist()
    proportions = proportions[cell_type_order]

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    positions = np.arange(len(proportions))
    bottom = np.zeros(len(proportions))
    for cell_type in cell_type_order:
        color = CELL_TYPE_COLORS.get(cell_type, '#999999')
        ax.bar(positions, proportions[cell_type], bottom=bottom,
               color=color, label=cell_type, width=0.8, edgecolor='white',
               linewidth=0.3 * MARK)
        bottom += proportions[cell_type].values

    ax.set_xticks(positions)
    ax.set_xticklabels(x_labels, rotation=90, ha='center')
    ax.set_ylabel('Proportion')
    ax.set_ylim(0, 1)

    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False,
              handlelength=1.2, handletextpad=0.4, labelspacing=0.35)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUTPUT_DIR / '01_C_sample_balance_stomach')
    print(f"\nSaved: {OUTPUT_DIR / '01_C_sample_balance_stomach'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")
    print(f"  Samples: {len(proportions)}")
    print("=" * 80)


if __name__ == '__main__':
    main()
