#!/usr/bin/env python3
"""
S4_F: Stacked cell type proportion bar chart per sample from spatial deconvolution.
Shows GraphST deconvolution results averaged per sample for all 10 GSE251950 samples.
Color mapping matches Figure 1C (scRNA-seq stacked bar).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import SPATIAL_SPOT_DATA

SCALE = 4
DPI = 300
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 16.0 * SCALE
PANEL_HEIGHT_CM = 6.0 * SCALE

OUTPUT_DIR = Path(__file__).parent

SAMPLE_INFO = [
    ('sample_01', 'GC1'), ('sample_02', 'GC2'), ('sample_03', 'GC3'),
    ('sample_04', 'GC4'), ('sample_05', 'GC5'), ('sample_06', 'GC6'),
    ('sample_07', 'GC6-PM'), ('sample_08', 'GC7'), ('sample_09', 'GC8'),
    ('sample_10', 'GC9'),
]

# Cell type colors matching Figure 1C + spatial-specific types
CELL_TYPE_COLORS = {
    'Epi CEACAM-high': '#d95f02',
    'Epi CEACAM-low': '#e5c494',
    'Monocytes/Macrophages': '#fc8d62',
    'CD4+ T cells': '#66c2a5',
    'CD8+ T cells': '#1b9e77',
    'NK cells': '#7570b3',
    'B cells': '#ffd92f',
    'Plasma cells': '#8da0cb',
    'DC cells': '#80b1d3',
    'Endothelial cells': '#a6d854',
    'Fibroblast': '#e78ac3',
    'Pericyte': '#bebada',
    'Neutrophils': '#b3b3b3',
    'Mast cells': '#fb8072',
}

# Order for stacking (epithelial first, then immune, then stroma)
CELL_TYPE_ORDER = [
    'Epi CEACAM-high', 'Epi CEACAM-low',
    'CD4+ T cells', 'CD8+ T cells', 'NK cells',
    'B cells', 'Plasma cells',
    'Monocytes/Macrophages', 'DC cells',
    'Neutrophils', 'Mast cells',
    'Endothelial cells', 'Fibroblast', 'Pericyte',
]


def main():
    print("S4_F: Stacked Cell Type Proportion Bar Chart")

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    df = pd.read_csv(SPATIAL_SPOT_DATA)

    # Get cell type columns
    ct_cols = [c for c in CELL_TYPE_ORDER if c in df.columns]
    print(f"  Cell types: {len(ct_cols)}")

    # Mean proportion per sample
    mean_by_sample = df.groupby('sample')[ct_cols].mean()

    # Reorder samples
    sample_order = [s[0] for s in SAMPLE_INFO]
    labels = [s[1] for s in SAMPLE_INFO]
    mean_by_sample = mean_by_sample.loc[sample_order]

    fig, ax = plt.subplots(figsize=(
        PANEL_WIDTH_CM * CM_TO_INCH,
        PANEL_HEIGHT_CM * CM_TO_INCH
    ))

    x = np.arange(len(sample_order))
    bottom = np.zeros(len(sample_order))

    for ct in ct_cols:
        vals = mean_by_sample[ct].values
        color = CELL_TYPE_COLORS.get(ct, '#999999')
        ax.bar(x, vals, bottom=bottom, color=color, label=ct,
               width=0.75, edgecolor='white', linewidth=0.3)
        bottom += vals

    ax.set_xticks(x)
    label_colors = ['#B2182B' if 'PM' in lbl else 'black' for lbl in labels]
    ax.set_xticklabels(labels, fontsize=5 * SCALE, rotation=45, ha='right')
    for tick, color in zip(ax.get_xticklabels(), label_colors):
        tick.set_color(color)
        if color == '#B2182B':
            tick.set_fontweight('bold')

    ax.set_ylabel('Mean proportion', fontsize=6 * SCALE)
    ax.set_ylim(0, 1)
    ax.set_title('GraphST Deconvolution: Cell Type Proportions per Sample',
                 fontsize=6.5 * SCALE)

    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
              fontsize=4.5 * SCALE, frameon=False, ncol=1)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5 * SCALE)
    ax.spines['bottom'].set_linewidth(0.5 * SCALE)
    ax.tick_params(axis='both', labelsize=5 * SCALE,
                   width=0.5 * SCALE, length=3 * SCALE)

    plt.tight_layout()

    out = OUTPUT_DIR / 'panel_S4_F'
    fig.savefig(str(out) + '.png', dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(str(out) + '.svg', format='svg', bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}.png / .svg")

    # Print summary
    print("\nMean proportions per sample:")
    print(mean_by_sample.round(3).to_string())
    plt.close()


if __name__ == '__main__':
    main()
