#!/usr/bin/env python3
"""
S4_E: Distance to immune spatial maps for all 10 GSE251950 samples.
2x5 grid, sample_07 (GC6-PM) labeled as Peritoneal Metastasis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import SPATIAL_SPOT_DATA

SCALE = 4
DPI = 300
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 18.0 * SCALE
PANEL_HEIGHT_CM = 8.0 * SCALE

OUTPUT_DIR = Path(__file__).parent

SAMPLE_INFO = [
    ('sample_01', 'GC1'), ('sample_02', 'GC2'), ('sample_03', 'GC3'),
    ('sample_04', 'GC4'), ('sample_05', 'GC5'), ('sample_06', 'GC6'),
    ('sample_07', 'GC6-PM*'), ('sample_08', 'GC7'), ('sample_09', 'GC8'),
    ('sample_10', 'GC9'),
]


def main():
    print("S4_E: Distance to Immune Spatial Maps (All 10 Samples)")

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    df = pd.read_csv(SPATIAL_SPOT_DATA)
    print(f"  Total spots: {len(df)}")

    fig, axes = plt.subplots(2, 5, figsize=(
        PANEL_WIDTH_CM * CM_TO_INCH,
        PANEL_HEIGHT_CM * CM_TO_INCH
    ))

    vmin = df['distance_to_immune'].quantile(0.02)
    vmax = df['distance_to_immune'].quantile(0.98)

    for idx, (sample_name, label) in enumerate(SAMPLE_INFO):
        row, col = divmod(idx, 5)
        ax = axes[row, col]

        sample_df = df[df['sample'] == sample_name].copy()

        if len(sample_df) == 0:
            ax.text(0.5, 0.5, f'{label}\nNo data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=5 * SCALE)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        sc = ax.scatter(sample_df['x'], sample_df['y'],
                        c=sample_df['distance_to_immune'],
                        cmap='magma_r', s=6, alpha=0.8,
                        vmin=vmin, vmax=vmax, rasterized=True)

        ax.set_title(label, fontsize=5.5 * SCALE, color='black',
                     fontweight='normal', pad=4 * SCALE)

        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.3 * SCALE)

        print(f"  {sample_name} ({label}): {len(sample_df)} spots")

    cbar_ax = fig.add_axes([0.92, 0.25, 0.012, 0.50])
    sm = plt.cm.ScalarMappable(cmap='magma_r',
                                norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Distance to\nimmune', fontsize=5 * SCALE)
    cbar.ax.tick_params(labelsize=4.5 * SCALE, width=0.5 * SCALE)

    fig.text(0.5, 0.01, '*GC6-PM = Peritoneal metastasis (paired with GC6 primary)',
             ha='center', fontsize=4.5 * SCALE, fontstyle='italic', color='#666666')

    plt.subplots_adjust(wspace=0.15, hspace=0.30, right=0.90)

    out = OUTPUT_DIR / 'panel_S4_E'
    fig.savefig(str(out) + '.png', dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(str(out) + '.svg', format='svg', bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}.png / .svg")
    plt.close()


if __name__ == '__main__':
    main()
