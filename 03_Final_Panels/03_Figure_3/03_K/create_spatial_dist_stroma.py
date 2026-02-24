#!/usr/bin/env python3
"""
Panel H: Spatial Distance to Stroma (single panel, sample_03)
Extracted from old 03_C composite. 4x scaling — ONLY canvas + fonts scaled.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA

DPI = 300
SCALE = 4
PANEL_WIDTH_CM = 5.0 * SCALE
PANEL_HEIGHT_CM = 5.0 * SCALE
CM_TO_INCH = 1 / 2.54

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_NAME = 'sample_03'


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("=" * 60)
    print("Panel H: Spatial Distance to Stroma")
    print("=" * 60)

    print("\n[1/2] Loading spot data...")
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    sample_df = df[df['sample'] == SAMPLE_NAME].copy()
    print(f"  {SAMPLE_NAME}: {len(sample_df)} spots")

    print("\n[2/2] Creating visualization...")
    fig_w = PANEL_WIDTH_CM * CM_TO_INCH
    fig_h = PANEL_HEIGHT_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sc = ax.scatter(sample_df['x'], sample_df['y'],
                    c=sample_df['distance_to_stroma'],
                    cmap='viridis', s=6, alpha=0.8)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.solids.set_rasterized(False)
    cbar.set_label('Dist to Stroma', fontsize=7 * SCALE)
    cbar.ax.tick_params(labelsize=6 * SCALE, width=1.0, length=4)
    cbar.outline.set_linewidth(1.0)

    ax.set_title('Distance to Stroma', fontsize=8 * SCALE, fontweight='normal')
    ax.set_xlabel('X coordinate', fontsize=7 * SCALE)
    ax.set_ylabel('Y coordinate', fontsize=7 * SCALE)
    ax.tick_params(axis='both', labelsize=6 * SCALE, width=1.0, length=4)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    ax.set_aspect('equal')

    plt.tight_layout()

    out = os.path.join(BASE_DIR, 'spatial_dist_stroma.png')
    plt.savefig(out, dpi=DPI, facecolor='white', bbox_inches='tight')
    plt.savefig(out.replace('.png', '.svg'), format='svg', facecolor='white', bbox_inches='tight')
    plt.savefig(out.replace('.png', '.pdf'), dpi=DPI, facecolor='white', bbox_inches='tight')
    print(f"\n  Saved: {out}")
    plt.close()

    print("\nDone!")


if __name__ == '__main__':
    main()
