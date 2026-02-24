#!/usr/bin/env python3
"""
Panel F: NF-κB NES Radar Plot — Pre vs Post treatment (NR vs R).
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
PANEL_SIZE_CM = 5.0 * SCALE  # square

INPUT_CSV = NFKB_RANKINGS_CSV
MAX_PATHWAYS = 50

# Pre vs Post colors
COLOR_PRE = '#2166AC'
COLOR_POST = '#B2182B'


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    df = pd.read_csv(INPUT_CSV)

    # NES values are already in NR vs R direction (from t-test DEGs)
    labels = df['label'].tolist()
    n_vars = len(labels)
    angles = np.linspace(0, 2 * np.pi, n_vars, endpoint=False).tolist()
    angles_closed = angles + [angles[0]]

    pre_nes = df['pre_nes'].fillna(0).tolist()
    post_nes = df['post_nes'].fillna(0).tolist()
    pre_nes_closed = pre_nes + [pre_nes[0]]
    post_nes_closed = post_nes + [post_nes[0]]

    fig_size = PANEL_SIZE_CM * CM_TO_INCH
    fig, ax = plt.subplots(figsize=(fig_size, fig_size), subplot_kw=dict(polar=True))

    ax.plot(angles_closed, pre_nes_closed, 'o--', linewidth=2,
            color=COLOR_PRE, label='Pre-treatment', markersize=6)

    ax.plot(angles_closed, post_nes_closed, 'o-', linewidth=2,
            color=COLOR_POST, label='Post-treatment', markersize=6)

    ax.set_xticks(angles)
    ax.set_xticklabels(labels, size=6 * SCALE)
    ax.set_ylim(-2, 2)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_yticklabels(['-2', '-1', '0', '1', '2'], size=5 * SCALE, color='gray')
    ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1.1), fontsize=5 * SCALE)
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

    # Save a copy of the rankings used
    df.to_csv(BASE_DIR / 'nfkb_rankings_used.csv', index=False)


if __name__ == "__main__":
    main()
