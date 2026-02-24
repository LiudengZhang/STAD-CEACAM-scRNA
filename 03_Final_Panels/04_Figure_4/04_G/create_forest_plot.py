#!/usr/bin/env python3
"""
Panel G: Forest Plot with Diverging Colors
Fixed version: Blue for FC<1 (NR-enriched), Red for FC>1 (R-enriched)
Includes x-axis label "Fold Change (R/NR)"
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *

# 4× scaling method
SCALE = 4
DPI = 300

# Target print size: ~4.5 x 6.4 cm (wider, matches Row 3 height = 64mm)
FIGURE_SIZE = (4.5 * SCALE / 2.54, 6.4 * SCALE / 2.54)

# Font sizes
LABEL_FONTSIZE = 5 * SCALE
TICK_FONTSIZE = 5 * SCALE
TITLE_FONTSIZE = 6 * SCALE

# Colors (matching Figure 2 Panel A diverging scheme)
COLOR_PROTECTIVE = '#0072B2'  # Blue for OR < 1
COLOR_RISK = '#D55E00'        # Red/Orange for OR > 1
COLOR_TOTAL = '#009E73'       # Green for Total MoMac
COLOR_NEUTRAL = '#666666'     # Gray for non-significant

def get_short_label(cell_state):
    """Get shortened label for display."""
    if cell_state == 'Total_MoMac':
        return 'Total MoMac'

    parts = cell_state.split('_')
    if len(parts) >= 3:
        short = f"{parts[0]}_{parts[1]}"
    else:
        short = cell_state

    # Strip the Cx_ prefix from the short label
    if '_' in short:
        prefix = short.split('_')[0]
        if len(prefix) >= 2 and prefix[0] == 'C' and prefix[1:].isdigit():
            short = short[len(prefix) + 1:]

    label_mappings = {
        'MoMac': 'Mac',
        'Plasma': 'PB'
    }

    return label_mappings.get(short, short)

def main():
    print("=" * 60)
    print("Panel F: Forest Plot (Diverging Colors)")
    print("=" * 60)

    script_dir = Path(__file__).parent

    # Load data from central config paths
    results_df = pd.read_csv(FIG4_MOMAC_FOLDCHANGE)
    results_df = results_df.sort_values('Fold_Change', ascending=False)

    print(f"\nLoaded {len(results_df)} cell states")

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': TICK_FONTSIZE,
        'axes.linewidth': 1.5 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    y_positions = np.arange(len(results_df))

    for idx, row in enumerate(results_df.itertuples()):
        is_total = row.Cell_State == 'Total_MoMac'
        fc = row.Fold_Change
        ci_l = row.CI_Lower
        ci_u = row.CI_Upper
        p_val = row.P_Value

        if np.isnan(fc) or np.isnan(ci_l) or np.isnan(ci_u):
            continue

        # Determine color based on fold change direction and Total MoMac
        if is_total:
            color = COLOR_TOTAL  # Green for Total MoMac
        else:
            color = COLOR_PROTECTIVE if fc < 1 else COLOR_RISK

        # Marker and line properties (same for all)
        line_width = 1.0 * SCALE
        marker_size = 35 * SCALE

        # Plot whisker (CI)
        ax.plot([ci_l, ci_u], [idx, idx], color=color, linewidth=line_width,
                alpha=0.85, zorder=2, solid_capstyle='butt')

        # Add caps
        cap_height = 0.15
        ax.plot([ci_l, ci_l], [idx - cap_height, idx + cap_height], color=color,
                linewidth=line_width, alpha=0.85, zorder=2)
        ax.plot([ci_u, ci_u], [idx - cap_height, idx + cap_height], color=color,
                linewidth=line_width, alpha=0.85, zorder=2)

        # Plot square marker (same style for all, edge = fill color)
        ax.scatter(fc, idx, s=marker_size, marker='s', color=color, alpha=0.85,
                  edgecolors=color, linewidths=1.2 * SCALE, zorder=4)

    # Reference line at FC=1.0
    ax.axvline(x=1, color='black', linestyle='--', linewidth=1.0 * SCALE, alpha=0.7, zorder=1)

    # Y-axis labels
    short_labels = [get_short_label(state) for state in results_df['Cell_State']]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(short_labels, fontsize=TICK_FONTSIZE)

    # No bold formatting - keep all labels consistent
    # for i, label in enumerate(ax.get_yticklabels()):
    #     if "Total" in label.get_text():
    #         label.set_fontweight('bold')

    # X-axis label
    ax.set_xlabel('Fold Change (R/NR)', fontsize=LABEL_FONTSIZE, labelpad=5 * SCALE)

    # X-axis settings
    ax.set_xscale('log')
    ax.set_xlim(0.04, 12)
    ax.tick_params(axis='x', labelsize=TICK_FONTSIZE, width=1.5 * SCALE, length=4 * SCALE)
    ax.tick_params(axis='y', width=1.5 * SCALE, length=4 * SCALE)

    # Grid - only horizontal lines, vertical reference line at x=1 is kept via axvline
    ax.grid(axis='y', alpha=0.1, linestyle=':', zorder=0, color='gray')

    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5 * SCALE)
    ax.spines['bottom'].set_linewidth(1.5 * SCALE)

    plt.tight_layout()

    # Save
    output_path = script_dir / 'momac_forest_plot_final.png'
    fig.savefig(output_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(output_path.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_path}")

    plt.close()

if __name__ == "__main__":
    main()
