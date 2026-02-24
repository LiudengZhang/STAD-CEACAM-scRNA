#!/usr/bin/env python3
"""
Panel Q: GSEA Summary — 4 separate dotplots (MoMac, Epithelial, Fibroblast, DC).
Style adapted from Figure 2 checkpoint receptor dotplot:
  - Pathways on X-axis (rotated 45°)
  - NES on Y-axis
  - Dot size: significance tier (p≤0.05, p≤0.10, p>0.10)
  - Dot color: NES value (PuOr_r colormap)
  - Reference line at NES=0
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / 'gsea_data'
PANEL_DIR = BASE_DIR
SCALE = 4
CM = 1 / 2.54
DPI = 300

# 8 pathways in display order
PATHWAYS = [
    'TNF-alpha Signaling via NF-kB',
    'Inflammatory Response',
    'IL-6/JAK/STAT3 Signaling',
    'Interferon Gamma Response',
    'Epithelial Mesenchymal Transition',
    'Angiogenesis',
    'Hypoxia',
    'Apoptosis',
]

PATHWAY_LABELS = [
    'TNF-\u03b1/NF-\u03baB',
    'Inflammatory\nResponse',
    'IL-6/JAK/\nSTAT3',
    'IFN-\u03b3\nResponse',
    'EMT',
    'Angio-\ngenesis',
    'Hypoxia',
    'Apoptosis',
]

# 4 cell types in order
CELL_TYPES = ['MoMac', 'Epithelial', 'Fibroblast', 'DC']
CELL_LABELS = {
    'MoMac': 'Monocytes/\nMacrophages',
    'Epithelial': 'Epithelial',
    'Fibroblast': 'Fibroblast',
    'DC': 'DC',
}


def assign_dot_size(p_value):
    if p_value <= 0.05:
        return 150 * SCALE
    elif p_value <= 0.10:
        return 100 * SCALE
    else:
        return 60 * SCALE


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * SCALE,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("=" * 60)
    print("Panel Q: 4× GSEA Dotplots (Fig 2 style)")
    print("=" * 60)

    # Load and combine data
    existing = pd.read_csv(DATA_DIR / 'gsea_combined_5types.csv')
    momac = pd.read_csv(DATA_DIR / 'gsea_momac.csv')
    df = pd.concat([existing, momac], ignore_index=True)
    df['NES'] = df['NES'].astype(float)
    df['NOM_pval'] = df['NOM_pval'].astype(float)

    # Compute global NES range for consistent color mapping
    subset = df[df['CellType'].isin(CELL_TYPES) & df['Pathway'].isin(PATHWAYS)]
    nes_abs_max = max(abs(subset['NES'].min()), abs(subset['NES'].max()))
    vmin, vmax = -nes_abs_max, nes_abs_max

    cmap = plt.cm.PuOr_r  # warm orange = positive NES (enriched in NR)

    # Short labels keyed by full pathway name
    LABEL_MAP = dict(zip(PATHWAYS, PATHWAY_LABELS))

    for ct in CELL_TYPES:
        ct_df = df[df['CellType'] == ct]
        print(f"\n  {ct}:")

        # Collect data for the 8 pathways
        rows = []
        for pw in PATHWAYS:
            row = ct_df[ct_df['Pathway'] == pw]
            if len(row) > 0:
                r = row.iloc[0]
                rows.append({'pathway': pw, 'nes': float(r['NES']), 'pval': float(r['NOM_pval'])})
            else:
                rows.append({'pathway': pw, 'nes': 0.0, 'pval': 1.0})

        # Sort by NES descending (each cell type gets its own order)
        rows.sort(key=lambda r: r['nes'], reverse=True)

        sorted_labels = [LABEL_MAP[r['pathway']] for r in rows]
        nes_vals = np.array([r['nes'] for r in rows])
        p_vals = [r['pval'] for r in rows]

        fig_w = 6.5 * SCALE * CM
        fig_h = 5.5 * SCALE * CM
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        x_positions = np.arange(len(rows))
        sizes = [assign_dot_size(p) for p in p_vals]

        scatter = ax.scatter(
            x_positions, nes_vals, s=sizes, c=nes_vals, cmap=cmap, alpha=0.85,
            edgecolors='black', linewidths=0.5, vmin=vmin, vmax=vmax, zorder=3
        )

        # Reference line at NES=0
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
        ax.grid(True, axis='y', alpha=0.3, linestyle=':', linewidth=0.3)

        # X-axis — sorted labels
        ax.set_xticks(x_positions)
        ax.set_xticklabels(sorted_labels, fontsize=4.5 * SCALE, rotation=45, ha='right')
        ax.set_xlim(-0.8, len(rows) - 0.2)

        # Y-axis
        ax.set_ylim(0.5, 2.2)
        ax.set_ylabel('NES', fontsize=6 * SCALE)

        # Title
        ax.set_title(CELL_LABELS[ct], fontsize=6 * SCALE, fontweight='normal', pad=8)

        # Spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.tick_params(axis='both', labelsize=5 * SCALE, width=0.5, length=4)

        # Print values (in sorted order)
        for i, r in enumerate(rows):
            sig = '*' if r['pval'] < 0.05 else ''
            short = r['pathway'][:25]
            print(f"    {short:28s} NES={r['nes']:+.3f}  p={r['pval']:.4f} {sig}")

        plt.tight_layout()

        stem = f'gsea_dotplot_{ct.lower().replace("/", "_")}'
        for ext in ['svg', 'png']:
            kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
            if ext == 'png':
                kw['dpi'] = DPI
            fig.savefig(PANEL_DIR / f'{stem}.{ext}', **kw)
        print(f"    Saved: {stem}")
        plt.close()

    # Also generate shared legend panel
    _make_legend(cmap, vmin, vmax)


def _make_legend(cmap, vmin, vmax):
    """Create a small standalone legend panel with size + colorbar."""
    fig, ax = plt.subplots(figsize=(3.0 * SCALE * CM, 4.5 * SCALE * CM))
    ax.set_axis_off()

    # Size legend
    y_start = 0.92
    ax.text(0.05, y_start + 0.05, 'Significance', fontsize=5.5 * SCALE,
            fontweight='bold', transform=ax.transAxes, va='top')
    for i, (label, size) in enumerate([
        ('p \u2264 0.05', 150 * SCALE),
        ('p \u2264 0.10', 100 * SCALE),
        ('p > 0.10', 60 * SCALE),
    ]):
        y = y_start - 0.12 * (i + 1)
        ax.scatter(0.15, y, s=size, c='gray', alpha=0.6, edgecolors='black',
                   linewidths=0.5, transform=ax.transAxes, zorder=3)
        ax.text(0.35, y, label, fontsize=4.5 * SCALE, va='center',
                transform=ax.transAxes)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cax = fig.add_axes([0.15, 0.08, 0.6, 0.04])
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label('NES', fontsize=5 * SCALE)
    cbar.ax.tick_params(labelsize=4.5 * SCALE)

    ax.text(0.05, 0.22, 'Enrichment', fontsize=5.5 * SCALE,
            fontweight='bold', transform=ax.transAxes, va='top')
    ax.text(0.05, 0.15, 'Positive = NR enriched', fontsize=4 * SCALE,
            fontstyle='italic', color='gray', transform=ax.transAxes, va='top')

    stem = 'gsea_legend'
    for ext in ['svg', 'png']:
        kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if ext == 'png':
            kw['dpi'] = DPI
        fig.savefig(PANEL_DIR / f'{stem}.{ext}', **kw)
    print(f"\n    Saved: {stem}")
    plt.close()


if __name__ == '__main__':
    main()
