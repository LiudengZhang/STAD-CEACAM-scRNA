#!/usr/bin/env python3
"""
Figure 5 panel N, RESTYLED (Version B) - GSEA summary, four dotplots (MoMac,
Epithelial, Fibroblast, DC) plus the shared legend panel.

ONE SCRIPT, FOUR DRAWINGS UNDER ONE PRINTED LETTER. `gsea_dotplot_momac`,
`gsea_dotplot_epithelial`, `gsea_dotplot_fibroblast` and `gsea_dotplot_dc` are
all printed panel N (PROVENANCE.csv). `gsea_legend` is the fifth drawing the
script makes; PROVENANCE.csv does not list it as a panel source, but Version A
produces it and nothing is dropped, so it is restyled with the rest.

Version A is
`03_Revised_Panels/Main_Figures/05_Figure_5/05_GSEA_Summary/create_gsea_4panels.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvases are the
millimetre boxes the panels print in, non-type point sizes are rescaled by MARK
(areas by AREA), margins are millimetres, and the save is `style.save_panel`.

The two input tables, the eight pathways, the four cell types, the per-cell-type
NES sort, the global symmetric colour limits, the PuOr_r colormap and the three
significance tiers are Version A's, unchanged.

`run_gsea_5celltypes.py` and `run_gsea_momac.py` next door are ANALYSIS
scripts: they compute the tables this one reads. They are not panel scripts and
are not restyled.

  printed panel  Figure 5 N        (PROVENANCE.csv)
  printed rect   156.0 x 30.3 mm   (panel_rects.csv, all four dotplots
                 together; ~39 x 30 mm each)
  Version B box  88.0 x 66.0 mm per dotplot; 44.0 x 58.0 mm for the legend

    Eight pathway names set at 45 degrees below the axis are what set the size.
    The longest two are two-line strings ("Inflammatory\\nResponse",
    "Angio-\\ngenesis"); at 7 pt they reach about 18 mm down and to the left of
    their tick, so the panel needs about 20 mm of bottom margin before a dot is
    drawn. 39 mm square cannot hold that.

MARK
    Version A drew every panel at SCALE = 4 and its smallest body type is the
    italic "Positive = NR enriched" note in the legend panel, set at
    `4 * SCALE`, so SMALL_PT = 4 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 16 = 0.4375
        AREA = MARK ** 2                                     = 0.1914

    The three dot areas, the dot edge width, the NES = 0 rule and the y grid
    are scaled by those. Tick widths/lengths and spine widths are not: those
    are style, and cnsplots sets them.

Judgement calls, stated plainly:
  - `fontweight='bold'` on the two legend headings is dropped, per PANEL_SPEC.
    They keep their emphasis by size instead: `body_pt` (8) against `tick_pt`
    (7) for the items, which is the same ordering Version A drew with 5.5 pt
    against 4.5 pt.
  - The dotplot title's explicit `fontweight='normal'` is dropped so cnsplots'
    bold axis title applies.
  - The legend panel's `ax.text` positions and the colorbar's `add_axes`
    rectangle are axes/figure fractions, i.e. layout, and are left alone except
    where the taller box needed them respaced; every string is unchanged.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
import panel_style_cns as style

BASE_DIR = Path(__file__).resolve().parent
# The two GSEA tables are Version A's and are read where they already live.
# They are inputs, not outputs: this script reads them and writes nothing
# back. Copying them under 11_Version_B_Restyle/ would create a second
# copy of a number that the manuscript quotes, which is exactly the kind
# of drift PROVENANCE.csv exists to prevent. Version A's line was
# `DATA_DIR = BASE_DIR / 'gsea_data'`; the directory it named is this one.
DATA_DIR = (Path(__file__).resolve().parents[4] / '03_Revised_Panels'
            / 'Main_Figures' / '05_Figure_5' / '05_GSEA_Summary' / 'gsea_data')
PANEL_DIR = BASE_DIR

SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 4.0                  # Version A's smallest body type, before * SCALE
MARK = style.tick_pt() / (SMALL_PT * SCALE)
AREA = MARK ** 2

PRINTED_MM = (156.0, 30.3)      # published rect, all four dotplots together
PANEL_W_MM, PANEL_H_MM = 88.0, 66.0
MARGIN = dict(left=13.0, right=3.0, top=8.0, bottom=24.0)
LEGEND_W_MM, LEGEND_H_MM = 44.0, 58.0

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
    'TNF-α/NF-κB',
    'Inflammatory\nResponse',
    'IL-6/JAK/\nSTAT3',
    'IFN-γ\nResponse',
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
        return 150 * SCALE * AREA
    elif p_value <= 0.10:
        return 100 * SCALE * AREA
    else:
        return 60 * SCALE * AREA


def main():
    family = style.apply()
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

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

        fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

        x_positions = np.arange(len(rows))
        sizes = [assign_dot_size(p) for p in p_vals]

        scatter = ax.scatter(
            x_positions, nes_vals, s=sizes, c=nes_vals, cmap=cmap, alpha=0.85,
            edgecolors='black', linewidths=0.5 * MARK, vmin=vmin, vmax=vmax, zorder=3
        )

        # Reference line at NES=0
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5 * MARK, alpha=0.7)
        ax.grid(True, axis='y', alpha=0.3, linestyle=':', linewidth=0.3 * MARK)

        # X-axis — sorted labels
        ax.set_xticks(x_positions)
        ax.set_xticklabels(sorted_labels, rotation=45, ha='right')
        ax.set_xlim(-0.8, len(rows) - 0.2)

        # Y-axis
        ax.set_ylim(0.5, 2.2)
        ax.set_ylabel('NES')

        # Title
        ax.set_title(CELL_LABELS[ct], pad=8 * MARK)

        # Spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Print values (in sorted order)
        for i, r in enumerate(rows):
            sig = '*' if r['pval'] < 0.05 else ''
            short = r['pathway'][:25]
            print(f"    {short:28s} NES={r['nes']:+.3f}  p={r['pval']:.4f} {sig}")

        style.margins_mm(fig, **MARGIN)
        over = style.overflow_mm(fig)
        if max(over) > 0.05:
            print(f"    WARNING ink outside the canvas (l,r,b,t mm): "
                  f"{tuple(round(v, 2) for v in over)}")

        stem = f'gsea_dotplot_{ct.lower().replace("/", "_")}'
        style.save_panel(fig, PANEL_DIR / stem)
        print(f"    Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")

    # Also generate shared legend panel
    _make_legend(cmap, vmin, vmax)


def _make_legend(cmap, vmin, vmax):
    """Create a small standalone legend panel with size + colorbar."""
    fig, ax = style.subplots_mm(LEGEND_W_MM, LEGEND_H_MM)
    ax.set_axis_off()

    # Size legend
    y_start = 0.92
    ax.text(0.05, y_start + 0.05, 'Significance', fontsize=style.body_pt(),
            transform=ax.transAxes, va='top')
    for i, (label, size) in enumerate([
        ('p ≤ 0.05', 150 * SCALE * AREA),
        ('p ≤ 0.10', 100 * SCALE * AREA),
        ('p > 0.10', 60 * SCALE * AREA),
    ]):
        y = y_start - 0.12 * (i + 1)
        ax.scatter(0.15, y, s=size, c='gray', alpha=0.6, edgecolors='black',
                   linewidths=0.5 * MARK, transform=ax.transAxes, zorder=3)
        ax.text(0.35, y, label, fontsize=style.tick_pt(), va='center',
                transform=ax.transAxes)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cax = fig.add_axes([0.15, 0.15, 0.6, 0.04])
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label('NES')

    ax.text(0.05, 0.38, 'Enrichment', fontsize=style.body_pt(),
            transform=ax.transAxes, va='top')
    ax.text(0.05, 0.30, 'Positive = NR enriched', fontsize=style.tick_pt(),
            fontstyle='italic', color='gray', transform=ax.transAxes, va='top')

    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"    WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    stem = 'gsea_legend'
    style.save_panel(fig, PANEL_DIR / stem)
    print(f"\n    Saved: {stem}.[svg|pdf|png] at {LEGEND_W_MM} x {LEGEND_H_MM} mm")


if __name__ == '__main__':
    main()
