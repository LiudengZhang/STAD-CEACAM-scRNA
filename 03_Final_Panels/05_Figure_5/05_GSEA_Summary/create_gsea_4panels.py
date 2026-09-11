#!/usr/bin/env python3
"""
Figure 5 panel N - GSEA summary: four dotplots (MoMac, epithelial, fibroblast,
DC) and the legend strip they share.

ONE SCRIPT, FIVE DRAWINGS UNDER ONE PRINTED LETTER. `gsea_dotplot_momac`,
`gsea_dotplot_epithelial`, `gsea_dotplot_fibroblast` and `gsea_dotplot_dc` are
the four boxes of printed panel N, in that order, and `gsea_legend` is its
fifth. Each is drawn at its own printed sub-box, read from
03_Final_Panels/slot_subrects.csv through 00_Config/slots.py - not at a fifth
of the panel rect, which is taller than the five boxes.

The two input tables, the eight pathways, the four cell types, the per-cell-type
NES sort, the global symmetric colour limits, the PuOr_r colormap and the three
significance tiers are unchanged.

EVERY DOTPLOT KEEPS ITS OWN PATHWAY LABELS. The eight pathways are sorted by
NES within each cell type, so the four boxes print them in four different
orders - MoMac runs TNF-alpha, EMT, inflammatory, hypoxia; fibroblast runs
TNF-alpha, inflammatory, angiogenesis, IL-6/JAK/STAT3. One shared set of tick
labels would name the wrong dot in three boxes out of four.

`run_gsea_5celltypes.py` and `run_gsea_momac.py` next door are ANALYSIS
scripts: they compute the tables this one reads. They are not panel scripts and
are not restyled.

MARK
    The earlier drawing used a canvas four times the printed size and its
    smallest body type is the italic "Positive = NR enriched" note in the
    legend, set at 4 * SCALE, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The three dot areas, the dot edge width, the NES = 0 rule and the y grid
    are scaled by those. Tick widths and lengths and spine widths are not:
    those are style, and cnsplots sets them. Both are recomputed once the
    style is applied, because `tick_pt` is only settled then.

Judgement calls, stated plainly:
  - `fontweight='bold'` on the two legend headings is dropped; cnsplots bolds
    axis titles and panel letters and nothing else. They keep their emphasis by
    size, `body_pt` against `tick_pt` for the items, which is the ordering the
    earlier drawing had at 5.5 pt against 4.5 pt.
  - The dotplot title's explicit `fontweight='normal'` is dropped so cnsplots'
    bold axis title applies.
  - The legend's `ax.text` positions and the colorbar's `add_axes` rectangle
    are axes and figure fractions, i.e. layout; every string is unchanged.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
import panel_style_cns as style
import slots

BASE_DIR = Path(__file__).resolve().parent
# The two GSEA tables are inputs, not outputs: this script reads them and
# writes nothing back. A second copy of a number the manuscript quotes is
# exactly the kind of drift PROVENANCE.csv exists to prevent, so they are read
# where they already live, beside this script.
DATA_DIR = BASE_DIR / 'gsea_data'
PANEL_DIR = BASE_DIR

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 4.0                  # the earlier smallest body type, before * SCALE
MARK = style.tick_pt() / (SMALL_PT * SCALE)   # recomputed in main()
AREA = MARK ** 2

PANEL_LETTER = "N"
# The four dotplots are sub-boxes 1 to 4 of the printed panel, in the order
# CELL_TYPES lists them; the legend strip is the fifth.
SUB_OF = {'MoMac': 1, 'Epithelial': 2, 'Fibroblast': 3, 'DC': 4}
LEGEND_SUB = 5

# The panel letter is measured from the corner of the whole panel rect, which
# starts above and to the left of the first box; only the part of the keep-out
# that falls inside that box has to be reserved.
_rect = slots.rect_mm(5, PANEL_LETTER)
_sub1 = slots.rect_mm(5, PANEL_LETTER, sub=1)
_cw, _ch = slots.letter_cell_mm(5, PANEL_LETTER)
LETTER_CELL = (max(0.0, _cw - (_sub1[0] - _rect[0])),
               max(0.0, _ch - (_sub1[1] - _rect[1])))

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

# The eight Hallmark sets in their conventional short forms, the same forms
# the pathway UMAPs of this figure print. Each is one line: a label set at 45
# degrees puts its second line across its neighbour's first, so a two-line name
# collides whatever its length.
PATHWAY_LABELS = [
    'TNFα/NF-κB',
    'Inflammation',
    'IL-6/STAT3',
    'IFN-γ',
    'EMT',
    'Angiogenesis',
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
    global MARK, AREA
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

    print("=" * 60)
    print("Panel N: four GSEA dotplots and their shared legend")
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

        panel_w_mm, panel_h_mm = slots.size_mm(5, PANEL_LETTER,
                                               sub=SUB_OF[ct])
        fig, ax = style.subplots_mm(panel_w_mm, panel_h_mm)

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
        # Steeper than 45 degrees. Eight ticks across this box are about
        # 3.1 mm apart, and two labels set at 45 degrees leave only
        # 3.1/sqrt(2) = 2.2 mm between their baselines - less than the
        # 2.3 mm a 6 pt line occupies, so they touch whatever their
        # length. At 60 degrees the gap is 2.7 mm.
        ax.set_xticklabels(sorted_labels, rotation=60, ha='right',
                           rotation_mode='anchor')
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

        first = SUB_OF[ct] == 1
        style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL,
                          reserve_letter=first)
        over = style.overflow_mm(fig)
        if max(over) > 0:
            raise RuntimeError(
                f"ink outside the {panel_w_mm} x {panel_h_mm} mm canvas "
                f"(l,r,b,t mm): {over}")
        intruders = style.letter_clear(fig, LETTER_CELL) if first else []
        if intruders:
            raise RuntimeError(f"ink under the panel letter cell: {intruders}")

        stem = f'gsea_dotplot_{ct.lower().replace("/", "_")}'
        style.save_panel(fig, PANEL_DIR / stem)
        print(f"    Saved: {stem}.[svg|pdf|png] at "
              f"{panel_w_mm} x {panel_h_mm} mm")

    # Also generate shared legend panel
    _make_legend(cmap, vmin, vmax)


def _make_legend(cmap, vmin, vmax):
    """Create a small standalone legend panel with size + colorbar."""
    legend_w_mm, legend_h_mm = slots.size_mm(5, PANEL_LETTER,
                                             sub=LEGEND_SUB)
    fig, ax = style.subplots_mm(legend_w_mm, legend_h_mm)
    ax.set_axis_off()
    # The strip is 12 mm wide, so all of it is used: the default subplot
    # margins would leave the entries about 9 mm to be set in. Entries are
    # stacked rather than set beside one another for the same reason.
    ax.set_position([0.0, 0.0, 1.0, 1.0])

    # Size legend. Every marker stays exactly where it was placed: the three
    # dots are drawn artists, and moving one would be a content change. Only
    # the strings beside them are set closer, which is what makes them fit.
    y_start = 0.92
    ax.text(0.05, y_start + 0.05, 'P value', fontsize=style.body_pt(),
            transform=ax.transAxes, va='top')
    for i, (label, size) in enumerate([
        ('p ≤ 0.05', 150 * SCALE * AREA),
        ('p ≤ 0.10', 100 * SCALE * AREA),
        ('p > 0.10', 60 * SCALE * AREA),
    ]):
        y = y_start - 0.12 * (i + 1)
        ax.scatter(0.15, y, s=size, c='gray', alpha=0.6, edgecolors='black',
                   linewidths=0.5 * MARK, transform=ax.transAxes, zorder=3)
        ax.text(0.30, y, label, fontsize=style.tick_pt(), va='center',
                transform=ax.transAxes)

    # Colour legend. The heading names the scale, so the colorbar is not
    # labelled a second time underneath it.
    ax.text(0.05, 0.38, 'NES', fontsize=style.body_pt(),
            transform=ax.transAxes, va='top')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cax = fig.add_axes([0.15, 0.15, 0.6, 0.04])
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')

    # The legend strip places its items itself, in axes fractions, and its
    # colorbar sits in an add_axes rectangle that subplots_adjust does not
    # move, so the margins are not fitted here; overflow still has to be zero.
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {legend_w_mm} x {legend_h_mm} mm legend strip "
            f"(l,r,b,t mm): {over}")
    stem = 'gsea_legend'
    style.save_panel(fig, PANEL_DIR / stem)
    print(f"\n    Saved: {stem}.[svg|pdf|png] at "
          f"{legend_w_mm} x {legend_h_mm} mm")


if __name__ == '__main__':
    main()
