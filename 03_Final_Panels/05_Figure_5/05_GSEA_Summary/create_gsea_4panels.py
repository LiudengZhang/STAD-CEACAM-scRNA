#!/usr/bin/env python3
"""Figure 5N: Hallmark GSEA summary, one dot plot per cell type (MoMac,
epithelial, fibroblast, DC) plus a shared legend strip.

Each dot plot sorts the eight pathways by its own NES. Colour is NES on one
symmetric scale (PuOr_r) across all four; dot size marks nominal P
(<= 0.05, <= 0.10, > 0.10).

Input:  gsea_data/gsea_combined_5types.csv (run_gsea_5celltypes.py) and
        gsea_data/gsea_momac.csv (run_gsea_momac.py)
Output: gsea_dotplot_{momac,epithelial,fibroblast,dc}.{svg,pdf,png} and
        gsea_legend.{svg,pdf,png} in this directory
Run:    python create_gsea_4panels.py
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
from cnsfig.layout import pin_frame_mm

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / 'gsea_data'
PANEL_DIR = BASE_DIR

SCALE = 4                       # title-pad scaling constants
SMALL_PT = 4.0
MARK = style.tick_pt() / (SMALL_PT * SCALE)   # recomputed in main()

PANEL_LETTER = "N"
# Sub-boxes 1-4 of the panel are the dot plots; 5 is the legend strip.
SUB_OF = {'MoMac': 1, 'Epithelial': 2, 'Fibroblast': 3, 'DC': 4}
LEGEND_SUB = 5

# The first box reserves only the part of the panel-letter cell inside it.
_rect = slots.rect_mm(5, PANEL_LETTER)
_sub1 = slots.rect_mm(5, PANEL_LETTER, sub=1)
_cw, _ch = slots.letter_cell_mm(5, PANEL_LETTER)
LETTER_CELL = (max(0.0, _cw - (_sub1[0] - _rect[0])),
               max(0.0, _ch - (_sub1[1] - _rect[1])))

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

#: Short pathway names, as in panel B.
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

CELL_TYPES = ['MoMac', 'Epithelial', 'Fibroblast', 'DC']
CELL_LABELS = {
    'MoMac': 'Monocytes/\nMacrophages',
    'Epithelial': 'Epithelial',
    'Fibroblast': 'Fibroblast',
    'DC': 'DC',
}


#: Largest dot 1.56 mm across; the three P classes have areas 150 : 100 : 30.
DOT_S_MAX = (1.56 * style.PT_PER_MM) ** 2


def assign_dot_size(p_value):
    if p_value <= 0.05:
        return DOT_S_MAX
    elif p_value <= 0.10:
        return DOT_S_MAX * 100.0 / 150.0
    else:
        return DOT_S_MAX * 30.0 / 150.0


def main():
    global MARK
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)

    existing = pd.read_csv(DATA_DIR / 'gsea_combined_5types.csv')
    momac = pd.read_csv(DATA_DIR / 'gsea_momac.csv')
    df = pd.concat([existing, momac], ignore_index=True)
    df['NES'] = df['NES'].astype(float)
    df['NOM_pval'] = df['NOM_pval'].astype(float)

    # One symmetric NES colour range for all four plots
    subset = df[df['CellType'].isin(CELL_TYPES) & df['Pathway'].isin(PATHWAYS)]
    nes_abs_max = max(abs(subset['NES'].min()), abs(subset['NES'].max()))
    vmin, vmax = -nes_abs_max, nes_abs_max

    cmap = plt.cm.PuOr_r  # warm orange = positive NES (enriched in NR)

    LABEL_MAP = dict(zip(PATHWAYS, PATHWAY_LABELS))

    for ct in CELL_TYPES:
        ct_df = df[df['CellType'] == ct]

        rows = []
        for pw in PATHWAYS:
            row = ct_df[ct_df['Pathway'] == pw]
            if len(row) > 0:
                r = row.iloc[0]
                rows.append({'pathway': pw, 'nes': float(r['NES']), 'pval': float(r['NOM_pval'])})
            else:
                rows.append({'pathway': pw, 'nes': 0.0, 'pval': 1.0})

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
            edgecolors='black', linewidths=style.EDGE_PT, vmin=vmin, vmax=vmax, zorder=3
        )

        ax.axhline(y=0, color='gray', linestyle='--', linewidth=style.RULE_PT, alpha=0.7)

        ax.set_xticks(x_positions)
        ax.set_xticklabels(sorted_labels, rotation=60, ha='right',
                           rotation_mode='anchor')
        ax.set_xlim(-0.8, len(rows) - 0.2)

        ax.set_ylim(0.7, 2.25)
        ax.set_ylabel('NES')

        ax.set_title(CELL_LABELS[ct], pad=8 * MARK)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        first = SUB_OF[ct] == 1
        style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL,
                          reserve_letter=first)
        # Align the four frame tops.
        _sy0 = slots.rect_mm(5, PANEL_LETTER, sub=SUB_OF[ct])[1]
        pin_frame_mm(fig, ax, top_mm=178.2 - _sy0,
                     bottom_mm=ax.get_position().y0 * panel_h_mm)
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

    _make_legend(cmap, vmin, vmax)
    print(f"  5 drawings saved: 4 dotplots at their sub-boxes + the legend strip")


def _make_legend(cmap, vmin, vmax):
    """Create a small standalone legend panel with size + colorbar."""
    legend_w_mm, legend_h_mm = slots.size_mm(5, PANEL_LETTER,
                                             sub=LEGEND_SUB)
    fig, ax = style.subplots_mm(legend_w_mm, legend_h_mm)
    ax.set_axis_off()
    ax.set_position([0.0, 0.0, 1.0, 1.0])

    y_start = 0.92
    ax.text(0.00, y_start + 0.05, 'Significance', fontsize=style.tick_pt(),
            transform=ax.transAxes, va='top')
    for i, (label, size) in enumerate([
        ('P ≤ 0.05', assign_dot_size(0.05)),
        ('P ≤ 0.10', assign_dot_size(0.10)),
        ('P > 0.10', assign_dot_size(0.50)),
    ]):
        y = y_start - 0.12 * (i + 1)
        ax.scatter(0.15, y, s=size, c='gray', alpha=0.6, edgecolors='black',
                   linewidths=style.EDGE_PT, transform=ax.transAxes, zorder=3)
        ax.text(0.30, y, label, fontsize=style.tick_pt(), va='center',
                transform=ax.transAxes)

    ax.text(0.00, 0.38, 'Enrichment', fontsize=style.tick_pt(),
            transform=ax.transAxes, va='top')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cax = fig.add_axes([0.05, 0.15, 0.90, 0.04])
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal', ticks=[-1, 0, 1])
    cax.tick_params(labelsize=style.tick_pt(), width=style.RULE_PT, length=1.5, pad=1.0)
    cbar.outline.set_linewidth(style.RULE_PT)

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {legend_w_mm} x {legend_h_mm} mm legend strip "
            f"(l,r,b,t mm): {over}")
    style.save_panel(fig, PANEL_DIR / 'gsea_legend')


if __name__ == '__main__':
    main()
