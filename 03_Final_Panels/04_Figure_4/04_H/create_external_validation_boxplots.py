#!/usr/bin/env python3
"""Figure 4H: MoMac_IL1B (Mac_3) proportion in tumour against normal tissue in
two external cohorts, GSE239676 (top) and GSE183904 (bottom), drawn as two
separate boxes. Each prints the exact one-sided Mann-Whitney P (tumour >
normal); the two-sided value is computed alongside for the legend.

Input:  FIG4_EXTERNAL_LINGHUA, FIG4_EXTERNAL_KUMAR (00_Config/paths.py)
Output: mac3_boxplot_linghua.{svg,pdf,png} and mac3_boxplot_kumar.{svg,pdf,png}
        in this directory
Run:    python create_external_validation_boxplots.py
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
from scipy.stats import mannwhitneyu  # noqa: F401
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                       # noqa: E402,F401,F403
import panel_style_cns as style           # noqa: E402
import slots                              # noqa: E402
from cnsfig.boxes import draw_boxes, bracket, ylim_above, assert_no_points, printed_and_two_sided  # noqa: E402
from cnsfig.layout import pin_frame_mm    # noqa: E402

PANEL_LETTER = "H"

# Each box reserves only the part of the panel-letter cell that falls inside it.
_rect = slots.rect_mm(4, PANEL_LETTER)
_cw, _ch = slots.letter_cell_mm(4, PANEL_LETTER)


#: Both frames have the same height, aligned with panel G (page mm).
FRAME_TOP_MM, FRAME_BOTTOM_MM = 159.0, 203.5
TICKS_MM, TITLE_MM = 3.92, 3.79
_s1, _s2 = slots.rect_mm(4, PANEL_LETTER, sub=1), slots.rect_mm(4, PANEL_LETTER, sub=2)
FRAME_H_MM = ((_s1[3] - TICKS_MM - FRAME_TOP_MM)
              + (FRAME_BOTTOM_MM - _s2[1] - TITLE_MM)) / 2.0


def letter_cell(sub):
    x0, y0, _x1, _y1 = slots.rect_mm(4, PANEL_LETTER, sub=sub)
    return (max(0.0, _cw - (x0 - _rect[0])), max(0.0, _ch - (y0 - _rect[1])))


COLOR_TUMOR = '#E377C2'
COLOR_NORMAL = '#17BECF'


def create_boxplot(data, title, stem, sub, reserve_letter,
                   y_label='MoMac_IL1B\nProportion'):
    """Create a single boxplot for T vs N comparison, at its printed sub-box."""

    w_mm, h_mm = slots.size_mm(4, PANEL_LETTER, sub=sub)
    cell = letter_cell(sub)

    fig, ax = style.subplots_mm(w_mm, h_mm)

    # One-sided P (tumour > normal) for the panel, two-sided for the legend.
    tumor_vals = data[data['Group'] == 'T']['Mac3_Proportion'].values
    normal_vals = data[data['Group'] == 'N']['Mac3_Proportion'].values
    p_value, p_two_sided = printed_and_two_sided(tumor_vals, normal_vals)

    bp = draw_boxes(ax, [tumor_vals, normal_vals], [0, 1],
                    [COLOR_TUMOR, COLOR_NORMAL], width=0.6)
    ax.set_xlim(-0.5, 1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['T', 'N'])

    y_max = data['Mac3_Proportion'].max()
    y_range = data['Mac3_Proportion'].max() - data['Mac3_Proportion'].min()
    if y_range == 0:
        y_range = 0.1
    bracket_top = y_max + y_range * 0.18
    _, p_txt, _ = bracket(fig, ax, 0, 1, y_max, y_range, p_value, kind="pair",
                          lift=0.10, arm=0.08, exact=True)

    ax.set_title(title)
    ax.set_ylabel(y_label)
    ax.set_xlabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(ax.get_ylim()[0], bracket_top + y_range * 0.25)
    ylim_above(ax, p_txt)
    assert_no_points(ax, bp)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=cell,
                      reserve_letter=reserve_letter)
    _sy0, _sy1 = slots.rect_mm(4, PANEL_LETTER, sub=sub)[1::2]
    if sub == 1:
        top = FRAME_TOP_MM
    else:
        top = FRAME_BOTTOM_MM - FRAME_H_MM
    pin_frame_mm(fig, ax, top_mm=top - _sy0,
                 bottom_mm=_sy1 - (top + FRAME_H_MM))
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {w_mm} x {h_mm} mm canvas "
            f"(l,r,b,t mm): {over}")
    if reserve_letter:
        intruders = style.letter_clear(fig, cell)
        if intruders:
            raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, stem)
    print(f"  Saved: {stem.name}.[svg|pdf|png] at box {sub}, "
          f"{w_mm} x {h_mm} mm; T n={len(tumor_vals)}, N n={len(normal_vals)}, "
          f"P={p_value:.4e}")


def main():
    script_dir = Path(__file__).parent

    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)

    df_linghua = pd.read_csv(FIG4_EXTERNAL_LINGHUA)

    # Primary tumour vs adjacent normal (peritoneal samples excluded)
    df_linghua_filtered = df_linghua[df_linghua['Tissue'].isin(['Primary_tumor', 'Adjacent_normal'])].copy()
    df_linghua_filtered['Group'] = df_linghua_filtered['Tissue'].map({
        'Primary_tumor': 'T',
        'Adjacent_normal': 'N'
    })

    create_boxplot(df_linghua_filtered, 'GSE239676',
                   script_dir / 'mac3_boxplot_linghua', sub=1,
                   reserve_letter=True)

    df_kumar = pd.read_csv(FIG4_EXTERNAL_KUMAR)

    df_kumar['Group'] = df_kumar['Tissue'].map({
        'Primary_tumor': 'T',
        'Normal': 'N'
    })

    create_boxplot(df_kumar, 'GSE183904',
                   script_dir / 'mac3_boxplot_kumar', sub=2,
                   reserve_letter=False)


if __name__ == "__main__":
    main()
