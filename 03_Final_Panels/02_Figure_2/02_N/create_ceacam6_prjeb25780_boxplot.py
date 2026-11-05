#!/usr/bin/env python3
"""
Figure 2L, left violin: BayesPrism epithelial-deconvolved CEACAM6 expression
in the PRJEB25780 bulk cohort, responders versus non-responders. draw() is
shared with the right violin (02_O/create_ceacam5_prjeb25780_boxplot.py,
CEACAM5).

Violins with an inner box and the per-sample points. The panel prints the
exact one-tailed Mann-Whitney P (non-responders > responders); the two-sided
P is reported in the legend.

Reads TIGER_BAYESPRISM_EPI and TIGER_META (00_Config/paths.py); normal
samples are excluded and expression is log2(x + 1).
Writes ceacam6_prjeb25780_boxplot.[svg|pdf|png] beside this script.

Run: python create_ceacam6_prjeb25780_boxplot.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TIGER_BAYESPRISM_EPI, TIGER_META       # noqa: E402
import panel_style_cns as style                          # noqa: E402
import slots                                             # noqa: E402
from cnsfig.boxes import box_xlim, frame_fixed, printed_and_two_sided  # noqa: E402

# Shared response/phase colours.
from shared.figure_config import GROUP_COLORS, PRE_PAIR, POST_PAIR  # noqa: E402
COLOR_R, COLOR_NR = PRE_PAIR          # pre-treatment cohort: the light pair

PANEL_LETTER = "L"


def draw(out_dir, sub, gene, stem, ylabel, left_mm, letter_cell=None):
    """One gene's R-versus-NR violin in sub-box `sub` of Figure 2L.

    The four plots of Figure 2K and 2L share fixed mm margins (frame_fixed;
    left 9.0 mm with a y label, 5.4 mm without), so their frames line up.
    printed_and_two_sided returns the one-tailed P the panel prints (as an
    exact value) and the two-sided P for the legend from one call.
    """
    panel_w_mm, panel_h_mm = slots.size_mm(2, PANEL_LETTER, sub=sub)

    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)

    epi_expr = pd.read_csv(TIGER_BAYESPRISM_EPI, sep='\t', index_col=0)
    meta_df = pd.read_csv(TIGER_META, sep='\t')
    meta_df = meta_df[meta_df['Treatment'] != 'Normal']

    common = [s for s in meta_df['sample_id'] if s in epi_expr.index]
    meta_df = meta_df[meta_df['sample_id'].isin(common)].set_index('sample_id')
    epi_expr = epi_expr.loc[common]

    r_vals = np.log2(epi_expr.loc[meta_df['response_NR'] == 'R', gene].values + 1)
    nr_vals = np.log2(epi_expr.loc[meta_df['response_NR'] == 'N', gene].values + 1)

    pval, pval_two_sided = printed_and_two_sided(nr_vals, r_vals)

    fig, ax = style.subplots_mm(panel_w_mm, panel_h_mm)

    data = [r_vals, nr_vals]
    positions = [1, 2]
    colors = [COLOR_R, COLOR_NR]

    vp = ax.violinplot(data, positions=positions, widths=0.7,
                       showmeans=False, showmedians=False, showextrema=False)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors[i])
        body.set_edgecolor('black')
        body.set_linewidth(style.EDGE_PT)
        body.set_alpha(1.0)

    bp = ax.boxplot(data, positions=positions, widths=0.15, patch_artist=True,
                    showfliers=False,
                    boxprops=dict(facecolor='white', linewidth=style.RULE_PT),
                    whiskerprops=dict(color='black', linewidth=style.RULE_PT),
                    capprops=dict(color='black', linewidth=style.RULE_PT),
                    medianprops=dict(color='black', linewidth=style.RULE_PT))

    np.random.seed(42)
    for i, (pos, vals, color) in enumerate(zip(positions, data, colors)):
        jitter = np.random.uniform(-0.08, 0.08, len(vals))
        ax.scatter(pos + jitter, vals, c=color,
                   s=(0.76 * style.PT_PER_MM) ** 2,  # 0.76 mm across
                   alpha=0.8,
                   edgecolors='white', linewidths=style.EDGE_PT, zorder=3)

    y_max = max(np.max(r_vals), np.max(nr_vals))
    y_bracket = y_max * 1.1
    ax.plot([1, 1, 2, 2],
            [y_bracket, y_bracket * 1.03, y_bracket * 1.03, y_bracket],
            'k-', linewidth=style.RULE_PT)
    # Exact P, not a star.
    p_text, p_kw = style.p_text_kw(pval, exact=True)
    ax.text(1.5, y_bracket * 1.05, p_text, ha='center', va='bottom', **p_kw)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Pre-R', 'Pre-NR'])
    ax.set_xlim(*box_xlim([1, 2], 0.7, clear=0.25))   # narrow frame: room for the group names
    ax.set_ylim(0, y_max * 1.70)   # headroom: the P string clears the title
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # The y label fits on one line, so the K and L left boxes share a margin.
    frame_fixed(fig, ax, title=f'*{gene}* (Epi)', ylabel=ylabel,
                left_mm=left_mm, right_mm=0.6, bottom_mm=3.8,
                panel_w_mm=panel_w_mm, panel_h_mm=panel_h_mm,
                letter_cell=letter_cell)

    return style.save_panel(fig, Path(out_dir) / stem)


def main():
    draw(Path(__file__).parent, 1, 'CEACAM6', 'ceacam6_prjeb25780_boxplot',
         'Expression (log2)', 9.0,
         letter_cell=slots.letter_cell_mm(2, PANEL_LETTER))


if __name__ == '__main__':
    main()

