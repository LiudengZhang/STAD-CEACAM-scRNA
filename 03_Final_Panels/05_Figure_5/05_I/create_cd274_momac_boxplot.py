#!/usr/bin/env python3
"""Figure 5J, box 1 of 4: CD274 (PD-L1) expression in monocytes/macrophages,
post-treatment responders against non-responders (stomach samples with at
least 20 cells, one mean per sample). The panel prints the exact one-sided
Mann-Whitney P (NR > R); the two-sided value is computed for the legend.

draw() here is shared by the other three boxes, whose entry points are
05_J/create_cd274_epithelial_boxplot.py, 05_K/create_cd274_fibroblast_boxplot.py
and 05_DC_CD274/create_cd274_dc_boxplot.py. Each box keeps its own y range;
the y axis is labelled once, on the leftmost box.

Input:  MOMAC_H5AD (00_Config/paths.py)
Output: cd274_momac_boxplot.{svg,pdf,png} in this directory
Run:    python create_cd274_momac_boxplot.py
"""
import warnings, numpy as np, scanpy as sc
import matplotlib; matplotlib.use('Agg')
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD
import panel_style_cns as style
import slots
from cnsfig.boxes import finish_two_group, draw_boxes, bracket, ylim_above, assert_no_points, printed_and_two_sided

warnings.filterwarnings('ignore')

H5AD = MOMAC_H5AD
OUT_DIR = Path(__file__).parent

PANEL_LETTER = "J"
SUB = 1

# Response/phase colours shared by all figures.
from shared.figure_config import GROUP_COLORS, PRE_PAIR, POST_PAIR  # noqa: E402
BOX_COLORS = {'R': POST_PAIR[0], 'NR': POST_PAIR[1]}   # post-treatment
MIN_CELLS = 20


def draw(h5ad, sub, title, ylabel_markup, top_extra_mm, stem, out_dir):
    """Draw one box of panel J into sub-box `sub`.

    `title` is the cell type, `ylabel_markup` the shared y-axis label (empty
    except on the leftmost box), `top_extra_mm` the room the title needs, and
    `out_dir` the calling script's directory.
    """
    panel_w_mm, panel_h_mm = slots.size_mm(5, PANEL_LETTER, sub=sub)

    # Reserve only the part of the panel-letter cell that falls inside this box.
    _rect = slots.rect_mm(5, PANEL_LETTER)
    _sub = slots.rect_mm(5, PANEL_LETTER, sub=sub)
    _cw, _ch = slots.letter_cell_mm(5, PANEL_LETTER)
    letter_cell = (max(0.0, _cw - (_sub[0] - _rect[0])),
                   max(0.0, _ch - (_sub[1] - _rect[1])))

    np.random.seed(0)

    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)

    adata = sc.read_h5ad(h5ad)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata = adata[adata.obs['Treatment phase'] == 'Post'].copy()
    adata.obs['response'] = adata.obs['stomach_post_grouping'].map({
        'Responsed': 'R', 'No-response': 'NR'}).astype(str)
    adata = adata[adata.obs['response'].isin(['R', 'NR'])].copy()

    gene = 'CD274'
    if adata.raw:
        idx = list(adata.raw.var_names).index(gene)
        expr = adata.raw.X[:, idx]
    else:
        idx = list(adata.var_names).index(gene)
        expr = adata.X[:, idx]
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    adata.obs['value'] = np.array(expr).flatten()

    df = adata.obs[['sample', 'response', 'value']].copy()
    df['sample'] = df['sample'].astype(str)
    counts = df.groupby('sample', observed=True).size()
    valid = counts[counts >= MIN_CELLS].index
    df = df[df['sample'].isin(valid)]
    sample_df = df.groupby(['sample', 'response'], observed=True)['value'].mean().reset_index()

    r_data = sample_df[sample_df['response'] == 'R']['value'].values
    nr_data = sample_df[sample_df['response'] == 'NR']['value'].values

    fig, ax = style.subplots_mm(panel_w_mm, panel_h_mm)
    bp = draw_boxes(ax, [r_data, nr_data], [0, 1],
                    [BOX_COLORS['R'], BOX_COLORS['NR']], width=0.5)

    # One-sided P (NR > R) for the panel, two-sided for the legend.
    pval, pval_two_sided = printed_and_two_sided(nr_data, r_data)
    y_max = max(np.max(r_data), np.max(nr_data))
    y_range = y_max - min(np.min(r_data), np.min(nr_data))
    bracket_y = y_max + 0.10 * y_range
    _, p_text, _ = bracket(fig, ax, 0, 1, y_max, y_range, pval, kind="pair",
                           lift=0.08, arm=0.02, exact=True)
    ax.set_ylim(top=bracket_y + 0.42 * y_range)
    ylim_above(ax, p_text)
    assert_no_points(ax, bp)

    # Shared frame layout for the boxes of panels J and L.
    finish_two_group(
        fig, ax, title=title,
        tick_labels=['Post-R', 'Post-NR'], tick_rotation=45,
        positions=[0, 1], width=0.5, ylabel_markup=ylabel_markup,
        letter_cell=letter_cell, panel_w_mm=panel_w_mm, panel_h_mm=panel_h_mm,
        top_extra_mm=top_extra_mm, bottom_mm=8.60)
    style.save_panel(fig, Path(out_dir) / stem)
    print(f"  Saved: {stem}.[svg|pdf|png] at {panel_w_mm} x {panel_h_mm} mm; "
          f"R n={len(r_data)}, NR n={len(nr_data)}, P={pval:.4f}")


def main():
    draw(H5AD, SUB, 'Monocytes/\nMacrophages',
         'PD-L1 (*CD274*) expression', 0.27, 'cd274_momac_boxplot', OUT_DIR)


if __name__ == '__main__':
    main()
