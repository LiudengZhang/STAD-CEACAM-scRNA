#!/usr/bin/env python3
"""
Figure 5 panel L - IL-6/JAK/STAT3 signalling score in CD4+ T cells, responders
versus non-responders after treatment, one point per sample.

The h5ad, the stomach/post filter, the R-versus-NR mapping, the twenty-cell
minimum, the sample-level mean, the two-sided Mann-Whitney test and the bracket
geometry are unchanged. `np.random.seed(0)` stays exactly where it was.

  printed panel  Figure 5 L       (PROVENANCE.csv - the directory is
                                   "05_IL6_CD4"; the letter was looked up)

The panel is 18.5 mm wide. The gene set is named in full in the axis label, on
one line rather than two, and the title carries the cell type alone, as the box
titles of panel J do. Both are declared in 00_Config/shared/labels.py.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the tick labels - at 5 * SCALE, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    Box, whisker, cap and median line widths, the bracket width and
    the jitter marker area and edge width are scaled by those. Tick widths and lengths and spine widths are
    not: those are style, and cnsplots sets them.

Judgement calls, stated plainly:
  - `fontweight='bold'` on the significant-P annotation is dropped; cnsplots
    bolds axis titles and panel letters and nothing else. The size distinction
    between a significant and a non-significant P is kept, taken from the
    system: body_pt versus tick_pt.
  - The axes title's explicit `fontweight='normal'` is dropped so cnsplots'
    bold axis title applies, for the same reason.
"""
import warnings, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD4_H5AD
import panel_style_cns as style
import slots
from cnsfig.boxes import finish_two_group, draw_boxes, bracket, ylim_above, assert_no_points

warnings.filterwarnings('ignore')

OUT_DIR = Path(__file__).parent
SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier smallest body type, before * SCALE

PANEL_LETTER = "L"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(5, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(5, PANEL_LETTER)
#: The letter cell, then two lines of body type, given up by the axes.

BOX_COLORS = {'R': '#2166AC', 'NR': '#B2182B'}
MIN_CELLS = 20

# IL-6/JAK/STAT3 Signaling gene list from MSigDB Hallmark
PATHWAY_GENES = [
    'INHBE','IL17RA','IRF9','IL17RB','MAP3K8','CCR1','FAS','CXCL3','A2M','CD38',
    'SOCS3','TYK2','GRB2','CXCL13','TNFRSF1B','CXCL1','CBL','PF4','CSF1','IFNGR1',
    'HMOX1','TNF','HAX1','IL12RB1','CSF2','IL2RG','JUN','ITGA4','IL18R1','IL6',
    'MYD88','CXCL11','LEPR','LTB','PDGFC','PTPN11','IFNAR1','DNTT','IL1B','SOCS1',
    'TNFRSF12A','PIK3R5','IL2RA','CSF2RA','STAT3','IL13RA1','BAK1','TLR2','CRLF2',
    'CXCL9','PIM1','TNFRSF21','PTPN2','OSMR','CSF3R','IL4R','IL6ST','STAM2','CSF2RB',
    'EBI3','STAT2','TNFRSF1A','IL1R2','STAT1','CCL7','CD14','TGFB1','IRF1','IL3RA',
    'IL10RB','IL1R1','CD44','ITGB3','ACVRL1','CXCL10','IL15RA','CNTFR','PLA2G2A',
    'ACVR1B','IL9R','LTBR','CD9','IFNGR2','PTPN1','CD36','REG1A','IL7',
]


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("=" * 60)
    print("Panel O: IL-6/JAK/STAT3 — CD4+ T cells")
    print("=" * 60)

    adata = sc.read_h5ad(TCD4_H5AD)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata = adata[adata.obs['Treatment phase'] == 'Post'].copy()
    adata.obs['response'] = adata.obs['stomach_post_grouping'].map({
        'Responsed': 'R', 'No-response': 'NR'}).astype(str)
    adata = adata[adata.obs['response'].isin(['R', 'NR'])].copy()
    print(f"  Post-treatment stomach CD4+ T cells: {adata.n_obs}")

    gene_names = list(adata.raw.var_names) if adata.raw else list(adata.var_names)
    avail = [g for g in PATHWAY_GENES if g in gene_names]
    print(f"  Genes available: {len(avail)}/{len(PATHWAY_GENES)}")
    # score_genes with use_raw=True reaches straight into adata.raw.var_names,
    # so it has to be told when the file carries no .raw - the clean deposit
    # holds the same log1p matrix in .X.
    sc.tl.score_genes(adata, gene_list=avail, score_name='value',
                     ctrl_size=min(50, len(avail)),
                     use_raw=adata.raw is not None)

    # Sample-level aggregation
    df = adata.obs[['sample', 'response', 'value']].copy()
    df['sample'] = df['sample'].astype(str)
    counts = df.groupby('sample', observed=True).size()
    valid = counts[counts >= MIN_CELLS].index
    df = df[df['sample'].isin(valid)]
    sample_df = df.groupby(['sample', 'response'], observed=True)['value'].mean().reset_index()

    r_data = sample_df[sample_df['response'] == 'R']['value'].values
    nr_data = sample_df[sample_df['response'] == 'NR']['value'].values
    print(f"  R samples: n={len(r_data)}, NR samples: n={len(nr_data)}")
    print(f"  R mean: {np.mean(r_data):.4f}, NR mean: {np.mean(nr_data):.4f}")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    # ONE BOX, NO POINTS (2026-09-16, the author's fifth reading: "box
    # plots should look alike throughout; don't show every point"). The
    # box is cnsfig.boxes.draw_boxes (the 2H/I box, 0.79 mm open fliers);
    # the jittered points (default_rng(42)) are no longer drawn - the
    # author's ruling, a declared departure from the published panel.
    bp = draw_boxes(ax, [r_data, nr_data], [0, 1],
                    [BOX_COLORS['R'], BOX_COLORS['NR']], width=0.5)
    bp['boxes'][0].set_alpha(0.6); bp['boxes'][1].set_alpha(0.6)

    _, pval = stats.mannwhitneyu(nr_data, r_data, alternative='two-sided')
    print(f"  P-value (two-sided): {pval:.4f}")

    y_max = max(np.max(r_data), np.max(nr_data))
    y_range = y_max - min(np.min(r_data), np.min(nr_data))
    # 0.22 of the range above the data, not 0.10: this panel is 18.5 mm wide,
    # its axes box about 10, and 'P = 0.030' sets 9.5 mm at 6 pt, so the
    # string overhangs the box on both sides. At 0.10 it sat level with the top
    # y tick and the sweep convicted it against '0.22' at 0.55 mm2; lifted, it
    # clears the last tick label. The bracket and the P value do not change.
    bracket_y = y_max + 0.22 * y_range
    # Same vertices (top line at +0.22, arms 0.02 down), the label's ink
    # 0.4 mm above the line (cnsfig.boxes.bracket, 2026-09-16).
    _, p_text, _ = bracket(fig, ax, 0, 1, y_max, y_range, pval, kind="pair",
                           lift=0.20, arm=0.02)
    # Headroom (2026-09-15): the bracket and star sit lower in the frame, so
    # the star clears the title's second line.
    ax.set_ylim(top=bracket_y + 0.42 * y_range)
    ylim_above(ax, p_text)
    assert_no_points(ax, bp)

    # ONE FAMILY, ONE FRAME LINE  (2026-09-14, evening)
    #   The four boxes of J and the box of L are framed by
    #   cnsfig.boxes.finish_two_group: the tick labels with their counts, the
    #   y label as a rich run at the tick size, the title on the canvas in a
    #   two-line band, the frame's top at 139.8 mm and its bottom at 158.7
    #   mm on the page (sweep_pages.ROW_ALIGN), the boxes off the spines.
    finish_two_group(
        fig, ax, title='IL-6/JAK/STAT3\nCD4+ T cells',
        tick_labels=[f'R\n(n={len(r_data)})', f'NR\n(n={len(nr_data)})'],
        positions=[0, 1], width=0.5, ylabel_markup='IL-6/JAK/STAT3 score',
        letter_cell=LETTER_CELL, panel_w_mm=PANEL_W_MM, panel_h_mm=PANEL_H_MM,
        top_extra_mm=3.32, bottom_mm=6.30)
    print(f"  tick labels R, NR; counts R n={len(r_data)}, "
          f"NR n={len(nr_data)} -> figure legend")
    stem = 'il6_stat3_cd4t_boxplot'
    style.save_panel(fig, OUT_DIR / stem)
    print(f"  Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
