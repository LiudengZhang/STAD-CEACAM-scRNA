#!/usr/bin/env python3
"""Figure 5L: Hallmark IL-6/JAK/STAT3 score (sc.tl.score_genes) in CD4+ T
cells, post-treatment responders against non-responders (stomach samples with
at least 20 cells, one mean per sample). The panel prints the exact one-sided
Mann-Whitney P (NR > R); the two-sided value is computed for the legend.

Input:  TCD4_H5AD (00_Config/paths.py)
Output: il6_stat3_cd4t_boxplot.{svg,pdf,png} in this directory
Run:    python create_il6_stat3_cd4t_boxplot.py
"""
import warnings, numpy as np, scanpy as sc
import matplotlib; matplotlib.use('Agg')
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD4_H5AD
from cnsfig import Panel, style
from cnsfig.boxes import finish_two_group, draw_boxes, bracket, ylim_above, assert_no_points, printed_and_two_sided

warnings.filterwarnings('ignore')

PANEL = Panel(5, "L", __file__, small_pt=5.0)

# Response/phase colours shared by all figures.
from shared.figure_config import GROUP_COLORS, PRE_PAIR, POST_PAIR  # noqa: E402
BOX_COLORS = {'R': POST_PAIR[0], 'NR': POST_PAIR[1]}   # post-treatment
MIN_CELLS = 20

# MSigDB Hallmark IL-6/JAK/STAT3 signaling genes
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
    PANEL.apply()

    adata = sc.read_h5ad(TCD4_H5AD)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata = adata[adata.obs['Treatment phase'] == 'Post'].copy()
    adata.obs['response'] = adata.obs['stomach_post_grouping'].map({
        'Responsed': 'R', 'No-response': 'NR'}).astype(str)
    adata = adata[adata.obs['response'].isin(['R', 'NR'])].copy()

    gene_names = list(adata.raw.var_names) if adata.raw else list(adata.var_names)
    avail = [g for g in PATHWAY_GENES if g in gene_names]
    # Use .raw when present, else .X.
    sc.tl.score_genes(adata, gene_list=avail, score_name='value',
                     ctrl_size=min(50, len(avail)),
                     use_raw=adata.raw is not None)

    df = adata.obs[['sample', 'response', 'value']].copy()
    df['sample'] = df['sample'].astype(str)
    counts = df.groupby('sample', observed=True).size()
    valid = counts[counts >= MIN_CELLS].index
    df = df[df['sample'].isin(valid)]
    sample_df = df.groupby(['sample', 'response'], observed=True)['value'].mean().reset_index()

    r_data = sample_df[sample_df['response'] == 'R']['value'].values
    nr_data = sample_df[sample_df['response'] == 'NR']['value'].values

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)
    bp = draw_boxes(ax, [r_data, nr_data], [0, 1],
                    [BOX_COLORS['R'], BOX_COLORS['NR']], width=0.5)

    # One-sided P (NR > R) for the panel, two-sided for the legend.
    pval, pval_two_sided = printed_and_two_sided(nr_data, r_data)

    y_max = max(np.max(r_data), np.max(nr_data))
    y_range = y_max - min(np.min(r_data), np.min(nr_data))
    # Bracket raised to keep the P label clear of the top y tick label.
    bracket_y = y_max + 0.22 * y_range
    _, p_text, _ = bracket(fig, ax, 0, 1, y_max, y_range, pval, kind="pair",
                           lift=0.20, arm=0.02, exact=True)
    ax.set_ylim(top=bracket_y + 0.42 * y_range)
    ylim_above(ax, p_text)
    assert_no_points(ax, bp)

    # Shared frame layout for the boxes of panels J and L.
    finish_two_group(
        fig, ax, title='IL-6/JAK/STAT3\nCD4+ T cells',
        tick_labels=['Post-R', 'Post-NR'], tick_rotation=45,
        positions=[0, 1], width=0.5, ylabel_markup='IL-6/JAK/STAT3 score',
        letter_cell=PANEL.letter_cell, panel_w_mm=PANEL.w_mm,
        panel_h_mm=PANEL.h_mm,
        top_extra_mm=3.32, bottom_mm=8.60)
    PANEL.finish(fig, 'il6_stat3_cd4t_boxplot')


if __name__ == '__main__':
    main()
