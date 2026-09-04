#!/usr/bin/env python3
# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:
# the test is two-sided and the annotation reports the exact P value.
# Original one-tailed version: Round_5/03_Final_Panels/05_Figure_5/05_IL6_CD4/create_il6_stat3_cd4t_boxplot.py
"""
Figure 5 panel L, RESTYLED (Version B) - IL-6/JAK/STAT3 signalling score in
CD4+ T cells, Post-R vs Post-NR.

Version A is
`03_Revised_Panels/Main_Figures/05_Figure_5/05_IL6_CD4/create_il6_stat3_cd4t_boxplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK
(areas by AREA), margins are millimetres, and the save is `style.save_panel`.

The 87-gene Hallmark list, the >= 20-cell sample filter, the sample-level
aggregation, the scoring call and the two-sided Mann-Whitney are Version A's,
unchanged.

  printed panel  Figure 5 L        (PROVENANCE.csv; the directory is not
                 letter-named, so this one could not have been inferred either
                 way - it is looked up, like all the others)
  printed rect   18.5 x 28.6 mm    (panel_rects.csv - the smallest panel in the
                 paper, printed at ~2.8 pt type)
  Version B box  44.0 x 48.0 mm

    18.5 mm cannot carry a two-line y axis label, two-line tick labels and a
    two-line title at 7/8 pt: about 14 mm of that width would be type before a
    single box is drawn. The growth is the point of the rebuild, not a
    deviation from it. Nothing plotted moves.

MARK
    Version A drew 3.5 x 5 cm at SCALE = 4 with its smallest body type at
    `5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2                                     = 0.1225

Judgement calls: the bold on a significant P is dropped, per PANEL_SPEC, and
its size emphasis taken from the system (body_pt vs tick_pt); the axes title's
explicit `fontweight='normal'` is dropped so cnsplots' bold axis title applies.
"""
import warnings, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import TCD4_H5AD
import panel_style_cns as style

warnings.filterwarnings('ignore')

OUT_DIR = Path(__file__).parent
SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 5.0                  # Version A's smallest body type, before * SCALE

PRINTED_MM = (18.5, 28.6)
PANEL_W_MM, PANEL_H_MM = 44.0, 48.0
MARGIN = dict(left=14.0, right=3.0, top=10.0, bottom=9.0)

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
    family = style.apply()
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
    bp = ax.boxplot([r_data, nr_data], positions=[0, 1], widths=0.5,
                   patch_artist=True, showfliers=False,
                   medianprops=dict(color='black', linewidth=1.5 * SCALE * MARK),
                   whiskerprops=dict(linewidth=1.0 * SCALE * MARK),
                   capprops=dict(linewidth=1.0 * SCALE * MARK),
                   boxprops=dict(linewidth=1.0 * SCALE * MARK))
    bp['boxes'][0].set_facecolor(BOX_COLORS['R']); bp['boxes'][0].set_alpha(0.6)
    bp['boxes'][1].set_facecolor(BOX_COLORS['NR']); bp['boxes'][1].set_alpha(0.6)

    rng = np.random.default_rng(42)
    for k, (data, color) in enumerate(zip([r_data, nr_data], [BOX_COLORS['R'], BOX_COLORS['NR']])):
        jitter = rng.uniform(-0.08, 0.08, len(data))
        ax.scatter([k] * len(data) + jitter, data, c=color, s=20 * SCALE * AREA,
                  edgecolors='white', linewidths=0.3 * SCALE * MARK, alpha=0.85, zorder=3)

    _, pval = stats.mannwhitneyu(nr_data, r_data, alternative='two-sided')
    print(f"  P-value (two-sided): {pval:.4f}")

    p_str = f'P = {pval:.3f}' if pval >= 0.001 else 'P < 0.001'
    is_star = pval < 0.05
    y_max = max(np.max(r_data), np.max(nr_data))
    y_range = y_max - min(np.min(r_data), np.min(nr_data))
    bracket_y = y_max + 0.10 * y_range
    ax.plot([0, 0, 1, 1], [bracket_y - 0.02 * y_range, bracket_y,
            bracket_y, bracket_y - 0.02 * y_range], color='black',
            linewidth=0.8 * SCALE * MARK)
    ax.text(0.5, bracket_y + 0.02 * y_range, p_str, ha='center',
            fontsize=style.body_pt() if is_star else style.tick_pt())

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"R\n(n={len(r_data)})", f"NR\n(n={len(nr_data)})"])
    ax.set_ylabel('IL-6/JAK/STAT3\nScore')
    ax.set_title('IL-6/JAK/STAT3\nCD4+ T cells')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(ax.get_ylim()[0], bracket_y + 0.15 * y_range)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    stem = 'il6_stat3_cd4t_boxplot'
    style.save_panel(fig, OUT_DIR / stem)
    print(f"  Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
