#!/usr/bin/env python3
"""
Figure 5 panel M - NF-kB score against the Th17 score in CD4+ T cells, one
point per sample.

The 200-gene Hallmark NF-kB list, the 24-gene Th17 list, the twenty-cell sample
filter, `sc.tl.score_genes`, the Spearman test and the P-value rounding are
unchanged.

  printed panel  Figure 5 M       (PROVENANCE.csv; NOT inferred from "05_Q")

The correlation and its P value take a line each: set on one line they are
wider than the plotting box, and a title wider than the box it is centred on
cannot be fitted at all - narrowing the box by a millimetre moves the title's
edge by half of one. The break is declared in 00_Config/shared/labels.py; the
statistics do not move.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the tick labels - at 5 * SCALE, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The marker area, the marker edge width and the dashed trend line are scaled by those. Tick widths and lengths and spine widths are
    not: those are style, and cnsplots sets them.

Judgement call: the title's explicit `fontweight='normal'` is dropped so
cnsplots' bold axis title applies. Its `linespacing=1.4` is kept - that is
layout of a multi-line string, not type.
"""
import warnings, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD4_H5AD
import panel_style_cns as style
import slots
from cnsfig import layout as cnslayout, corr_stats, rich_xlabel, rich_ylabel

warnings.filterwarnings('ignore')

H5AD = TCD4_H5AD
OUT_DIR = Path(__file__).parent
SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier smallest body type, before * SCALE

GROUP_COLORS = {
    'Pre-R': '#bde0fe', 'Pre-NR': '#a2d2ff',
    'Post-R': '#ffcfd2', 'Post-NR': '#f1c0e8',
    'Other': '#cccccc',
}
MIN_CELLS = 20

PANEL_LETTER = "M"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(5, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(5, PANEL_LETTER)
#: The y label, its ticks and the box; the box is 18.9 mm tall (the row's
#: frame line) and as wide as the slot allows.
LEFT_MM = 7.8
BOX_W_MM = 13.8

HALLMARK_NFKB = [
    'ABCA1','ACKR3','AREG','ATF3','ATP2B1','B4GALT1','B4GALT5','BCL2A1','BCL3','BCL6',
    'BHLHE40','BIRC2','BIRC3','BMP2','BTG1','BTG2','BTG3','CCL2','CCL20','CCL4',
    'CCL5','CCND1','CCNL1','CCR7','CCRL2','CD44','CD69','CD80','CD83','CDKN1A',
    'CEBPB','CEBPD','CFLAR','CLCF1','CSF1','CSF2','CXCL1','CXCL10','CXCL11','CXCL2',
    'CXCL3','CXCL6','CXCL8','DENND5A','DUSP1','DUSP2','DUSP4','DUSP5','EDN1','EFNA1',
    'EGR1','EGR2','EGR3','EHD1','EIF1','ETS2','F2RL1','F3','FJX1','FOS',
    'FOSB','FOSL1','FOSL2','FUT4','G0S2','GADD45A','GADD45B','GCH1','GEM','GFPT2',
    'GPR183','HBEGF','HES1','ICAM1','ICOSLG','ID2','IER2','IER3','IER5','IFIH1',
    'IFNGR2','IL12B','IL15RA','IL18','IL1A','IL1B','IL23A','IL6','IL6ST','IL7R',
    'INHBA','IRF1','IRS2','JAG1','JUN','JUNB','KDM6B','KLF10','KLF2','KLF4',
    'KLF6','KLF9','KYNU','LAMB3','LIF','LITAF','MAFF','MAP2K3','MAP3K8','MARCKS',
    'MCL1','MSC','MXD1','MYC','NAMPT','NFAT5','NFE2L2','NFIL3','NFKB1','NFKB2',
    'NFKBIA','NFKBIE','NIN','NR4A1','NR4A2','NR4A3','OLR1','PANX1','PDE4B','PDLIM5',
    'PER1','PFKFB3','PHLDA1','PHLDA2','PIK3R1','PLAU','PLAUR','PLEK','PLK2','PLPP3',
    'PMEPA1','PNRC1','PPP1R15A','PTGER4','PTGS2','PTX3','RCAN1','REL','RELA','RELB',
    'RHOB','RIPK2','RNF19B','SAT1','SDC4','SERPINB2','SERPINB8','SERPINE1','SGK1','SIK1',
    'SLC16A6','SLC2A3','SLC2A6','SMAD3','SNN','SOCS3','SOD2','SPHK1','SQSTM1','TANK',
    'TGIF1','TIPARP','TLR2','TNC','TNF','TNFAIP2','TNFAIP3','TNFAIP6','TNFAIP8','TNFRSF9',
    'TNFSF9','TNIP1','TNIP2','TRAF1','TRIB1','TRIP10','TSC22D1','TUBB2A','VEGFA','YRDC',
    'ZC3H12A','ZFP36',
]

# The Th17 signature behind the published panel: rho = 0.397, P = 0.0243, which
# is the "rho = 0.40, P = 0.02" printed in Figure 5M and quoted in the Results.
#
# A 28-gene variant sat here instead - IL26, IL6R, TGFBR2, CSF2, IFNG, TNF, IL2,
# CTLA4, CD44, IL4I1, LGALS3 and CXCR3 added, CXCR6, CTSH, PTPN13, TMEM176A,
# TMEM176B, CAPG, LGMN and FKBP5 dropped. It was written on 23 Feb 2026, one day
# after the version that made the figure, and no figure was ever regenerated
# from it. It scores rho = 0.655, P = 5e-5: a different published claim.
#
# Verified by running both lists against TCD4.h5ad in one pass, unchanged since
# Oct 2025 - 24 genes give 0.397 and a Th17 range of -0.014..0.211, matching the
# printed y axis of 0.00..0.20; 28 genes give 0.655 over 0.047..0.248.
STATE_GENES = [
    'IL17A','IL17F','RORC','CCR6','IL23R','IL22','AHR','BATF','IRF4','STAT3',
    'CCL20','CXCR6','KLRB1','IL21','IL1R1','RORA','CTSH','PTPN13','TMEM176A',
    'TMEM176B','CAPG','LGMN','FKBP5','ICOS',
]


def assign_group(row):
    treat = row.get('Treatment phase', '')
    if treat == 'Pre':
        g = row.get('stomach_pre_grouping', '')
        if g == 'Responsed': return 'Pre-R'
        if g == 'No-response': return 'Pre-NR'
    elif treat == 'Post':
        g = row.get('stomach_post_grouping', '')
        if g == 'Responsed': return 'Post-R'
        if g == 'No-response': return 'Post-NR'
    return 'Other'


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    adata = sc.read_h5ad(H5AD)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata.obs['group'] = adata.obs.apply(assign_group, axis=1).astype(str)

    gene_names = list(adata.raw.var_names) if adata.raw else list(adata.var_names)
    nfkb_avail = [g for g in HALLMARK_NFKB if g in gene_names]
    state_avail = [g for g in STATE_GENES if g in gene_names]
    print(f"NF-kB: {len(nfkb_avail)}/{len(HALLMARK_NFKB)}, State: {len(state_avail)}/{len(STATE_GENES)}")

    # score_genes with use_raw=True reaches straight into adata.raw.var_names,
    # so it has to be told when the file carries no .raw - the clean deposit
    # holds the same log1p matrix in .X.
    use_raw = adata.raw is not None
    sc.tl.score_genes(adata, gene_list=nfkb_avail, score_name='nfkb',
                     ctrl_size=50, use_raw=use_raw)
    sc.tl.score_genes(adata, gene_list=state_avail, score_name='state',
                     ctrl_size=min(50, len(state_avail)), use_raw=use_raw)

    df = adata.obs[['sample', 'group', 'nfkb', 'state']].copy()
    df['sample'] = df['sample'].astype(str)
    df['group'] = df['group'].astype(str)
    sample_df = df.groupby(['sample', 'group'], observed=True).agg(
        nfkb=('nfkb', 'mean'), state=('state', 'mean'), n_cells=('nfkb', 'size'),
    ).reset_index()
    sample_df = sample_df[sample_df['n_cells'] >= MIN_CELLS]

    rho, pval = stats.spearmanr(sample_df['nfkb'], sample_df['state'])
    print(f"Spearman: rho={rho:.3f}, P={pval:.4f} (n={len(sample_df)})")

    # THE BOX AT MILLIMETRES, RHO INSIDE, P IN THE LEGEND  (2026-09-14, evening)
    #   The same rule as Figure 3's scatter pairs: the title is the state
    #   name alone, rho is drawn inside the box, the P value is written by
    #   cnsfig.corr_stats and the legend is generated from it. The box is
    #   placed by cnsfig.layout.scatter_pair_mm with its top at 139.8 mm and
    #   its bottom at 158.7 mm on the page, the frame line of J and L.
    fig = style.figure_mm(PANEL_W_MM, PANEL_H_MM)
    ax, = cnslayout.scatter_pair_mm(fig, left_mm=LEFT_MM, box_mm=BOX_W_MM,
                                    top_mm=139.8 - 130.0, n=1)
    _pos = ax.get_position()
    ax.set_position([_pos.x0, (166.5 - 158.7) / PANEL_H_MM, _pos.width,
                     (158.7 - 139.8) / PANEL_H_MM])

    for g in ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']:
        mask = sample_df['group'] == g
        if mask.sum() == 0:
            continue
        ax.scatter(sample_df.loc[mask, 'nfkb'], sample_df.loc[mask, 'state'],
                  c=GROUP_COLORS[g], edgecolors='white',
                  linewidths=style.EDGE_PT,
                  s=(0.48 * style.PT_PER_MM) ** 2, alpha=0.85,  # published 0.47-0.49 mm
                  zorder=3 if g != 'Other' else 2,
                  label=f"{g} (n={mask.sum()})")

    x = sample_df['nfkb'].values
    y = sample_df['state'].values
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() > 2:
        z = np.polyfit(x[valid], y[valid], 1)
        x_line = np.linspace(x[valid].min(), x[valid].max(), 100)
        ax.plot(x_line, np.polyval(z, x_line), 'k--', linewidth=style.RULE_PT,
                alpha=0.6, zorder=2)

    p_str = corr_stats.p_string(pval)
    print(f"    Th17 (CD4+): rho = {rho:.2f}, {p_str}")
    corr_stats.write(OUT_DIR, [('NF-\u03baB score', rho, pval, len(sample_df))])
    # THE TITLE ON THE ROW'S TITLE LINE  (2026-09-16, the author's fifth
    #   reading: "align the 5K and 5M subtitles with the rest of their row").
    #   As an axes title it hung 1.4 mm over the plotting box, 2.87 mm below
    #   J's and L's titles. J and L set theirs through finish_two_group: the
    #   text's top at 0.4 mm + top_extra_mm below their canvas top, which
    #   for J's first box (slot top 133.05, top_extra 0.27) is 133.72 mm on
    #   the page - the same page line for all four J boxes and for L. The
    #   same top, the same 7 pt and linespacing, drawn on this canvas
    #   (sweep_pages.check_title_rows holds the row to 0.5 mm).
    _row_title_top_mm = slots.rect_mm(5, 'J', sub=1)[1] + 0.4 + 0.27
    _own_top_mm = slots.rect_mm(5, PANEL_LETTER)[1]
    _pos = ax.get_position()
    fig.text((_pos.x0 + _pos.x1) / 2,
             1.0 - (_row_title_top_mm - _own_top_mm) / PANEL_H_MM,
             'Th17 (CD4+)', ha='center', va='top', fontsize=style.body_pt(),
             linespacing=1.15)
    cnslayout.corr_annotate(ax, rho)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    rich_xlabel(ax, 'NF-\u03baB Score')
    rich_ylabel(ax, 'Th17 Score')
    from cnsfig.layout import recover_x
    recover_x(fig, ax)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    stem = 'th17_nfkb_scatter_cd4'
    style.save_panel(fig, OUT_DIR / stem)
    print(f"Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
