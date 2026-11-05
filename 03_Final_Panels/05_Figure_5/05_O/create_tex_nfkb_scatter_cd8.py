#!/usr/bin/env python3
"""Figure 5K: Hallmark TNF-alpha/NF-kB score against the Tex exhaustion score
in CD8+ T cells, one point per stomach sample (samples with at least 20
cells), with Spearman rho and P and a least-squares line.

draw() is shared with panel M (05_Q/create_th17_nfkb_scatter_cd4.py). rho and
P are also written to a table by cnsfig.corr_stats.

Input:  TCD8_H5AD (00_Config/paths.py)
Output: tex_nfkb_scatter_cd8.{svg,pdf,png} in this directory
Run:    python create_tex_nfkb_scatter_cd8.py
"""
import warnings, numpy as np, scanpy as sc
import matplotlib; matplotlib.use('Agg')
from scipy import stats
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD8_H5AD
from cnsfig import Panel, style  # noqa: E402
from cnsfig import layout as cnslayout, corr_stats, rich_xlabel, rich_ylabel
# slots gives panel J's rectangle, used to align this panel's title with J's.
import slots  # noqa: E402

warnings.filterwarnings('ignore')

H5AD = TCD8_H5AD

# Response/phase colours shared by all figures.
from shared.figure_config import GROUP_COLORS as _SCHEME  # noqa: E402
GROUP_COLORS = {**_SCHEME, 'Other': '#cccccc'}
MIN_CELLS = 20

PANEL = Panel(5, "K", __file__, small_pt=5.0)
#: Left margin for the y label and ticks, and the plot box width (mm).
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

STATE_GENES = [
    'PDCD1','CTLA4','HAVCR2','LAG3','TIGIT','CD244','CD160','TOX','TOX2','EOMES',
    'PRDM1','NR4A1','NR4A2','NR4A3','ENTPD1','CD38','CXCL13','LAYN','TNFRSF9',
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


def draw(panel, h5ad, state_genes, title, ylabel, stem):
    """One of the two scatter panels, at the slot `panel` names."""
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)

    adata = sc.read_h5ad(h5ad)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata.obs['group'] = adata.obs.apply(assign_group, axis=1).astype(str)

    gene_names = list(adata.raw.var_names) if adata.raw else list(adata.var_names)
    nfkb_avail = [g for g in HALLMARK_NFKB if g in gene_names]
    state_avail = [g for g in state_genes if g in gene_names]

    # Use .raw when present, else .X.
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

    # Square plot box level with the frames of panels J and L (page mm).
    fig = style.figure_mm(panel.w_mm, panel.h_mm)
    ax, = cnslayout.scatter_pair_mm(fig, left_mm=LEFT_MM, box_mm=BOX_W_MM,
                                    top_mm=139.8 - 130.0, n=1)
    _pos = ax.get_position()
    ax.set_position([_pos.x0, (166.5 - 156.4) / panel.h_mm, _pos.width,
                     (156.4 - 139.8) / panel.h_mm])

    for g in ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']:
        mask = sample_df['group'] == g
        if mask.sum() == 0:
            continue
        ax.scatter(sample_df.loc[mask, 'nfkb'], sample_df.loc[mask, 'state'],
                  c=GROUP_COLORS[g], edgecolors='white',
                  linewidths=style.EDGE_PT,
                  s=(0.48 * style.PT_PER_MM) ** 2, alpha=0.85,  # 0.48 mm dots
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

    corr_stats.write(panel.dir, [('NF-κB score', rho, pval, len(sample_df))])
    # Title on the same line as the titles of panels J and L.
    _row_title_top_mm = slots.rect_mm(5, 'J', sub=1)[1] + 0.4 + 0.27
    _own_top_mm = slots.rect_mm(5, panel.letter)[1]
    _pos = ax.get_position()
    fig.text((_pos.x0 + _pos.x1) / 2,
             1.0 - (_row_title_top_mm - _own_top_mm) / panel.h_mm,
             title, ha='center', va='top', fontsize=style.body_pt(),
             linespacing=1.15)
    cnslayout.corr_annotate(ax, rho, p=pval)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    rich_xlabel(ax, 'NF-κB Score')
    rich_ylabel(ax, ylabel)
    from cnsfig.layout import recover_x
    recover_x(fig, ax)
    panel.finish(fig, stem)
    print(f"  {title}: rho = {rho:.2f}, {corr_stats.p_string(pval)} "
          f"(n={len(sample_df)})")


def main():
    draw(PANEL, H5AD, STATE_GENES, 'Tex (CD8+)', 'Tex (CD8+) Score',
         'tex_nfkb_scatter_cd8')


if __name__ == '__main__':
    main()
