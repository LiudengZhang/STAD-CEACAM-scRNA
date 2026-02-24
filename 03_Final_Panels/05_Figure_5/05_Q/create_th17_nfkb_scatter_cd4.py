#!/usr/bin/env python3
"""Panel R: NF-kB vs Th17 score in CD4+ T cells"""
import warnings, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD4_H5AD

warnings.filterwarnings('ignore')

H5AD = TCD4_H5AD
OUT_DIR = Path(__file__).parent
SCALE = 4; CM = 1 / 2.54

GROUP_COLORS = {
    'Pre-R': '#bde0fe', 'Pre-NR': '#a2d2ff',
    'Post-R': '#ffcfd2', 'Post-NR': '#f1c0e8',
    'Other': '#cccccc',
}
MIN_CELLS = 20

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
    'IL17A','IL17F','RORC','CCR6','IL23R','IL22','AHR','BATF','IRF4','STAT3',
    'CCL20','CXCR6','KLRB1','IL21','IL1R1','RORA','CTSH','PTPN13','TMEM176A','TMEM176B',
    'CAPG','LGMN','FKBP5','ICOS',
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
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'svg.fonttype': 'none',
        'pdf.fonttype': 42, 'ps.fonttype': 42,
    })

    adata = sc.read_h5ad(H5AD)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    adata.obs['group'] = adata.obs.apply(assign_group, axis=1).astype(str)

    gene_names = list(adata.raw.var_names) if adata.raw else list(adata.var_names)
    nfkb_avail = [g for g in HALLMARK_NFKB if g in gene_names]
    state_avail = [g for g in STATE_GENES if g in gene_names]
    print(f"NF-kB: {len(nfkb_avail)}/{len(HALLMARK_NFKB)}, State: {len(state_avail)}/{len(STATE_GENES)}")

    sc.tl.score_genes(adata, gene_list=nfkb_avail, score_name='nfkb',
                     ctrl_size=50, use_raw=True)
    sc.tl.score_genes(adata, gene_list=state_avail, score_name='state',
                     ctrl_size=min(50, len(state_avail)), use_raw=True)

    df = adata.obs[['sample', 'group', 'nfkb', 'state']].copy()
    df['sample'] = df['sample'].astype(str)
    df['group'] = df['group'].astype(str)
    sample_df = df.groupby(['sample', 'group'], observed=True).agg(
        nfkb=('nfkb', 'mean'), state=('state', 'mean'), n_cells=('nfkb', 'size'),
    ).reset_index()
    sample_df = sample_df[sample_df['n_cells'] >= MIN_CELLS]

    rho, pval = stats.spearmanr(sample_df['nfkb'], sample_df['state'])
    print(f"Spearman: rho={rho:.3f}, P={pval:.4f} (n={len(sample_df)})")

    fig, ax = plt.subplots(figsize=(5 * SCALE * CM, 5 * SCALE * CM))

    for g in ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']:
        mask = sample_df['group'] == g
        if mask.sum() == 0:
            continue
        ax.scatter(sample_df.loc[mask, 'nfkb'], sample_df.loc[mask, 'state'],
                  c=GROUP_COLORS[g], edgecolors='white', linewidths=0.3 * SCALE,
                  s=30 * SCALE, alpha=0.85, zorder=3 if g != 'Other' else 2,
                  label=f"{g} (n={mask.sum()})")

    x = sample_df['nfkb'].values
    y = sample_df['state'].values
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() > 2:
        z = np.polyfit(x[valid], y[valid], 1)
        x_line = np.linspace(x[valid].min(), x[valid].max(), 100)
        ax.plot(x_line, np.polyval(z, x_line), 'k--', linewidth=0.8 * SCALE, alpha=0.6, zorder=2)

    p_str = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'

    ax.set_xlabel('NF-\u03baB Score', fontsize=6 * SCALE)
    ax.set_ylabel('Th17 Score', fontsize=6 * SCALE)
    ax.set_title(f'Th17 (CD4+)\n\u03c1 = {rho:.2f}, {p_str}', fontsize=6.5 * SCALE, fontweight='normal', linespacing=1.4)
    ax.tick_params(labelsize=5 * SCALE, width=0.5 * SCALE, length=3 * SCALE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(0.5 * SCALE)
    ax.set_box_aspect(1)

    plt.tight_layout()
    stem = 'th17_nfkb_scatter_cd4'
    for ext in ['svg', 'png']:
        kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if ext == 'png':
            kw['dpi'] = 300
        fig.savefig(OUT_DIR / f'{stem}.{ext}', **kw)
    print(f"Saved: {stem}")
    plt.close()


if __name__ == '__main__':
    main()
