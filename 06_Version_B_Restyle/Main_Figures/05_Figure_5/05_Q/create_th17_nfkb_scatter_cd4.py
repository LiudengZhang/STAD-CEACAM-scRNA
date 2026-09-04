#!/usr/bin/env python3
"""
Figure 5 panel M, RESTYLED (Version B) - NF-kB score against the Th17 score in
CD4+ T cells, one point per sample.

Version A is
`03_Revised_Panels/Main_Figures/05_Figure_5/05_Q/create_th17_nfkb_scatter_cd4.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK
(areas by AREA), margins are millimetres, and the save is `style.save_panel`.

The 200-gene Hallmark NF-kB list and the 24-gene Th17 list are Version A's,
character for character - including the comment recording why the 24-gene list
and not the 28-gene variant is the one that made the published panel. Not one
gene, filter or statistic is touched here.

  printed panel  Figure 5 M        (PROVENANCE.csv; NOT inferred from "05_Q")
  printed rect   23.6 x 27.5 mm    (panel_rects.csv)
  Version B box  50.0 x 55.0 mm

MARK
    Version A drew 5 x 5 cm at SCALE = 4 and set its smallest body type - the
    tick labels - at `5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2                                     = 0.1225

Judgement call: as 05_O - the title's explicit `fontweight='normal'` is
dropped so cnsplots' bold axis title applies; `linespacing=1.4` is kept.
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

H5AD = TCD4_H5AD
OUT_DIR = Path(__file__).parent
SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 5.0                  # Version A's smallest body type, before * SCALE

GROUP_COLORS = {
    'Pre-R': '#bde0fe', 'Pre-NR': '#a2d2ff',
    'Post-R': '#ffcfd2', 'Post-NR': '#f1c0e8',
    'Other': '#cccccc',
}
MIN_CELLS = 20

PRINTED_MM = (23.6, 27.5)
PANEL_W_MM, PANEL_H_MM = 50.0, 55.0
MARGIN = dict(left=13.0, right=2.0, top=9.5, bottom=10.0)

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
    family = style.apply()
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

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    for g in ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']:
        mask = sample_df['group'] == g
        if mask.sum() == 0:
            continue
        ax.scatter(sample_df.loc[mask, 'nfkb'], sample_df.loc[mask, 'state'],
                  c=GROUP_COLORS[g], edgecolors='white',
                  linewidths=0.3 * SCALE * MARK,
                  s=30 * SCALE * AREA, alpha=0.85,
                  zorder=3 if g != 'Other' else 2,
                  label=f"{g} (n={mask.sum()})")

    x = sample_df['nfkb'].values
    y = sample_df['state'].values
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() > 2:
        z = np.polyfit(x[valid], y[valid], 1)
        x_line = np.linspace(x[valid].min(), x[valid].max(), 100)
        ax.plot(x_line, np.polyval(z, x_line), 'k--', linewidth=0.8 * SCALE * MARK,
                alpha=0.6, zorder=2)

    import math
    _e = math.floor(math.log10(pval)); _c = int(pval / 10**_e)
    if _e >= -3:
        p_str = f'P = {_c * 10**_e:.{-_e}f}'
    else:
        p_str = f'P = {_c}' + r'$\times 10^{' + str(_e) + r'}$'

    ax.set_xlabel('NF-\u03baB Score')
    ax.set_ylabel('Th17 Score')
    ax.set_title(f'Th17 (CD4+)\n\u03c1 = {rho:.2f}, {p_str}', linespacing=1.4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_box_aspect(1)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    stem = 'th17_nfkb_scatter_cd4'
    style.save_panel(fig, OUT_DIR / stem)
    print(f"Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
