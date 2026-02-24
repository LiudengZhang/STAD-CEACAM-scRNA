#!/usr/bin/env python3
"""Panel O: IL-6/JAK/STAT3 Signaling Score — CD4+ T cells, Post-R vs Post-NR."""
import warnings, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD4_H5AD

warnings.filterwarnings('ignore')

OUT_DIR = Path(__file__).parent
SCALE = 4; CM = 1 / 2.54

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
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'svg.fonttype': 'none',
        'pdf.fonttype': 42, 'ps.fonttype': 42,
    })

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
    sc.tl.score_genes(adata, gene_list=avail, score_name='value',
                     ctrl_size=min(50, len(avail)), use_raw=True)

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

    fig, ax = plt.subplots(figsize=(3.5 * SCALE * CM, 5 * SCALE * CM))
    bp = ax.boxplot([r_data, nr_data], positions=[0, 1], widths=0.5,
                   patch_artist=True, showfliers=False,
                   medianprops=dict(color='black', linewidth=1.5 * SCALE),
                   whiskerprops=dict(linewidth=1.0 * SCALE),
                   capprops=dict(linewidth=1.0 * SCALE),
                   boxprops=dict(linewidth=1.0 * SCALE))
    bp['boxes'][0].set_facecolor(BOX_COLORS['R']); bp['boxes'][0].set_alpha(0.6)
    bp['boxes'][1].set_facecolor(BOX_COLORS['NR']); bp['boxes'][1].set_alpha(0.6)

    rng = np.random.default_rng(42)
    for k, (data, color) in enumerate(zip([r_data, nr_data], [BOX_COLORS['R'], BOX_COLORS['NR']])):
        jitter = rng.uniform(-0.08, 0.08, len(data))
        ax.scatter([k] * len(data) + jitter, data, c=color, s=20 * SCALE,
                  edgecolors='white', linewidths=0.3 * SCALE, alpha=0.85, zorder=3)

    _, pval = stats.mannwhitneyu(nr_data, r_data, alternative='greater')
    print(f"  P-value (one-sided NR>R): {pval:.4f}")

    p_str = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
    is_star = pval < 0.05
    y_max = max(np.max(r_data), np.max(nr_data))
    y_range = y_max - min(np.min(r_data), np.min(nr_data))
    bracket_y = y_max + 0.10 * y_range
    ax.plot([0, 0, 1, 1], [bracket_y - 0.02 * y_range, bracket_y,
            bracket_y, bracket_y - 0.02 * y_range], color='black', linewidth=0.8 * SCALE)
    ax.text(0.5, bracket_y + 0.02 * y_range, p_str, ha='center',
            fontsize=(7 if is_star else 5) * SCALE,
            fontweight='bold' if is_star else 'normal')

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"R\n(n={len(r_data)})", f"NR\n(n={len(nr_data)})"],
                       fontsize=5 * SCALE)
    ax.set_ylabel('IL-6/JAK/STAT3\nScore', fontsize=5 * SCALE)
    ax.set_title('IL-6/JAK/STAT3\nCD4+ T cells', fontsize=6 * SCALE, fontweight='normal')
    ax.tick_params(labelsize=5 * SCALE, width=1.0 * SCALE, length=4 * SCALE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(1.0 * SCALE)

    ax.set_ylim(ax.get_ylim()[0], bracket_y + 0.15 * y_range)

    plt.tight_layout()
    stem = 'il6_stat3_cd4t_boxplot'
    for ext in ['svg', 'png']:
        kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if ext == 'png':
            kw['dpi'] = 300
        fig.savefig(OUT_DIR / f'{stem}.{ext}', **kw)
    print(f"  Saved: {stem}")
    plt.close()


if __name__ == '__main__':
    main()
