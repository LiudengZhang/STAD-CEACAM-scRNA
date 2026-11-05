#!/usr/bin/env python3
"""Supplementary Figure S3 F: heat map of the top-10 consensus genes of the five
stomach metaprograms (rank-weighted votes over member NMF programs, scaled to
the metaprogram maximum).

Reads data/loadings.csv; the per-sample NMF gene lists are read only if that
table is missing or with --recompute. Run: python create_S3_F_metaprogram_heatmap.py"""

import sys
from collections import Counter
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import NMF_PER_SAMPLE, NMF_INTERMEDIATE        # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_F"
MP_COLORS = {'MP1': '#E64B35', 'MP2': '#4DBBD5', 'MP3': '#00A087',
             'MP4': '#3C5488', 'MP5': '#F39B7F'}
ALL_MPS = ['MP1', 'MP2', 'MP3', 'MP4', 'MP5']
TOP_N_GENES = 10

# Panel size in mm.
W, H = 171.0, 36.0
MARGIN = dict(left=12.0, right=14.0, top=1.5, bottom=11.5)
CBAR_W_MM = 2.0
# Metaprogram colour strip left of the matrix (gene-column units) and the
# minimum gap in mm between it and the row labels.
STRIP_X, STRIP_W = -0.32, 0.18
STRIP_GAP_MM = 0.8

def _program_genes(program_id):
    sample, k, p_num = program_id.rsplit('_', 2)
    k, p_num = int(k[1:]), int(p_num[1:])
    f = NMF_PER_SAMPLE / sample / f'{sample}_nmf_k{k}_genes.csv'
    if f.exists():
        df = pd.read_csv(f)
        return list(df[df.columns[p_num - 1]].dropna().tolist()[:50])
    return []

def compute_loadings():
    assignments = pd.read_csv(NMF_INTERMEDIATE / 'panel_A2_stomach_all_mp_assignments.csv')
    top, counts = {}, {}
    for mp in ALL_MPS:
        c = Counter()
        for pid in assignments[assignments['metaprogram_id'] == mp]['program_id']:
            for i, g in enumerate(_program_genes(pid)):
                c[g] += (50 - i)
        top[mp] = [g for g, _ in c.most_common(50)][:TOP_N_GENES]
        counts[mp] = c
    genes, gene_mp = [], {}
    for mp in ALL_MPS:
        for g in top[mp]:
            if g not in gene_mp:
                genes.append(g); gene_mp[g] = mp
    rows = []
    for mp in ALL_MPS:
        mx = max(counts[mp].values()) if counts[mp] else 1
        for j, g in enumerate(genes):
            rows.append({'mp': mp, 'gene': g, 'gene_order': j, 'gene_mp': gene_mp[g],
                         'loading': counts[mp].get(g, 0) / mx})
    return pd.DataFrame(rows)

def draw(tab):
    import seaborn as sns
    import matplotlib.pyplot as plt
    base.apply_style()
    genes = tab.drop_duplicates('gene').sort_values('gene_order')['gene'].tolist()
    gene_mp = dict(zip(tab['gene'], tab['gene_mp']))
    mat = tab.pivot(index='mp', columns='gene', values='loading').loc[ALL_MPS, genes]
    mat.index = [f'S-{mp}' for mp in ALL_MPS]
    fig, ax = style.subplots_mm(W, H)
    style.margins_mm(fig, **MARGIN)
    cax = fig.add_axes([(W - MARGIN['right'] + 1.5) / W, 1 - (MARGIN['top'] + 20.0) / H,
                        CBAR_W_MM / W, 18.0 / H])
    sns.heatmap(mat, cmap='YlOrRd', vmin=0, vmax=1, ax=ax, xticklabels=True,
                yticklabels=True, linewidths=0.3, linecolor='lightgray',
                cbar_ax=cax, cbar_kws={'label': 'Normalized Loading'})
    for coll in ax.collections:
        coll.set_rasterized(False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center', style='italic',
                       fontsize=style.tick_pt())
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=style.tick_pt())
    ax.tick_params(axis='x', length=1.5, width=style.RULE_PT)
    col_mm = (W - MARGIN['left'] - MARGIN['right']) / len(genes)
    strip_left_mm = -STRIP_X * col_mm                 # from the matrix edge
    ax.tick_params(axis='y', length=0,
                   pad=(strip_left_mm + STRIP_GAP_MM) * style.PT_PER_MM)
    ax.set_xlabel(''); ax.set_ylabel('')      # seaborn names the pivot's axes
    cax.tick_params(labelsize=style.tick_pt(), length=1.5, width=style.RULE_PT)
    cax.yaxis.label.set_size(style.tick_pt())
    cax.yaxis.label.set_rotation(90)
    cum = 0
    for mp in ALL_MPS:
        cum += sum(1 for g in genes if gene_mp[g] == mp)
        if cum < len(genes):
            ax.axvline(cum, color='black', linewidth=style.RULE_PT)
    strips = []
    for i, mp in enumerate(ALL_MPS):
        strips.append(ax.add_patch(plt.Rectangle(
            (STRIP_X, i), STRIP_W, 1, facecolor=MP_COLORS[mp],
            edgecolor='none', clip_on=False)))
    # Check that every row label ends at least STRIP_GAP_MM left of the strip.
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    strip_x0 = min(p_.get_window_extent(r).x0 for p_ in strips)
    for t in ax.get_yticklabels():
        gap_mm = (strip_x0 - t.get_window_extent(r).x1) / fig.dpi * 25.4
        if gap_mm < STRIP_GAP_MM - 0.05:
            raise RuntimeError(f"row label {t.get_text()!r} ends {gap_mm:.2f} mm "
                               f"from the colour strip; {STRIP_GAP_MM} mm wanted")
    for sp in ax.spines.values():
        sp.set_visible(False)
    return fig

if __name__ == "__main__":
    base.save(draw(cache.table(HERE, "loadings", compute_loadings)),
              FIG, PANEL, "panel_S3_F")
