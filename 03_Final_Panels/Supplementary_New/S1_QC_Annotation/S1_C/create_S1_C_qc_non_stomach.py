#!/usr/bin/env python3
"""
S1 panel C - per-sample QC box plots of the 38 non-stomach samples (genes
per cell, total counts, mitochondrial %), coloured by tissue site, drawn at
the size they print at.

  printed panel  Supplementary Figure S1 C   (PROVENANCE.csv)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/S1_QC_Annotation/
                 S1_C/create_S1_C_qc_non_stomach.py  (three boxes side by
                 side on a 24 x 5 cm canvas at 4x; its 38 sample labels
                 printed at 2 pt)

THE THREE METRICS ARE STACKED, one above the other, sharing the sample axis,
as in panel B above it: 38 rotated 6 pt sample labels need 95 mm. The
sample order, the colours, the tissue separators, the 25 % mitochondrial
line, the legend and every string are the predecessor's.

DRAWING READS A TABLE (cnsfig.cache): data/box_stats.csv holds, per sample
and metric, the five numbers a box is made of (matplotlib's own
boxplot_stats with whis = 1.5, fliers not shown) plus the sample's study ID
and group; the h5ad is opened only when the table is absent or with
--recompute. ax.bxp() draws those numbers, which is what ax.boxplot() did
after computing them.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import FULL_DATASET_H5AD                       # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S1_QC_Annotation", "S1_C"

# The printed box, millimetres: a row to itself (38 samples at 4.2 mm).
# 46 mm tall since 2026-09-16 (S1 back on one page): the strips gave back
# 10 mm between them, the longer rotated labels ("P28-LN1_2") kept theirs.
W, H = 171.0, 46.0

COLORS = {
    'Liver': '#e6ab02',
    'Peripheral blood': '#66a61e',
    'Lymph node': '#7570b3',
    'Ovary': '#e7298a',
    'Metastatic lymph node': '#a6761d',
}
GROUP_ORDER = ['Liver', 'Peripheral blood', 'Lymph node', 'Ovary',
               'Metastatic lymph node']
METRICS = [('n_genes_by_counts', 'Genes per Cell', 'Number of genes'),
           ('total_counts', 'Total Counts', 'Total UMI counts'),
           ('pct_counts_mt', 'Mitochondrial %', '% mitochondrial')]
STAT_COLS = ('med', 'q1', 'q3', 'whislo', 'whishi')
#: Paper round the three rows: tick labels at the left, the rotated sample
#: labels, "Sample" and the one-row legend below, titles above.
MARGIN = dict(left=8.0, right=1.0, top=3.5, bottom=18.0)
HSPACE = 0.6


def compute_stats():
    """The computing half: the predecessor's per-sample boxes, as numbers."""
    import scanpy as sc
    from matplotlib.cbook import boxplot_stats
    adata = sc.read_h5ad(FULL_DATASET_H5AD, backed='r')
    mask = adata.obs['Sample site'] != 'Stomach'
    df = adata.obs.loc[mask, ['sample', 'Sample ID', 'Sample site']
                       + [m for m, _, _ in METRICS]]
    info = df[['sample', 'Sample ID', 'Sample site']].drop_duplicates('sample').copy()
    info['group'] = info['Sample site'].astype(str)
    info['sort_key'] = info['group'].map({g: i for i, g in enumerate(GROUP_ORDER)})
    info = info.sort_values(['sort_key', 'sample'])
    rows = []
    for order, (_, r) in enumerate(info.iterrows()):
        label = str(r['Sample ID'])
        if not re.fullmatch(r'P\d+-\w+', label):
            raise SystemExit(f'not a study sample ID: {r["sample"]} -> {label}')
        sub = df[df['sample'] == r['sample']]
        for metric, _, _ in METRICS:
            s = boxplot_stats(sub[metric].to_numpy(dtype=float))[0]
            rows.append({'order': order, 'label': label, 'group': r['group'],
                         'metric': metric, **{k: float(s[k]) for k in STAT_COLS}})
    return pd.DataFrame(rows)


def draw(stats):
    import matplotlib.patches as mpatches
    base.apply_style()
    fig, axes = style.subplots_mm(W, H, 3, 1, sharex=True)
    style.margins_mm(fig, **MARGIN)
    fig.subplots_adjust(hspace=HSPACE)
    samples = stats.drop_duplicates('order').sort_values('order')
    n = len(samples)
    labels = samples['label'].tolist()
    groups = samples['group'].tolist()
    boundaries = [i - 0.5 for i in range(1, n) if groups[i] != groups[i - 1]]
    for ax, (metric, title, ylabel) in zip(axes, METRICS):
        sub = stats[stats['metric'] == metric].sort_values('order')
        bxp = [dict(med=r.med, q1=r.q1, q3=r.q3, whislo=r.whislo, whishi=r.whishi,
                    fliers=[]) for r in sub.itertuples()]
        art = ax.bxp(bxp, positions=list(range(n)), widths=0.7, patch_artist=True,
                     showfliers=False,
                     medianprops=dict(color='black', linewidth=style.RULE_PT),
                     whiskerprops=dict(color='black', linewidth=style.EDGE_PT),
                     capprops=dict(color='black', linewidth=style.EDGE_PT),
                     boxprops=dict(linewidth=style.EDGE_PT, edgecolor='black'))
        for patch, g in zip(art['boxes'], groups):
            patch.set_facecolor(COLORS[g])
            patch.set_alpha(0.8)
        for b in boundaries:
            ax.axvline(b, color='gray', linestyle=':', linewidth=style.RULE_PT, alpha=0.7)
        if metric == 'pct_counts_mt':
            ax.axhline(25, color='#666666', linestyle='--', linewidth=style.RULE_PT, zorder=10)
        ax.set_title(title, pad=1.5)
        # No y label: the title names the metric and a 7 pt rotated label
        # is longer than the 12 mm row (labels.py REMOVALS_SUPPLEMENTARY).
        ax.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=style.EDGE_PT)
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(-0.5, n - 0.5)
    axes[-1].set_xticks(list(range(n)))
    axes[-1].set_xticklabels(labels, rotation=90, ha='center')
    axes[-1].set_xlabel('Sample')
    present = set(groups)
    handles = [mpatches.Patch(facecolor=COLORS[g], edgecolor='black', alpha=0.8,
                              linewidth=style.EDGE_PT, label=g)
               for g in GROUP_ORDER if g in present]
    fig.legend(handles=handles, loc='lower left', ncol=5, frameon=False,
               bbox_to_anchor=(MARGIN['left'] / W, 0.0), borderaxespad=0.0,
               columnspacing=1.0, handlelength=1.2, handletextpad=0.4)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"ink outside the {W} x {H} mm canvas: {over}")
    return fig


def main():
    stats = cache.table(HERE, "box_stats", compute_stats)
    base.save(draw(stats), FIG, PANEL, "panel_S1_C")
    return 0


if __name__ == "__main__":
    sys.exit(main())
