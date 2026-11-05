#!/usr/bin/env python3
"""Per-sample QC box plots: genes per cell, total counts and mitochondrial %
stacked on one sample axis, drawn from a cached table of box statistics
(matplotlib boxplot_stats, whis = 1.5, no fliers). Shared by Supplementary
panels S1 B (stomach samples) and S1 C (non-stomach samples).
"""

import re

import pandas as pd

import panel_style_cns as style

FIG = "S1_QC_Annotation"
METRICS = [('n_genes_by_counts', 'Genes per Cell', 'Number of genes'),
           ('total_counts', 'Total Counts', 'Total UMI counts'),
           ('pct_counts_mt', 'Mitochondrial %', '% mitochondrial')]
METRIC_COLS = [m for m, _, _ in METRICS]
STAT_COLS = ('med', 'q1', 'q3', 'whislo', 'whishi')
HSPACE = 0.6

def box_stats(df, info, group_order):
    """Per-sample boxes, as numbers. `info` carries one row per sample with a
    'group' column; the samples are ordered by group and then by specimen.
    """
    from matplotlib.cbook import boxplot_stats
    info = info.copy()
    info['sort_key'] = info['group'].map({g: i for i, g in enumerate(group_order)})
    info = info.sort_values(['sort_key', 'sample'])
    rows = []
    for order, (_, r) in enumerate(info.iterrows()):
        label = str(r['Sample ID'])
        if not re.fullmatch(r'P\d+-\w+', label):
            raise SystemExit(f'not a study sample ID: {r["sample"]} -> {label}')
        sub = df[df['sample'] == r['sample']]
        for metric in METRIC_COLS:
            s = boxplot_stats(sub[metric].to_numpy(dtype=float))[0]
            rows.append({'order': order, 'label': label, 'group': r['group'],
                         'metric': metric, **{k: float(s[k]) for k in STAT_COLS}})
    return pd.DataFrame(rows)

def draw(stats, *, w_mm, h_mm, margin, colors, group_order, legend_all):
    """The three stacked strips. `legend_all` keys every group of
    `group_order`; otherwise only the groups the data carry (S1 C)."""
    import matplotlib.patches as mpatches
    import _driver_base as base
    base.apply_style()
    fig, axes = style.subplots_mm(w_mm, h_mm, 3, 1, sharex=True)
    style.margins_mm(fig, **margin)
    fig.subplots_adjust(hspace=HSPACE)
    samples = stats.drop_duplicates('order').sort_values('order')
    n = len(samples)
    labels = samples['label'].tolist()
    groups = samples['group'].tolist()
    boundaries = [i - 0.5 for i in range(1, n) if groups[i] != groups[i - 1]]
    for ax, (metric, title, _) in zip(axes, METRICS):
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
            patch.set_facecolor(colors[g])
            patch.set_alpha(1.0)   # the group colour as defined
        for b in boundaries:
            ax.axvline(b, color='gray', linestyle=':', linewidth=style.RULE_PT, alpha=0.7)
        if metric == 'pct_counts_mt':
            ax.axhline(25, color='#666666', linestyle='--', linewidth=style.RULE_PT, zorder=10)
        ax.set_title(title, pad=1.5)
        ax.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=style.EDGE_PT)
        ax.set_axisbelow(True)
        base.despine(ax)
        ax.set_xlim(-0.5, n - 0.5)
    axes[-1].set_xticks(list(range(n)))
    axes[-1].set_xticklabels(labels, rotation=90, ha='center')
    axes[-1].set_xlabel('Sample')
    keyed = group_order if legend_all else [g for g in group_order if g in set(groups)]
    handles = [mpatches.Patch(facecolor=colors[g], edgecolor='black', alpha=1.0,
                              linewidth=style.EDGE_PT, label=g) for g in keyed]
    fig.legend(handles=handles, loc='lower left', ncol=5, frameon=False,
               bbox_to_anchor=(margin['left'] / w_mm, 0.0), borderaxespad=0.0,
               columnspacing=1.0, handlelength=1.2, handletextpad=0.4)
    return fig

def panel(here, panel_dir, compute, *, w_mm, h_mm, margin, colors,
          group_order, legend_all):
    """Draw and save one S1 QC panel."""
    import _driver_base as base
    from cnsfig import cache
    stats = cache.table(here, "box_stats", compute)
    fig = draw(stats, w_mm=w_mm, h_mm=h_mm, margin=margin, colors=colors,
               group_order=group_order, legend_all=legend_all)
    base.save(fig, FIG, panel_dir, f"panel_{panel_dir}")
    return 0
