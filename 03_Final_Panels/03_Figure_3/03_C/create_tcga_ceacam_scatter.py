#!/usr/bin/env python3
"""
Figure 3B: CEACAM5 and CEACAM6 against CD274 in TCGA-STAD, using
BayesPrism-deconvolved epithelial expression (log2(x + 1)), one point per
tumour, filled (dead) or open (alive) by vital status. Tumours with zero
expression of any of the three genes are excluded, and each variable is
adjusted for ABSOLUTE purity (linear regression, mean added back). Spearman
rho and P are printed in each box and written to correlation_stats.csv.

draw(), purity_adjust() and vital_status() are shared with Figure 3E
(03_D/create_tcga_ceacam_scatter.py, a different file with the same name).

Reads data/tcga_samples.csv beside this script; with --recompute, computes it
from TCGA_BAYESPRISM_EPI, the sample list of TCGA_BULK_TUMOR_ONLY,
TCGA_ABSOLUTE_PURITY and TCGA_CLINICAL (vital_status).
Writes tcga_ceacam_cd274.[svg|pdf|png] beside this script.

Run: python create_tcga_ceacam_scatter.py [--recompute]
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy import stats
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import (TCGA_BAYESPRISM_EPI, TCGA_BULK_TUMOR_ONLY, TCGA_CLINICAL,  # noqa: E402
                   TCGA_ABSOLUTE_PURITY)
from cnsfig import layout as cnslayout, corr_stats, rich_xlabel, rich_ylabel  # noqa: E402
from cnsfig import cache, group_key, Panel, style  # noqa: E402

PANEL = Panel(3, "B", __file__, small_pt=4.5)
SCALE = Panel.SCALE
# Left margin, shared with Figure 3A above so the boxes align.
LEFT_MM = 12.3
# Marker area (about 0.42 mm across).
MARKER_AREA = 12.8

TCGA_EPI_EXPR = TCGA_BAYESPRISM_EPI

# Survival dot colors
COLOR_DEAD = '#c0392b'      # dark red, filled
COLOR_ALIVE = '#888888'     # gray edge, unfilled
COLOR_MISSING = '#cccccc'   # light gray edge, unfilled


def absolute_purity():
    """TCGA PanCanAtlas ABSOLUTE purity, one value per patient (mean over aliquots)."""
    ab = pd.read_csv(TCGA_ABSOLUTE_PURITY, sep='\t')
    ab['patient'] = ab['array'].astype(str).str[:12]
    return ab.dropna(subset=['purity']).groupby('patient')['purity'].mean()


def purity_adjust(df, cols):
    """Keep tumours with an ABSOLUTE call; regress each column on purity (OLS)
    and add its mean back, so values stay on their own scale."""
    pur = absolute_purity()
    df = df[df.index.isin(pur.index)].copy()
    X = np.column_stack([np.ones(len(df)), pur.loc[df.index].values])
    for c in cols:
        y = df[c].values.astype(float)
        df[c] = y - X @ np.linalg.lstsq(X, y, rcond=None)[0] + y.mean()
    return df


def vital_status(index):
    clin = pd.read_csv(TCGA_CLINICAL, sep='\t', usecols=['submitter_id', 'vital_status'])
    clin = clin.drop_duplicates(subset='submitter_id').set_index('submitter_id')
    return clin['vital_status'].reindex(index).fillna('Missing').values


def load_data_3b():
    """Figure 3B: log2(CPM + 1) of deconvolved epithelial CEACAM5, CEACAM6 and
    CD274; tumours with all three above zero; purity-adjusted."""
    epi = pd.read_csv(TCGA_EPI_EXPR, sep='\t', index_col=0)
    tum = pd.read_csv(TCGA_BULK_TUMOR_ONLY, sep='\t', index_col=0, nrows=1).columns
    ids = sorted(set(epi.index) & set(tum))
    df = np.log2(epi.loc[ids, ['CEACAM5', 'CEACAM6', 'CD274']] + 1)
    df = df[(df > 0).all(axis=1)]
    df = purity_adjust(df, ['CEACAM5', 'CEACAM6', 'CD274'])
    df['vital_status'] = vital_status(df.index)
    return df.reset_index(drop=True)


def plot_scatter(ax, x, y, vital, title, xlabel, ylabel, fontscale, mark,
                 area, marker_area, ylabel_y, left=True):
    r_val, p_val = stats.spearmanr(x, y)

    # Dead filled red, Alive open gray, Missing open light gray.
    for status, fc, ec, label in [
        ('Dead',    COLOR_DEAD,  COLOR_DEAD,    'Dead'),
        ('Alive',   'none',      COLOR_ALIVE,   'Alive'),
        ('Missing', 'none',      COLOR_MISSING, None),
    ]:
        mask = vital == status
        if mask.sum() == 0:
            continue
        ax.scatter(x[mask], y[mask], facecolors=fc, edgecolors=ec,
                   s=marker_area * area, alpha=0.85,
                   linewidths=style.EDGE_PT, zorder=3)

    # Regression line (above dots)
    valid = np.isfinite(x) & np.isfinite(y)
    xv, yv = np.array(x)[valid], np.array(y)[valid]
    if len(xv) >= 3 and np.std(xv) > 0:
        slope, intercept = np.polyfit(xv, yv, 1)
        x_line = np.linspace(xv.min(), xv.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=style.RULE_PT, alpha=0.6, zorder=5)

    ax.set_title(title, fontsize=style.tick_pt())
    # Headroom above the points for the rho and P label.
    lo, hi = np.nanmin(y), np.nanmax(y)
    ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.45 * (hi - lo))
    cnslayout.corr_annotate(ax, r_val, p=p_val, loc="upper left")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(style.RULE_PT)
    ax.tick_params(axis='both', width=style.RULE_PT, length=3*fontscale*mark)
    # Axis labels last: they are placed from the tick labels' extent.
    rich_xlabel(ax, xlabel)
    if left:
        if ylabel_y is None:
            rich_ylabel(ax, ylabel)
        else:
            rich_ylabel(ax, ylabel, y=ylabel_y)

    # Return the key handles; draw() places one shared key outside the boxes.
    key = 1.0 / plt.rcParams['legend.markerscale']
    h_dead = mlines.Line2D([], [], marker='o', color='none', markerfacecolor=COLOR_DEAD,
                           markeredgecolor=COLOR_DEAD,
                           markersize=4*fontscale/SCALE*mark*key,
                           label='Dead')
    h_alive = mlines.Line2D([], [], marker='o', color='none', markerfacecolor='none',
                            markeredgecolor=COLOR_ALIVE,
                            markersize=4*fontscale/SCALE*mark*key,
                            markeredgewidth=style.EDGE_PT, label='Alive')
    return r_val, p_val, [h_dead, h_alive]


def draw(panel, left_mm, combos, name, marker_area, loader, ylabel_y=None):
    """One TCGA-STAD scatter pair drawn from the table `loader` builds.

    `combos` is (x column, y column, x label, y label) per axes.
    """
    df = cache.table(panel.dir, 'tcga_samples', loader)

    results = []
    fig = panel.canvas()
    axes = cnslayout.scatter_pair_mm(fig, left_mm=left_mm)
    for col_idx, (xcol, ycol, xlabel, ylabel) in enumerate(combos):
        ax = axes[col_idx]
        sub = df
        r_val, p_val, keys = plot_scatter(ax, sub[xcol].values, sub[ycol].values,
                                          sub['vital_status'].values,
                                          'TCGA-STAD',
                                          xlabel, ylabel, SCALE,
                                          panel.mark, panel.area,
                                          marker_area, ylabel_y,
                                          left=(col_idx == 0))
        results.append({'x': xcol, 'y': ycol, 'r': r_val, 'p': p_val,
                        'n': len(sub)})
    corr_stats.write(panel.dir, [(r['x'], r['r'], r['p'], r['n'])
                                 for r in results])

    # Key to the right of the second box.
    group_key(fig, [(keys[0].get_label(), COLOR_DEAD),
                    (keys[1].get_label(), None, COLOR_ALIVE)],
              x_mm=left_mm + 2 * cnslayout.SCATTER_BOX_MM + cnslayout.SCATTER_GAP_MM + 0.8,
              y_mm=cnslayout.SCATTER_TOP_MM)

    return panel.finish(fig, name)


def main():
    draw(PANEL, LEFT_MM,
         [('CEACAM5', 'CD274', '*CEACAM5*', 'Epi. *CD274*\n(purity-adj.)'),
          ('CEACAM6', 'CD274', '*CEACAM6*', 'Epi. *CD274*\n(purity-adj.)')],
         'tcga_ceacam_cd274', MARKER_AREA, load_data_3b)


if __name__ == '__main__':
    main()

