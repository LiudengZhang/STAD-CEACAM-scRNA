#!/usr/bin/env python3
"""Supplementary Figure S1 B - per-sample QC box plots of the 32 stomach samples
(genes per cell, total counts, mitochondrial %), coloured by treatment phase
and response. Drawing is done by _drivers/_qc_boxstrip.py.

Reads data/box_stats.csv; the full-dataset h5ad is opened only when that table
is absent or with --recompute.

Run: python create_S1_B_qc_metrics.py [--recompute]"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import FULL_DATASET_H5AD                       # noqa: E402
import _qc_boxstrip as strip                              # noqa: E402

from shared.figure_config import GROUP_COLORS  # noqa: E402
COLORS = {**GROUP_COLORS, 'Not selected': '#d9d9d9'}
GROUP_ORDER = ['Pre-R', 'Pre-NR', 'Post-R', 'Post-NR', 'Not selected']

# Panel size and margins in mm.
W, H = 171.0, 44.0
MARGIN = dict(left=8.0, right=1.0, top=3.5, bottom=15.5)

def _group(row):
    pre, post = str(row['stomach_pre_grouping']), str(row['stomach_post_grouping'])
    if pre == 'Responsed':
        return 'Pre-R'
    if pre == 'No-response':
        return 'Pre-NR'
    if post == 'Responsed':
        return 'Post-R'
    if post == 'No-response':
        return 'Post-NR'
    return 'Not selected'

def compute_stats():
    import scanpy as sc
    adata = sc.read_h5ad(FULL_DATASET_H5AD, backed='r')
    mask = adata.obs['Sample site'] == 'Stomach'
    df = adata.obs.loc[mask, ['sample', 'Sample ID', 'stomach_pre_grouping',
                              'stomach_post_grouping'] + strip.METRIC_COLS]
    info = df[['sample', 'Sample ID', 'stomach_pre_grouping',
               'stomach_post_grouping']].drop_duplicates('sample').copy()
    info['group'] = info.apply(_group, axis=1)
    return strip.box_stats(df, info, GROUP_ORDER)

if __name__ == "__main__":
    sys.exit(strip.panel(HERE, "S1_B", compute_stats, w_mm=W, h_mm=H,
                         margin=MARGIN, colors=COLORS,
                         group_order=GROUP_ORDER, legend_all=True))
