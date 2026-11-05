#!/usr/bin/env python3
"""Supplementary Figure S1 C - per-sample QC box plots of the 38 non-stomach
samples (genes per cell, total counts, mitochondrial %), coloured by tissue
site. Drawing is done by _drivers/_qc_boxstrip.py.

Reads data/box_stats.csv; the full-dataset h5ad is opened only when that table
is absent or with --recompute.

Run: python create_S1_C_qc_non_stomach.py [--recompute]"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import FULL_DATASET_H5AD                       # noqa: E402
import _qc_boxstrip as strip                              # noqa: E402

COLORS = {
    'Liver': '#e6ab02',
    'Peripheral blood': '#66a61e',
    'Lymph node': '#7570b3',
    'Ovary': '#e7298a',
    'Metastatic lymph node': '#a6761d',
}
GROUP_ORDER = ['Liver', 'Peripheral blood', 'Lymph node', 'Ovary',
               'Metastatic lymph node']

# Panel size and margins in mm.
W, H = 171.0, 46.0
MARGIN = dict(left=8.0, right=1.0, top=3.5, bottom=18.0)

def compute_stats():
    import scanpy as sc
    adata = sc.read_h5ad(FULL_DATASET_H5AD, backed='r')
    mask = adata.obs['Sample site'] != 'Stomach'
    df = adata.obs.loc[mask, ['sample', 'Sample ID', 'Sample site']
                       + strip.METRIC_COLS]
    info = df[['sample', 'Sample ID', 'Sample site']].drop_duplicates('sample').copy()
    info['group'] = info['Sample site'].astype(str)
    return strip.box_stats(df, info, GROUP_ORDER)

if __name__ == "__main__":
    sys.exit(strip.panel(HERE, "S1_C", compute_stats, w_mm=W, h_mm=H,
                         margin=MARGIN, colors=COLORS,
                         group_order=GROUP_ORDER, legend_all=False))
