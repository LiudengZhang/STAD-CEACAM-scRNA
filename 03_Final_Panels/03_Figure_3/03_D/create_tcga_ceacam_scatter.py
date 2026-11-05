#!/usr/bin/env python3
"""
Figure 3E: deconvolved epithelial CEACAM5 and CEACAM6 against CD8 T-cell
scarcity in TCGA-STAD (minus the mean z-score of log2(CPM + 1) of CD8A and
CD8B in the tumour bulk), one point per tumour, filled (dead) or open (alive)
by vital status. Tumours with at least 10% epithelium; every variable is
purity-adjusted. Drawn by draw() in 03_C/create_tcga_ceacam_scatter.py (a
different file with the same name, Figure 3B).

Reads data/tcga_samples.csv beside this script; with --recompute, computes it
from TCGA_BAYESPRISM_EPI, TCGA_BAYESPRISM_FRAC, TCGA_TUMOUR_CD8_COUNTS, the
sample list of TCGA_BULK_TUMOR_ONLY, ABSOLUTE purity and vital status.
Writes tcga_ceacam_tex.[svg|pdf|png] beside this script.

Run: python create_tcga_ceacam_scatter.py [--recompute]
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "03_C"))
from cnsfig import Panel                                       # noqa: E402
from create_tcga_ceacam_scatter import (draw, purity_adjust, vital_status,  # noqa: E402
                                        TCGA_EPI_EXPR)
import numpy as np                                             # noqa: E402
import pandas as pd                                            # noqa: E402
from paths import (TCGA_BAYESPRISM_FRAC, TCGA_BULK_TUMOR_ONLY,  # noqa: E402
                   TCGA_TUMOUR_CD8_COUNTS)

PANEL = Panel(3, "E", __file__, small_pt=4.5)
# Left margin, shared with Figure 3D above so the boxes align.
LEFT_MM = 12.5
# Marker area (about 0.42 mm across).
MARKER_AREA = 13.0


def load_data_3e():
    """Figure 3E: deconvolved epithelial CEACAM5/6 (log2 CPM + 1) against CD8
    T-cell scarcity, minus the mean z of log2(CPM + 1) of CD8A and CD8B in the
    tumour bulk; tumours with at least 10% epithelium; purity-adjusted."""
    epi = pd.read_csv(TCGA_EPI_EXPR, sep='\t', index_col=0)
    tum = pd.read_csv(TCGA_BULK_TUMOR_ONLY, sep='\t', index_col=0, nrows=1).columns
    frac = pd.read_csv(TCGA_BAYESPRISM_FRAC, sep='\t', index_col=0)['Epithelial cells']
    ids = [i for i in sorted(set(epi.index) & set(tum)) if frac.get(i, 0) >= 0.10]
    cd8 = pd.read_csv(TCGA_TUMOUR_CD8_COUNTS, sep='\t', index_col=0).loc[ids]
    cpm = np.log2(cd8[['CD8A', 'CD8B']].div(cd8['library_size'], axis=0) * 1e6 + 1)
    df = np.log2(epi.loc[ids, ['CEACAM5', 'CEACAM6']] + 1)
    df['cd8_scarcity'] = -((cpm - cpm.mean()) / cpm.std(ddof=1)).mean(axis=1)
    df = purity_adjust(df, ['CEACAM5', 'CEACAM6', 'cd8_scarcity'])
    df['vital_status'] = vital_status(df.index)
    return df.reset_index(drop=True)


def main():
    draw(PANEL, LEFT_MM,
         [('CEACAM5', 'cd8_scarcity', '*CEACAM5*', 'CD8 scarcity\n(−CD8 score)'),
          ('CEACAM6', 'cd8_scarcity', '*CEACAM6*', 'CD8 scarcity\n(−CD8 score)')],
         'tcga_ceacam_tex', MARKER_AREA, load_data_3e, ylabel_y=0.22)


if __name__ == '__main__':
    main()
