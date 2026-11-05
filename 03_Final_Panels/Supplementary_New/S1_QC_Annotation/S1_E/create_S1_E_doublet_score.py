#!/usr/bin/env python3
"""Supplementary Figure S1 E - the Scrublet doublet-score distribution, singlets
against predicted doublets.

Reads data/histogram.csv: two density histograms over 60 bins on [0, 0.3] and
the two cell counts. The h5ad is read only when that table is absent or with
--recompute. ax.hist() is given the bin densities as weights.

Run: python create_S1_E_doublet_score.py [--recompute]"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import FULL_DATASET_H5AD                       # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S1_QC_Annotation", "S1_E"
BINS = np.linspace(0, 0.3, 61)

# Panel size and margins in mm.
W, H = 50.0, 40.0
MARGIN = dict(left=9.0, right=2.5, top=4.0, bottom=8.0)

def compute_hist():
    import scanpy as sc
    adata = sc.read_h5ad(FULL_DATASET_H5AD, backed="r")
    score = adata.obs["doublet_score"].to_numpy(dtype=float)
    doublet = adata.obs["predicted_doublet"].to_numpy(dtype=bool)
    d_single, _ = np.histogram(score[~doublet], bins=BINS, density=True)
    d_double, _ = np.histogram(score[doublet], bins=BINS, density=True)
    return pd.DataFrame({
        "bin_left": BINS[:-1], "bin_right": BINS[1:],
        "singlet_density": d_single, "doublet_density": d_double,
        "n_singlet": int((~doublet).sum()), "n_doublet": int(doublet.sum()),
    })

def draw(hist):
    base.apply_style()
    fig, ax = style.subplots_mm(W, H)
    style.margins_mm(fig, **MARGIN)
    centres = (hist["bin_left"] + hist["bin_right"]).to_numpy() / 2
    ax.hist(centres, bins=BINS, weights=hist["singlet_density"], alpha=0.7,
            color="#4daf4a")
    ax.hist(centres, bins=BINS, weights=hist["doublet_density"], alpha=0.7,
            color="#e41a1c")
    ax.set_xlabel("Doublet Score")
    ax.set_ylabel("Density")
    ax.set_title("Doublet Score Distribution")
    base.despine(ax)
    ax.set_xlim(0.0, 0.3)
    return fig

if __name__ == "__main__":
    base.save(draw(cache.table(HERE, "histogram", compute_hist)),
              FIG, PANEL, "panel_S1_E")
