#!/usr/bin/env python3
"""
S1 panel E - the Scrublet doublet-score distribution, singlets against
predicted doublets, drawn at the size it prints at.

  printed panel  Supplementary Figure S1 E   (PROVENANCE.csv)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/S1_QC_Annotation/
                 S1_E/create_S1_E_doublet_score.py  (8 x 6 cm at 4x)

DRAWING READS A TABLE (cnsfig.cache): data/histogram.csv holds the two
density histograms over the predecessor's 60 bins on [0, 0.3] and the two
cell counts (the key that printed them was dropped on 2026-09-16; the
figure legend carries the count); the h5ad's 542,121 doublet scores are read
only when the table is absent or with --recompute. ax.hist() is given the
bin densities as weights, which draws the same rectangles it drew from the
scores.
"""

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

FIG, PANEL = "S1_QC_Annotation", "S1_E"
BINS = np.linspace(0, 0.3, 61)

# The printed box, millimetres: beside D's two maps. 40 mm tall since
# 2026-09-16 (S1 back on one page).
W, H = 50.0, 40.0
MARGIN = dict(left=9.0, right=2.5, top=4.0, bottom=8.0)


def compute_hist():
    """The computing half: the two density histograms and their counts."""
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
    n_single, n_double = int(hist["n_singlet"][0]), int(hist["n_doublet"][0])
    # NO KEY, since 2026-09-16 (the author's fifth reading). The predecessor
    # keyed "Singlets (n=542,121)" and "Doublets (n=0)"; the published S1 E
    # carried no key, and a red "Doublets (n=0)" beside an empty histogram
    # invites a question the legend answers (the object holds no predicted
    # doublets). The two histograms are drawn as before; the strings are
    # declared removed in labels.py:REMOVALS_S1_S6. n_single/n_double stay
    # in the table for the legend's count.
    del n_single, n_double
    ax.hist(centres, bins=BINS, weights=hist["singlet_density"], alpha=0.7,
            color="#4daf4a")
    ax.hist(centres, bins=BINS, weights=hist["doublet_density"], alpha=0.7,
            color="#e41a1c")
    ax.set_xlabel("Doublet Score")
    ax.set_ylabel("Density")
    ax.set_title("Doublet Score Distribution")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0.0, 0.3)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"ink outside the {W} x {H} mm canvas: {over}")
    return fig


def main():
    hist = cache.table(HERE, "histogram", compute_hist)
    base.save(draw(hist), FIG, PANEL, "panel_S1_E")
    return 0


if __name__ == "__main__":
    sys.exit(main())
