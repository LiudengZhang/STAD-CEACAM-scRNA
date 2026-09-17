#!/usr/bin/env python3
"""
S2 panel A - epithelial sub-cluster marker dot plot (stomach), drawn at the
size it prints at.

  printed panel  Supplementary Figure S3 A   (PROVENANCE.csv)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/
                 S2_CEACAM_Metaprogram_Validation/S2_A/create_S2_A_epithelial_dotplot.py

DRAWING READS A TABLE (cnsfig.cache): data/dot_frames.csv (scanpy's own
standard_scale='var' means and fractions); the epithelial h5ad is opened only
when the table is absent or with --recompute. The predecessor read
Epithelial_tumor_scored.h5ad, which is Epithelial.h5ad with obs['tumor_score']
added and exists in the input tree only; this panel needs no score and reads
the deposited object (EPITHELIAL_DEPOSIT_H5AD), the same cells and the same
.raw - the comparator against the predecessor says so.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import EPITHELIAL_DEPOSIT_H5AD                 # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402
import _dotplot                                           # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_A"

MARKERS = {
    'PTMA': ['PTMA', 'STMN1'], 'KRT19': ['KRT19', 'KRT8'],
    'CEACAM5/6': ['CEACAM6', 'CEACAM5'], 'Chief-like': ['PGC', 'LIPF'],
    'MUC5AC': ['MUC5AC', 'TFF1'], 'Stem TPX2': ['TPX2', 'MKI67'],
    'Stem SPINK4': ['SPINK4', 'SPINK1'], 'MT1E': ['MT1E', 'MT2A'],
    'CD74': ['CD74', 'HLA-DRA'],
}
#: The twelve of the eighteen markers the printed panel carries. The
#: predecessor kept the markers present in Epithelial_tumor_scored.h5ad's
#: processed matrix (a highly-variable subset), which dropped PTMA, KRT19,
#: KRT8, SPINK1, MT1E and CD74; the deposited object's matrix is wider, so the
#: printed set is named here rather than rediscovered from whichever matrix is
#: on disk (00_GROUND_TRUTH/figures/S3_*.pdf is the authority).
GENES = ['STMN1', 'CEACAM6', 'CEACAM5', 'PGC', 'LIPF', 'MUC5AC', 'TFF1', 'TPX2',
         'MKI67', 'SPINK4', 'MT2A', 'HLA-DRA']
LABEL_MAP = {
    'C0_Epi_PTMA': 'PTMA', 'C1_Epi_KRT19': 'KRT19', 'C2_Epi_CEACAM6': 'CEACAM5/6',
    'C3_Epi_Chief_Like_PGC': 'Chief-like', 'C4_Epi_MUC5AC': 'MUC5AC',
    'C5_Epi_Stem_Like_TPX2': 'Stem TPX2', 'C6_Epi_Stem_Like_SPINK4': 'Stem SPINK4',
    'C7_Epi_MT1E': 'MT1E', 'C8_Epi_CD74': 'CD74',
}
ORDER = list(LABEL_MAP.values())

# The printed box, millimetres: a row to itself since 2026-09-16 (the
# author's fifth reading: "A can take almost a whole row"; G moved down into
# letter order). 105 x 48 until then.
W, H = 160.0, 46.0


def compute_frames():
    import scanpy as sc
    adata = sc.read_h5ad(EPITHELIAL_DEPOSIT_H5AD)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    states = adata.obs['minor_cell_state'].astype(str)
    adata.obs['cluster_label'] = states.map(LABEL_MAP).fillna(states)
    return _dotplot.frames(adata, GENES, 'cluster_label', ORDER)


def main():
    base.apply_style()
    fr = cache.table(HERE, "dot_frames", compute_frames)
    fig, _ = _dotplot.draw(fr, ORDER, GENES, w_mm=W, h_mm=H, left_mm=18.0,
                           bottom_mm=11.5, top_mm=2.0, largest_mm=2.4, panel="S3 A")
    base.save(fig, FIG, PANEL, "panel_S3_A")
    return 0


if __name__ == "__main__":
    sys.exit(main())
