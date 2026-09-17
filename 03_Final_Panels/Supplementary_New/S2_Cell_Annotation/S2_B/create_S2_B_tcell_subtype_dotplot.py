#!/usr/bin/env python3
"""
S2 panel B (S1 G until 2026-09-16, when the author split S1) - T-cell subtype marker dot plot (CD4+ T, CD8+ T, NK), drawn at
the size it prints at.

  printed panel  Supplementary Figure S2 B   (PROVENANCE.csv)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/S1_QC_Annotation/
                 S1_G/create_S1_G_tcell_subtype_dotplot.py

DRAWING READS A TABLE (cnsfig.cache): data/dot_frames.csv, as S2 A.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import FULL_DATASET_H5AD                       # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402
import _dotplot                                           # noqa: E402

FIG, PANEL = "S2_Cell_Annotation", "S2_B"

MARKERS = {
    'Pan-T': ['CD3D', 'CD3E'],
    'CD4': ['CD4', 'IL7R'],
    'CD8': ['CD8A', 'CD8B'],
    'NK': ['NKG7', 'GNLY', 'NCAM1'],
    'Treg': ['FOXP3', 'IL2RA'],
    'Cytotoxic': ['GZMB', 'PRF1'],
    'Exhaustion': ['PDCD1', 'HAVCR2', 'LAG3'],
}
GENES = [g for gs in MARKERS.values() for g in gs]
ORDER = ['CD4+ T cells', 'CD8+ T cells', 'NK cells']

# The printed box, millimetres: 16 genes at a 2.47 mm pitch, beside F on one
# row since the evening of 2026-09-15 (93 + 5 + 73 = 171 mm), with the key
# that decodes both plots (a panel that encodes a fraction in dot area
# decodes it; F reads this one - _dotplot.draw, shared_key). Until then G
# had a row to itself at 74 mm with a 2.5 mm pitch. H is F's height, so the
# two plots' gene names sit on one line.
W, H = 73.0, 50.0
#: The three rows sit low, with 24 mm of paper above them: the key column
#: beside them (size key over colour bar) is 25 mm tall and hangs from the
#: matrix top, so the rows cannot sit lower; and F's twelve rows fill the
#: same 50 mm, so the gene names of the two plots share a baseline.
TOP_MM = 24.0


def compute_frames():
    import scanpy as sc
    adata = sc.read_h5ad(FULL_DATASET_H5AD)
    adata = adata[adata.obs['major_cell_type'].isin(ORDER)].copy()
    return _dotplot.frames(adata, GENES, 'major_cell_type', ORDER)


def main():
    base.apply_style()
    fr = cache.table(HERE, "dot_frames", compute_frames)
    fig, _ = _dotplot.draw(fr, ORDER, GENES, w_mm=W, h_mm=H, left_mm=14.5,
                           bottom_mm=10.5, top_mm=TOP_MM, largest_mm=2.2,
                           panel="S2 B")
    base.save(fig, FIG, PANEL, "panel_S2_B")
    return 0


if __name__ == "__main__":
    sys.exit(main())
