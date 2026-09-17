#!/usr/bin/env python3
"""
S3 panel B - the stomach CD8+ T-cell UMAP coloured by PDCD1 (PD-1)
expression, drawn at the size it prints at.

  printed panel  Supplementary Figure S4 B   (PROVENANCE.csv)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/S3_CD8_TCells/
                 S4_B/create_S3_B_cd8_pdcd1_umap.py  (sc.pl.umap on an 8 x 7 cm
                 canvas at 4x)

DRAWING READS A TABLE (cnsfig.cache): data/cells.csv (x, y, value per stomach
CD8 cell; expression from .raw as scanpy reads it); the CD8 h5ad is opened
only when the table is absent or with --recompute.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import TCD8_H5AD                               # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402
import _feature_umap                                      # noqa: E402

FIG, PANEL = "S4_CD8_TCells", "S4_B"
GENE, TITLE, CMAP, VMIN = 'PDCD1', 'PDCD1 (PD-1) expression', 'Reds', 0

# The printed box, millimetres: one of three in row 2.
W, H = 53.0, 46.0
MAP_MM, LEFT_MM, TOP_MM = 36.0, 4.5, 3.6


def compute_cells():
    import scanpy as sc
    adata = sc.read_h5ad(TCD8_H5AD)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    src = adata.raw if adata.raw is not None else adata
    col = src[:, GENE].X
    values = col.toarray().ravel() if hasattr(col, 'toarray') else col.ravel()
    return _feature_umap.frame(adata, values)


def main():
    base.apply_style()
    fr = cache.table(HERE, "cells", compute_cells)
    fig, _ = _feature_umap.draw(fr, title=TITLE, cmap=CMAP, vmin=VMIN, w_mm=W, h_mm=H,
                                map_mm=MAP_MM, left_mm=LEFT_MM, top_mm=TOP_MM)
    base.save(fig, FIG, PANEL, f"panel_{PANEL}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
