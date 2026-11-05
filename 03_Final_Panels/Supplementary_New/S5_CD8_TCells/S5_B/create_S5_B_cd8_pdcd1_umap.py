#!/usr/bin/env python3
"""Supplementary Figure S5 B: stomach CD8+ T-cell UMAP coloured by PDCD1 (PD-1)
expression.

Reads data/cells.csv (UMAP x, y and expression per cell); TCD8_H5AD is opened
only if that table is missing or with --recompute. Drawing: _drivers/_feature_umap.py.
Run: python create_S5_B_cd8_pdcd1_umap.py [--recompute]"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import TCD8_H5AD                               # noqa: E402
import _feature_umap as fu                                # noqa: E402

GENE, TITLE, CMAP, VMIN = "PDCD1", "PDCD1 (PD-1) expression", "Reds", 0

def compute_cells():
    import scanpy as sc
    adata = fu.stomach(sc.read_h5ad(TCD8_H5AD))
    return fu.frame(adata, fu.gene_values(adata, GENE, path=TCD8_H5AD))

if __name__ == "__main__":
    sys.exit(fu.panel(HERE, "S5_B", compute_cells, title=TITLE, cmap=CMAP,
                      vmin=VMIN))
