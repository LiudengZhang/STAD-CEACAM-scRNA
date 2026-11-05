#!/usr/bin/env python3
"""Supplementary Figure S5 D: stomach CD8+ T-cell UMAP coloured by the Tex
exhaustion score (scanpy score_genes over TEX_GENES, default random_state).

Reads data/cells.csv; TCD8_H5AD is opened only if that table is missing or with
--recompute. Drawing: _drivers/_feature_umap.py.
Run: python create_S5_D_cd8_tex_score_umap.py [--recompute]"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import TCD8_H5AD                               # noqa: E402
import _feature_umap as fu                                # noqa: E402

TITLE, CMAP, VMIN = 'Tex exhaustion score', 'RdYlBu_r', None
TEX_GENES = ['PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'CTLA4',
             'TOX', 'ENTPD1', 'LAYN', 'CXCL13', 'BATF']

def compute_cells():
    import scanpy as sc
    adata = fu.stomach(sc.read_h5ad(TCD8_H5AD))
    available = [g for g in TEX_GENES if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=available, score_name='Tex_score')
    return fu.frame(adata, adata.obs['Tex_score'].to_numpy())

if __name__ == "__main__":
    sys.exit(fu.panel(HERE, "S5_D", compute_cells, title=TITLE, cmap=CMAP,
                      vmin=VMIN))
