#!/usr/bin/env python3
"""Supplementary Figure S2 A: canonical marker dot plot across twelve major cell
types (CD4+ T, CD8+ T and NK merged into T/NK).

Reads data/dot_frames.csv (scanpy per-group means, standard_scale='var', and
fractions); the atlas h5ad is opened only if that table is missing or with
--recompute. Writes the panel through _driver_base.save.
Run: python create_S2_A_marker_dotplot_12types.py [--recompute]

"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import FULL_DATASET_H5AD, SOUND_SOURCE       # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402
import _dotplot                                           # noqa: E402

FIG, PANEL = "S2_Cell_Annotation", "S2_A"

SOURCE = SOUND_SOURCE.get(FULL_DATASET_H5AD.name, FULL_DATASET_H5AD)


MARKERS = {
    'T/NK': ['CD3D', 'CD3E', 'NKG7', 'GNLY'],
    'B cells': ['CD79A', 'MS4A1'],
    'Plasma': ['JCHAIN', 'MZB1'],
    'DC': ['FCER1A', 'CLEC10A'],
    'Mono/Mac': ['CD14', 'LYZ'],
    'Neutrophils': ['FCGR3B', 'CSF3R'],
    'Mast cells': ['TPSAB1', 'KIT'],
    'Epithelial': ['EPCAM', 'KRT18'],
    'Hepatocyte': ['ALB', 'APOA1'],
    'Endothelial': ['PECAM1', 'VWF'],
    'Fibroblast': ['COL1A1', 'DCN'],
    'Pericyte': ['RGS5', 'ACTA2'],
}
GENES = [g for gs in MARKERS.values() for g in gs]
CELL_TYPE_MAP = {
    'CD4+ T cells': 'T/NK cells', 'CD8+ T cells': 'T/NK cells',
    'NK cells': 'T/NK cells', 'B cells': 'B cells', 'Plasma cells': 'Plasma cells',
    'DC cells': 'DC cells', 'Monocytes/Macrophages': 'Monocytes/Macrophages',
    'Neutrophils': 'Neutrophils', 'Mast cells': 'Mast cells',
    'Epithelial cells': 'Epithelial cells', 'Hepatocyte': 'Hepatocyte',
    'Endothelial cells': 'Endothelial cells', 'Fibroblast': 'Fibroblast',
    'Pericyte': 'Pericyte',
}
ORDER = ['T/NK cells', 'B cells', 'Plasma cells', 'DC cells',
         'Monocytes/Macrophages', 'Neutrophils', 'Mast cells', 'Epithelial cells',
         'Hepatocyte', 'Endothelial cells', 'Fibroblast', 'Pericyte']

# Panel size in mm.
W, H = 160.0, 50.0


def compute_frames():
    import scanpy as sc
    adata = sc.read_h5ad(SOURCE)
    adata.obs['cell_type_12'] = adata.obs['major_cell_type'].map(CELL_TYPE_MAP)
    adata = adata[adata.obs['cell_type_12'].notna()].copy()
    adata.obs['cell_type_12'] = adata.obs['cell_type_12'].astype('category')
    return _dotplot.frames(adata, GENES, 'cell_type_12', ORDER,
                           path=SOURCE)

if __name__ == "__main__":
    base.apply_style()
    fr = cache.table(HERE, "dot_frames", compute_frames)
    fig, _ = _dotplot.draw(fr, ORDER, GENES, w_mm=W, h_mm=H, left_mm=28.5,
                           bottom_mm=10.5, top_mm=2.0, largest_mm=2.2,
                           panel="S2 A")
    base.save(fig, FIG, PANEL, "panel_S2_A")
