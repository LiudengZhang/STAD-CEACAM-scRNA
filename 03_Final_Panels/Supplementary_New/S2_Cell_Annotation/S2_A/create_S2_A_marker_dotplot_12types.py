#!/usr/bin/env python3
"""
S2 panel A (S1 F until 2026-09-16, when the author split S1) - canonical marker dot plot across the twelve major cell types
(CD4+ T, CD8+ T and NK merged into T/NK), drawn at the size it prints at.

  printed panel  Supplementary Figure S2 A   (PROVENANCE.csv)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/S1_QC_Annotation/
                 S1_F/create_S1_F_marker_dotplot_12types.py  (a 4x scanpy
                 dot plot, fitted into the page at ~0.25)

DRAWING READS A TABLE (cnsfig.cache): data/dot_frames.csv holds scanpy's own
per-group means (standard_scale='var') and fractions for the 24 markers; the
full-dataset h5ad is opened only when the table is absent or with
--recompute. _drivers/_dotplot.py draws the frames through scanpy's
dot_color_df / dot_size_df.
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

FIG, PANEL = "S2_Cell_Annotation", "S2_A"

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

# The printed box, millimetres: 26 genes at a 2.48 mm pitch (rotated 6 pt
# gene names stand under check_restyled_panel's 0.5 pt clearance below 2.45) plus the row names ("Monocytes/Macrophages",
# 27.6 mm). No key column: F shares its row with G since the evening of
# 2026-09-15 (S1 on two pages, each under the main-figure page height), and
# G's key at the right of the row decodes both plots - the two are scaled
# alike (_dotplot.draw, shared_key). Until then F had a row to itself at
# 116 x 50 mm with a key of its own.
W, H = 93.0, 50.0

#: The panel whose key decodes this plot. It is the dot plot to the RIGHT on
#: the same row - today's S2 B. On 2026-09-16 the author split S1 in two
#: (QC A-E stays S1; the annotation panels F, G, H become S2 A, B, C), so the
#: printed name of that neighbour is "S2 B"; Wave 2 of that round renamed
#: the directories to match. The name is only printed to the console by
#: _dotplot.draw - nothing on the page carries it.
SHARED_KEY_WITH = "S2 B"


def compute_frames():
    import scanpy as sc
    adata = sc.read_h5ad(FULL_DATASET_H5AD)
    adata.obs['cell_type_12'] = adata.obs['major_cell_type'].map(CELL_TYPE_MAP)
    adata = adata[adata.obs['cell_type_12'].notna()].copy()
    adata.obs['cell_type_12'] = adata.obs['cell_type_12'].astype('category')
    return _dotplot.frames(adata, GENES, 'cell_type_12', ORDER)


def main():
    base.apply_style()
    fr = cache.table(HERE, "dot_frames", compute_frames)
    fig, _ = _dotplot.draw(fr, ORDER, GENES, w_mm=W, h_mm=H, left_mm=28.5,
                           bottom_mm=10.5, top_mm=2.0, largest_mm=2.2, panel="S2 A",
                           shared_key=SHARED_KEY_WITH)
    base.save(fig, FIG, PANEL, "panel_S2_A")
    return 0


if __name__ == "__main__":
    sys.exit(main())
