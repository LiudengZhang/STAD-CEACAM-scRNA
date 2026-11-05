#!/usr/bin/env python3
"""Supplementary Figure S2 B: T-cell subtype marker dot plot (CD4+ T, CD8+ T, NK).

Reads data/dot_frames.csv; the atlas h5ad is opened only if that table is
missing or with --recompute. Writes the panel through _driver_base.save.
Run: python create_S2_B_tcell_subtype_dotplot.py [--recompute]

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

FIG, PANEL = "S2_Cell_Annotation", "S2_B"

SOURCE = SOUND_SOURCE.get(FULL_DATASET_H5AD.name, FULL_DATASET_H5AD)


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

# Panel size in mm.
W, H = 160.0, 28.0
TOP_MM = 2.0

def compute_frames():
    import scanpy as sc
    adata = sc.read_h5ad(SOURCE)
    adata = adata[adata.obs['major_cell_type'].isin(ORDER)].copy()
    return _dotplot.frames(adata, GENES, 'major_cell_type', ORDER,
                           path=SOURCE)

if __name__ == "__main__":
    base.apply_style()
    fr = cache.table(HERE, "dot_frames", compute_frames)
    fig, _ = _dotplot.draw(fr, ORDER, GENES, w_mm=W, h_mm=H, left_mm=14.5,
                           bottom_mm=10.5, top_mm=TOP_MM, largest_mm=2.2,
                           panel="S2 B")
    base.save(fig, FIG, PANEL, "panel_S2_B")
