#!/usr/bin/env python3
"""
S3 panel A - CD8+ T-cell sub-cluster marker dot plot (stomach), drawn at the
size it prints at.

  printed panel  Supplementary Figure S4 A   (PROVENANCE.csv)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/S3_CD8_TCells/
                 S4_A/create_S3_A_cd8_dotplot.py

DRAWING READS A TABLE (cnsfig.cache): data/dot_frames.csv, as S2 A.
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
import _dotplot                                           # noqa: E402

FIG, PANEL = "S4_CD8_TCells", "S4_A"

MARKERS = {
    'Cytotoxic CCL': ['CCL4', 'CCL5'], 'Cytotoxic DUSP1': ['DUSP1', 'JUN'],
    'MAIT': ['KLRB1', 'SLC4A10'], 'Tcm': ['CCR7', 'SELL'],
    'TEMRA': ['KLRG1', 'CX3CR1'], 'Proliferating': ['MKI67', 'TOP2A'],
    'Tex': ['PDCD1', 'HAVCR2', 'LAG3'], 'ISG': ['ISG15', 'MX1'],
}
GENES = [g for gs in MARKERS.values() for g in gs]
LABEL_MAP = {
    'C0_CD8_Cytotoxic_CCL': 'Cytotoxic CCL', 'C1_CD8_Cytotoxic_DUSP1': 'Cytotoxic DUSP1',
    'C2_CD8_MAIT_KLRB1': 'MAIT', 'C3_CD8_Tcm_CCR7': 'Tcm', 'C4_CD8_Temra_KLRG1': 'TEMRA',
    'C5_CD8_Prolif_MKI67': 'Proliferating', 'C6_CD8_Tex_PDCD1': 'Tex',
    'C7_CD8_ISG_ISG15': 'ISG',
}
ORDER = list(LABEL_MAP.values())

# The printed box, millimetres: alone in row 1, and since 2026-09-16 across
# it (the author's fifth reading: "S3A can also stretch across a row";
# 100 mm until then, 17 genes at 3 mm).
W, H = 160.0, 46.0


def compute_frames():
    import scanpy as sc
    adata = sc.read_h5ad(TCD8_H5AD)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    states = adata.obs['minor_cell_state'].astype(str)
    adata.obs['cluster_label'] = states.map(LABEL_MAP).fillna(states)
    return _dotplot.frames(adata, GENES, 'cluster_label', ORDER)


def main():
    base.apply_style()
    fr = cache.table(HERE, "dot_frames", compute_frames)
    fig, _ = _dotplot.draw(fr, ORDER, GENES, w_mm=W, h_mm=H, left_mm=22.0,
                           bottom_mm=11.0, top_mm=2.0, largest_mm=2.4, panel="S4 A")
    base.save(fig, FIG, PANEL, "panel_S4_A")
    return 0


if __name__ == "__main__":
    sys.exit(main())
