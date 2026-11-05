#!/usr/bin/env python3
"""Supplementary Figure S5 A: marker dot plot of the CD8+ T-cell sub-clusters
(stomach samples).

Reads data/dot_frames.csv; TCD8_H5AD is opened only if that table is missing or
with --recompute. Run: python create_S5_A_cd8_dotplot.py [--recompute]

"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import TCD8_H5AD                               # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402
import _dotplot                                           # noqa: E402

FIG, PANEL = "S5_CD8_TCells", "S5_A"

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

# Row names as in Figure 3C; LABEL_MAP keys the cached table.
SHOWN = {
    'Cytotoxic CCL': 'TCD8_CCL', 'Cytotoxic DUSP1': 'TCD8_DUSP1',
    'MAIT': 'TCD8_MAIT', 'Tcm': 'TCD8_Tcm', 'TEMRA': 'TCD8_Temra',
    'Proliferating': 'TCD8_Prolif', 'Tex': 'TCD8_Tex', 'ISG': 'TCD8_ISG',
}

# Panel size in mm.
W, H = 160.0, 46.0

def compute_frames():
    import scanpy as sc
    adata = sc.read_h5ad(TCD8_H5AD)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    states = adata.obs['minor_cell_state'].astype(str)
    adata.obs['cluster_label'] = states.map(LABEL_MAP).fillna(states)
    return _dotplot.frames(adata, GENES, 'cluster_label', ORDER,
                           path=TCD8_H5AD)

if __name__ == "__main__":
    base.apply_style()
    fr = cache.table(HERE, "dot_frames", compute_frames)
    fr = fr.assign(group=fr["group"].map(SHOWN))
    fig, _ = _dotplot.draw(fr, [SHOWN[o] for o in ORDER], GENES, w_mm=W, h_mm=H, left_mm=22.0,
                           bottom_mm=11.0, top_mm=2.0, largest_mm=2.4, panel="S5 A")
    base.save(fig, FIG, PANEL, "panel_S5_A")
