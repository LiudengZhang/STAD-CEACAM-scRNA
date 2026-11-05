#!/usr/bin/env python3
"""Supplementary Figure S3 A: marker dot plot of the epithelial sub-clusters
(stomach samples).

Reads data/dot_frames.csv (scanpy per-group means, standard_scale='var', and
fractions); EPITHELIAL_DEPOSIT_H5AD is opened only if that table is missing or
with --recompute. Writes the panel through _driver_base.save.
Run: python create_S3_A_epithelial_dotplot.py [--recompute]

"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import EPITHELIAL_DEPOSIT_H5AD                 # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402
import _dotplot                                           # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_A"

GENES = ['STMN1', 'CEACAM6', 'CEACAM5', 'PGC', 'LIPF', 'MUC5AC', 'TFF1', 'TPX2',
         'MKI67', 'SPINK4', 'MT2A', 'HLA-DRA']
LABEL_MAP = {
    'C0_Epi_PTMA': 'PTMA', 'C1_Epi_KRT19': 'KRT19', 'C2_Epi_CEACAM6': 'CEACAM5/6',
    'C3_Epi_Chief_Like_PGC': 'Chief-like', 'C4_Epi_MUC5AC': 'MUC5AC',
    'C5_Epi_Stem_Like_TPX2': 'Stem TPX2', 'C6_Epi_Stem_Like_SPINK4': 'Stem SPINK4',
    'C7_Epi_MT1E': 'MT1E', 'C8_Epi_CD74': 'CD74',
}
ORDER = list(LABEL_MAP.values())

# Row names as in Figure 2A; LABEL_MAP keys the cached table.
SHOWN = {
    'PTMA': 'Epi_PTMA', 'KRT19': 'Epi_KRT19', 'CEACAM5/6': 'Epi_CEACAM5/6',
    'Chief-like': 'Epi_Chief_Like', 'MUC5AC': 'Epi_MUC5AC',
    'Stem TPX2': 'Epi_TPX2', 'Stem SPINK4': 'Epi_SPINK4',
    'MT1E': 'Epi_MT1E', 'CD74': 'Epi_CD74',
}

# Panel size in mm.
W, H = 160.0, 46.0

def compute_frames():
    import scanpy as sc
    adata = sc.read_h5ad(EPITHELIAL_DEPOSIT_H5AD)
    adata = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    states = adata.obs['minor_cell_state'].astype(str)
    adata.obs['cluster_label'] = states.map(LABEL_MAP).fillna(states)
    return _dotplot.frames(adata, GENES, 'cluster_label', ORDER,
                           path=EPITHELIAL_DEPOSIT_H5AD)

if __name__ == "__main__":
    base.apply_style()
    fr = cache.table(HERE, "dot_frames", compute_frames)
    fr = fr.assign(group=fr["group"].map(SHOWN))
    fig, _ = _dotplot.draw(fr, [SHOWN[o] for o in ORDER], GENES, w_mm=W, h_mm=H, left_mm=18.0,
                           bottom_mm=11.5, top_mm=2.0, largest_mm=2.4, panel="S3 A")
    base.save(fig, FIG, PANEL, "panel_S3_A")
