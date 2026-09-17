#!/usr/bin/env python3
"""
S4 panel C - the ten GSE251950 sections coloured by Total_Epi, drawn at the
size they print at.

  printed panel  Supplementary Figure S5 C   (PROVENANCE.csv)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/S4_Spatial_Validation/
                 S5_C/create_S4_C_epithelial_density_all_samples.py

Drawn by _drivers/_spatial_gallery.py from SPATIAL_SPOT_DATA, the per-spot
table the predecessor read; see that module for the labels.
"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import SPATIAL_SPOT_DATA                       # noqa: E402
import _driver_base as base                               # noqa: E402
import _spatial_gallery                                   # noqa: E402

FIG, PANEL = "S5_Spatial_Validation", "S5_C"
COLUMN, CMAP, CBAR_LABEL = 'Total_Epi', 'Greens', 'Total epithelial\nfraction'

# The printed box, millimetres: two galleries to a row.
W, H = 83.0, 36.0


def main():
    base.apply_style()
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    fig = _spatial_gallery.draw(df, COLUMN, cmap=CMAP, cbar_label=CBAR_LABEL,
                                w_mm=W, h_mm=H, panel="S5 C")
    base.save(fig, FIG, PANEL, f"panel_{PANEL}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
