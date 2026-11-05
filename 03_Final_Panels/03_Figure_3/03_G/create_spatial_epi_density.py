#!/usr/bin/env python3
"""
Figure 3G: neighbourhood epithelial density across the spots of one Visium
section (sample_03), drawn with cnsfig.spatial.spot_map.

Reads SPATIAL_SPOT_DATA (00_Config/paths.py), column neighborhood_epi_density.
Writes spatial_epi_density.[svg|pdf|png] beside this script.

Run: python create_spatial_epi_density.py
"""

import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA
from cnsfig import Panel                                  # noqa: E402
from cnsfig.spatial import spot_map                       # noqa: E402

PANEL = Panel(3, "G", __file__, small_pt=6.0, letter_gutter_mm=1.0)
SAMPLE_NAME = 'sample_03'


def main():
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    sample_df = df[df['sample'] == SAMPLE_NAME].copy()
    spot_map(PANEL, sample_df, 'neighborhood_epi_density', 'YlOrRd',
             'Epithelial Density', 'spatial_epi_density',
             )


if __name__ == '__main__':
    main()
