#!/usr/bin/env python3
"""
Figure 3F: CEACAM5/6 ratio across the spots of one Visium section
(sample_03), drawn with cnsfig.spatial.spot_map (shared with Figures 3G and
3K-M).

Reads SPATIAL_SPOT_DATA (00_Config/paths.py), column CEACAM_ratio.
Writes spatial_ceacam_ratio.[svg|pdf|png] beside this script.

Run: python create_spatial_ceacam_ratio.py
"""

import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA
from cnsfig import Panel                                  # noqa: E402
from cnsfig.spatial import spot_map                       # noqa: E402

# A 1 mm gutter keeps the centred title clear of the panel letter.
PANEL = Panel(3, "F", __file__, small_pt=6.0, letter_gutter_mm=1.0)
SAMPLE_NAME = 'sample_03'


def main():
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    sample_df = df[df['sample'] == SAMPLE_NAME].copy()
    spot_map(PANEL, sample_df, 'CEACAM_ratio', 'RdBu_r',
             'CEACAM5/6 Ratio', 'spatial_ceacam_ratio',
             )


if __name__ == '__main__':
    main()
