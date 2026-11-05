#!/usr/bin/env python3
"""
Figure 3L: distance to immune cells across the spots of one Visium section
(sample_03), in micrometres (cohort-wide factor from 00_Config/spatial_scale.py).

Reads SPATIAL_SPOT_DATA (00_Config/paths.py), column distance_to_immune.
Writes spatial_dist_immune.[svg|pdf|png] beside this script.

Run: python create_spatial_dist_immune.py
"""

import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA
from spatial_scale import cohort_um_per_unit               # noqa: E402
from cnsfig import Panel                                  # noqa: E402
from cnsfig.spatial import spot_map                       # noqa: E402

PANEL = Panel(3, "L", __file__, small_pt=6.0, letter_gutter_mm=1.0)
SAMPLE_NAME = 'sample_03'


def main():
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    sample_df = df[df['sample'] == SAMPLE_NAME].copy()
    um_per_unit = cohort_um_per_unit(df)
    sample_df['distance_to_immune'] = (
        sample_df['distance_to_immune'] * um_per_unit)
    spot_map(PANEL, sample_df, 'distance_to_immune', 'plasma',
             'Distance to Immune', 'spatial_dist_immune',
             under='G', cbar_label='µm')


if __name__ == '__main__':
    main()
