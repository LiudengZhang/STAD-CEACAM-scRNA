#!/usr/bin/env python3
"""
Figure 3K: distance to stroma across the spots of one Visium section
(sample_03), drawn with cnsfig.spatial.spot_map. Values are converted to
micrometres with the cohort-wide factor from 00_Config/spatial_scale.py, as
in Figures 3I and 3J; the colour bar is labelled in um.

Reads SPATIAL_SPOT_DATA (00_Config/paths.py), column distance_to_stroma.
Writes spatial_dist_stroma.[svg|pdf|png] beside this script.

Run: python create_spatial_dist_stroma.py
"""

import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA
from spatial_scale import cohort_um_per_unit               # noqa: E402
from cnsfig import Panel                                  # noqa: E402
from cnsfig.spatial import spot_map                       # noqa: E402

PANEL = Panel(3, "K", __file__, small_pt=6.0, letter_gutter_mm=1.0)
SAMPLE_NAME = 'sample_03'


def main():
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    sample_df = df[df['sample'] == SAMPLE_NAME].copy()
    um_per_unit = cohort_um_per_unit(df)
    sample_df['distance_to_stroma'] = (
        sample_df['distance_to_stroma'] * um_per_unit)
    spot_map(PANEL, sample_df, 'distance_to_stroma', 'viridis',
             # Unit on the colour bar.
             'Distance to Stroma', 'spatial_dist_stroma',
             under='F', cbar_label='µm')


if __name__ == '__main__':
    main()
