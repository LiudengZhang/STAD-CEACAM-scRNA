#!/usr/bin/env python3
"""Supplementary Figure S6 D: the ten GSE251950 sections coloured by
distance to stroma (µm).

Reads SPATIAL_SPOT_DATA (per-spot table); drawing: _drivers/_spatial_gallery.py.
Run: python create_S6_D_distance_to_stroma.py"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import SPATIAL_SPOT_DATA                       # noqa: E402
from spatial_scale import cohort_um_per_unit               # noqa: E402
import _spatial_gallery as gallery                        # noqa: E402

if __name__ == "__main__":
    # Convert spot-grid distances to micrometres (spatial_scale.py).
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    df['distance_to_stroma'] = df['distance_to_stroma'] * cohort_um_per_unit(df)
    sys.exit(gallery.panel(df, "S6_D", 'distance_to_stroma', 'viridis', 'Distance to\nstroma (µm)'))
