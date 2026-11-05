#!/usr/bin/env python3
"""Supplementary Figure S6 E: the ten GSE251950 sections coloured by
distance to immune cells (µm).

Reads SPATIAL_SPOT_DATA (per-spot table); drawing: _drivers/_spatial_gallery.py.
Run: python create_S6_E_distance_to_immune.py"""

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
    df['distance_to_immune'] = df['distance_to_immune'] * cohort_um_per_unit(df)
    sys.exit(gallery.panel(df, "S6_E", 'distance_to_immune', 'plasma', 'Distance to\nimmune (µm)'))
