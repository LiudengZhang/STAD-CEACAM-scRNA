#!/usr/bin/env python3
"""Supplementary Figure S6 A: the ten GSE251950 sections coloured by
the CEACAM5/6 ratio.

Reads SPATIAL_SPOT_DATA (per-spot table); drawing: _drivers/_spatial_gallery.py.
Run: python create_S6_A_ceacam_ratio_all_samples.py"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import SPATIAL_SPOT_DATA                       # noqa: E402
import _spatial_gallery as gallery                        # noqa: E402

if __name__ == "__main__":
    sys.exit(gallery.panel(pd.read_csv(SPATIAL_SPOT_DATA), "S6_A",
                           'CEACAM_ratio', 'RdBu_r', 'CEACAM5/6 ratio\n(high/total epi)'))
