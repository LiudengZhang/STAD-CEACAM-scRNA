#!/usr/bin/env python3
"""
Figure 3J: distance to immune cells from CEACAM5/6-high and CEACAM5/6-low
spots, one mean per section, paired within section. Drawn with
cnsfig.spatial.paired_boxplot (exact one-tailed Wilcoxon signed-rank P on the
panel, two-sided in the legend). Distances are converted to micrometres with
the cohort-wide factor from 00_Config/spatial_scale.py (see Figure 3I).

Reads SPATIAL_SPOT_DATA (00_Config/paths.py), column distance_to_immune.
Writes ceacam_spatial_boxplot.[svg|pdf|png] beside this script.

Run: python create_spatial_boxplot_immune.py
"""

import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA
from spatial_scale import cohort_um_per_unit              # noqa: E402
from cnsfig import Panel                                  # noqa: E402
from cnsfig.spatial import paired_boxplot, section_means  # noqa: E402

PANEL = Panel(3, "J", __file__, small_pt=7.0)


def main():
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    sample_data = section_means(df, 'distance_to_immune')
    um_per_unit = cohort_um_per_unit(df)
    sample_data['distance_to_immune'] = (
        sample_data['distance_to_immune'] * um_per_unit)
    paired_boxplot(PANEL, sample_data, 'distance_to_immune',
                   'Distance to Immune', 'Distance (µm)',
                   'ceacam_spatial_boxplot',
                   y_max_limit=1200 * um_per_unit,
                   exact=True)  # exact P, not a star


if __name__ == '__main__':
    main()
