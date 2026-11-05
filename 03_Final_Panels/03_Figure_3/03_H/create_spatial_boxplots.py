#!/usr/bin/env python3
"""
Figure 3H: neighbourhood epithelial proportion in CEACAM5/6-high versus
CEACAM5/6-low spots, one mean per section, paired within section. Drawn with
cnsfig.spatial.paired_boxplot (shared with Figures 3I and 3J), which prints the
exact one-tailed Wilcoxon signed-rank P (high > low); the two-sided P is
reported in the legend, with the test and sample count.

Reads SPATIAL_SPOT_DATA (00_Config/paths.py), column neighborhood_epi_density.
Writes ceacam_spatial_boxplot.[svg|pdf|png] beside this script (Figures 3I and
3J write the same file name in their own directories).

Run: python create_spatial_boxplots.py
"""

import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA
from cnsfig import Panel                                  # noqa: E402
from cnsfig.spatial import paired_boxplot, section_means  # noqa: E402

PANEL = Panel(3, "H", __file__, small_pt=7.0)


def main():
    df = pd.read_csv(SPATIAL_SPOT_DATA)
    sample_data = section_means(df, 'neighborhood_epi_density')
    paired_boxplot(PANEL, sample_data, 'neighborhood_epi_density',
                   'Epithelial Density',
                   'Neighborhood\nEpithelial Proportion',
                   'ceacam_spatial_boxplot', y_max_limit=1.0,
                   tight_y=True,
                   # The title band already clears the panel letter.
                   fit_kw=dict(pad_mm=0.3, reserve_letter=False),
                   exact=True)  # exact P, not a star


if __name__ == '__main__':
    main()
