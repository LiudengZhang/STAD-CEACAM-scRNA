#!/usr/bin/env python3
"""Figure 5J, box 4 of 4: CD274 (PD-L1) expression in dendritic cells, post-treatment
responders against non-responders. The drawing is draw() in
05_I/create_cd274_momac_boxplot.py; only the arguments below differ. The y
axis is labelled once, on the leftmost box.

Output: cd274_dc_boxplot.{svg,pdf,png} in this directory
Run:    python create_cd274_dc_boxplot.py
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import DC_CELLS_H5AD                             # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "05_I"))
from create_cd274_momac_boxplot import draw                 # noqa: E402

H5AD = DC_CELLS_H5AD
OUT_DIR = Path(__file__).parent
SUB = 4


def main():
    draw(H5AD, SUB, 'DC cells', '', -0.34, 'cd274_dc_boxplot', OUT_DIR)


if __name__ == '__main__':
    main()
