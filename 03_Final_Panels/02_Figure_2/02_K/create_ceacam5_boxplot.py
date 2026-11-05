#!/usr/bin/env python3
"""
Figure 2K, right box: sample-mean CEACAM5 expression in pre-treatment stomach
epithelium, responders versus non-responders, drawn by draw() in
02_J/create_ceacam6_boxplot.py (one-tailed Mann-Whitney P on the panel,
two-sided in the legend). No y label, so the left margin is narrower.

Reads the epithelial h5ad (EPITHELIAL_H5AD) through that shared drawing.
Writes ceacam5_pre_boxplot.[svg|pdf|png] beside this script.

Run: python create_ceacam5_boxplot.py
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "02_J"))
from create_ceacam6_boxplot import draw                  # noqa: E402


def main():
    draw(Path(__file__).parent, 2, 'CEACAM5', 'ceacam5_pre_boxplot',
         None, 5.4)


if __name__ == '__main__':
    main()
