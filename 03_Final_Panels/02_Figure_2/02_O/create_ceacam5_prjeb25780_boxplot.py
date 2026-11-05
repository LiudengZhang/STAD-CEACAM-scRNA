#!/usr/bin/env python3
"""
Figure 2L, right violin: BayesPrism epithelial-deconvolved CEACAM5 expression
in the PRJEB25780 bulk cohort, responders versus non-responders, drawn by
draw() in 02_N/create_ceacam6_prjeb25780_boxplot.py (one-tailed Mann-Whitney
P on the panel, two-sided in the legend). No y label, so the left margin is
narrower.

Reads TIGER_BAYESPRISM_EPI and TIGER_META through that shared drawing.
Writes ceacam5_prjeb25780_boxplot.[svg|pdf|png] beside this script.

Run: python create_ceacam5_prjeb25780_boxplot.py
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "02_N"))
from create_ceacam6_prjeb25780_boxplot import draw        # noqa: E402


def main():
    draw(Path(__file__).parent, 2, 'CEACAM5', 'ceacam5_prjeb25780_boxplot',
         None, 5.4)


if __name__ == '__main__':
    main()

