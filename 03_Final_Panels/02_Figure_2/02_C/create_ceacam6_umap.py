#!/usr/bin/env python3
"""
Figure 2C: CEACAM6 expression on the epithelial UMAP, drawn by draw() in
02_B/create_ceacam5_umap.py.

Reads the epithelial h5ad through that shared drawing.
Writes ceacam6_expression_umap.[svg|pdf|png] beside this script.

Run: python create_ceacam6_umap.py
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "02_B"))
from cnsfig import Panel                                   # noqa: E402
from create_ceacam5_umap import draw, SMALL_PT              # noqa: E402

PANEL = Panel(2, "C", __file__, small_pt=SMALL_PT)


def main():
    draw(PANEL, 'CEACAM6', "ceacam6_expression_umap")


if __name__ == '__main__':
    main()
