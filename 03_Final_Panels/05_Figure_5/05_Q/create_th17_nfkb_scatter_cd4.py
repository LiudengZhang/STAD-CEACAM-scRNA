#!/usr/bin/env python3
"""Figure 5M: Hallmark TNF-alpha/NF-kB score against a 24-gene Th17 score in
CD4+ T cells, one point per stomach sample. The drawing and statistics are
draw() in 05_O/create_tex_nfkb_scatter_cd8.py; only the arguments differ.

Input:  TCD4_H5AD (00_Config/paths.py)
Output: th17_nfkb_scatter_cd4.{svg,pdf,png} in this directory
Run:    python create_th17_nfkb_scatter_cd4.py
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD4_H5AD                                  # noqa: E402
from cnsfig import Panel                                     # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "05_O"))
from create_tex_nfkb_scatter_cd8 import draw                 # noqa: E402

H5AD = TCD4_H5AD
PANEL = Panel(5, "M", __file__, small_pt=5.0)

# 24-gene Th17 signature.
STATE_GENES = [
    'IL17A','IL17F','RORC','CCR6','IL23R','IL22','AHR','BATF','IRF4','STAT3',
    'CCL20','CXCR6','KLRB1','IL21','IL1R1','RORA','CTSH','PTPN13','TMEM176A',
    'TMEM176B','CAPG','LGMN','FKBP5','ICOS',
]


def main():
    draw(PANEL, H5AD, STATE_GENES, 'Th17 (CD4+)', 'Th17 Score',
         'th17_nfkb_scatter_cd4')


if __name__ == '__main__':
    main()
