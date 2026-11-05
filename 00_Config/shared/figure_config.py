#!/usr/bin/env python3
"""
Shared matplotlib settings and colour palettes.

use_panel_style() sets rcParams for scripts that draw at 4x the printed size.
Panels drawn at 1:1 use panel_style_cns instead. The palettes
(RESPONSE_COLORS, GROUP_COLORS, CELL_TYPE_COLORS) are used by both.

Usage:
    from shared.figure_config import use_panel_style
"""

import matplotlib.pyplot as plt

DPI = 300
MM_TO_INCH = 1 / 25.4

# Double-column width and maximum height (mm)
NC_DOUBLE_COL_MM = 183
NC_MAX_HEIGHT_MM = 170


# Scripts using use_panel_style() draw at four times the printed size.

SCALE = 4
CM_TO_INCH = 1 / 2.54
FONT_STACK = ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']

# Okabe-Ito vermillion and blue (colour-blind safe): responders red,
# non-responders blue. The darker pair is for median lines over the box fill.
RESPONSE_COLORS = {'Responsed': '#D55E00', 'No-response': '#0072B2'}
RESPONSE_MEDIAN_COLORS = {'Responsed': '#A34700', 'No-response': '#005689'}

# Treatment-response group colours: responders red, non-responders blue;
# pre-treatment light, post-treatment saturated. A two-group panel takes the
# pair for its phase.
GROUP_COLORS = {
    'Pre-R':   '#e17055',
    'Pre-NR':  '#74b9ff',
    'Post-R':  '#d63031',
    'Post-NR': '#0984e3',
}
PRE_PAIR = (GROUP_COLORS['Pre-R'], GROUP_COLORS['Pre-NR'])

# Major cell-type palette (ColorBrewer Set2 + Set3).
CELL_TYPE_COLORS = {
    'Epithelial cells': '#e5c494',
    'T/NK cells': '#66c2a5',
    'Monocytes/Macrophages': '#fc8d62',
    'Plasma cells': '#8da0cb',
    'B cells': '#ffd92f',
    'Endothelial cells': '#a6d854',
    'Fibroblasts': '#e78ac3',
    'Neutrophils': '#b3b3b3',
    'Mast cells': '#8c564b',
    'Pericytes': '#bebada',
    'Dendritic cells': '#80b1d3',
    'Hepatocytes': '#bc80bd',
}
POST_PAIR = (GROUP_COLORS['Post-R'], GROUP_COLORS['Post-NR'])


def use_panel_style(font_pt=None, scale=SCALE, **extra):
    """rcParams for a panel drawn at `scale` times its printed size."""
    params = {
        'font.family': 'sans-serif',
        'font.sans-serif': FONT_STACK,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    }
    if font_pt is not None:
        params['font.size'] = font_pt * scale
    params.update(extra)
    plt.rcParams.update(params)
