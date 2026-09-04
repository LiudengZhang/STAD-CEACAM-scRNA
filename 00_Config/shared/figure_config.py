#!/usr/bin/env python3
"""
Nature Cancer figure configuration.

Sets matplotlib rcParams for publication-quality figures.
All panel scripts should call apply_nature_style() before plotting.

Usage:
    from shared.figure_config import apply_nature_style, create_figure, save_panel
"""

import matplotlib.pyplot as plt
import matplotlib
import os

# =============================================================================
# Constants
# =============================================================================

DPI = 300
MM_TO_INCH = 1 / 25.4

# Nature Cancer column widths (mm)
NC_SINGLE_COL_MM = 88
NC_ONE_HALF_COL_MM = 120
NC_DOUBLE_COL_MM = 183
NC_MAX_HEIGHT_MM = 170

# In inches (for matplotlib)
NC_SINGLE_COL = NC_SINGLE_COL_MM * MM_TO_INCH
NC_DOUBLE_COL = NC_DOUBLE_COL_MM * MM_TO_INCH
NC_MAX_HEIGHT = NC_MAX_HEIGHT_MM * MM_TO_INCH


def apply_nature_style():
    """Apply Nature Cancer rcParams globally."""
    plt.rcParams.update({
        # Font
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6,
        'axes.labelsize': 6,
        'axes.titlesize': 7,
        'xtick.labelsize': 5,
        'ytick.labelsize': 5,
        'legend.fontsize': 5,

        # Line weights
        'axes.linewidth': 0.5,
        'lines.linewidth': 0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 2,
        'ytick.major.size': 2,

        # Font embedding — keep text as editable <text> in SVG, Type 42 in PDF
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,

        # Resolution
        'savefig.dpi': DPI,
        'figure.dpi': DPI,

        # Clean defaults
        'axes.spines.top': False,
        'axes.spines.right': False,
        'legend.frameon': False,
    })


def create_figure(width_mm, height_mm):
    """
    Create a matplotlib figure at exact target size.

    Args:
        width_mm: Panel width in millimeters (from layout engine)
        height_mm: Panel height in millimeters (from layout engine)

    Returns:
        (fig, ax) tuple
    """
    fig, ax = plt.subplots(
        figsize=(width_mm * MM_TO_INCH, height_mm * MM_TO_INCH)
    )
    return fig, ax


def save_panel(fig, output_path, formats=('svg', 'png')):
    """
    Save panel in specified formats.

    Args:
        fig: matplotlib Figure
        output_path: Path without extension
        formats: Tuple of formats to save ('svg', 'png', 'pdf')
    """
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)

    for fmt in formats:
        fig.savefig(
            f"{output_path}.{fmt}",
            dpi=DPI,
            bbox_inches='tight',
            facecolor='white',
            pad_inches=0.01,
            transparent=False,
        )
        print(f"Saved: {output_path}.{fmt}")

    plt.close(fig)


# =============================================================================
# Label helpers
# =============================================================================

import re

def strip_cx(label):
    """Remove Cx_ prefix from cell state labels for display.

    'C0_Epi_PTMA' → 'Epi_PTMA'
    'C2_CEACAM5/6' → 'CEACAM5/6'
    'C6_Tex' → 'Tex'
    """
    return re.sub(r'^C\d+_', '', str(label))


def format_pval_stars(p):
    """Convert p-value to 3-tier star notation.

    *** : P < 0.001
    **  : P < 0.01
    *   : P < 0.05
    ns  : P >= 0.05
    """
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'


# Immune module display names
IMMUNE_MODULE_NAMES = {
    1: 'IM-T/NK/DC',
    2: 'IM-MoMac',
    3: 'IM-Mixed',
    4: 'IM-Neutrophil',
    5: 'IM-B/Plasma',
}


# =============================================================================
# The 4x scaling convention
# =============================================================================
#
# apply_nature_style() above draws at print size. Fifty-nine of the panel
# scripts do the opposite: they draw at four times the printed size and let the
# journal scale the figure down, which keeps text crisp and line weights true.
# Both conventions are in use and neither is wrong, so the 4x helpers live
# beside the print-size ones rather than replacing them.
#
# The scripts are not identical in what they set, and a shared helper written as
# if they were would restyle a third of the figures. Measured across the tree:
#
#     font.sans-serif   identical everywhere
#     svg/pdf/ps type   identical everywhere
#     font.size         7 x SCALE in 22 scripts, 8 x in 6, 6 x in 4, 5 x in 3,
#                       9 x in 3, absent in 7, unscaled 6 pt in 4
#     SCALE             4 in 50 scripts, undefined in 9
#     panel size        fifteen different values
#
# So only the genuinely common part moves here - the font stack and the three
# font-type settings, which are what keep SVG text editable and PDF fonts
# embedded. Everything that varies stays an argument.

SCALE = 4
CM_TO_INCH = 1 / 2.54
FONT_STACK = ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']

# Okabe-Ito blue and vermillion, so the figures stay readable in the commonest
# forms of colour blindness. The darker pair is for median lines drawn over the
# box fill.
RESPONSE_COLORS = {'Responsed': '#0072B2', 'No-response': '#D55E00'}
RESPONSE_MEDIAN_COLORS = {'Responsed': '#005689', 'No-response': '#A34700'}


def use_panel_style(font_pt=None, scale=SCALE, **extra):
    """rcParams for a panel drawn at `scale` times its printed size.

    font_pt is the printed point size and is multiplied by scale. Pass None to
    leave the font size alone, which is what the seven scripts that never set
    one expect. Anything else - axes.linewidth, axes.labelsize, a size that is
    deliberately not scaled - goes through **extra and wins.
    """
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


def panel_figsize(width_cm, height_cm, scale=SCALE):
    """Figure size in inches for a panel of the given printed size in cm."""
    return (width_cm * scale * CM_TO_INCH, height_cm * scale * CM_TO_INCH)


def save_panel_scaled(fig, stem, dpi=DPI, formats=('png', 'svg', 'pdf')):
    """Write a 4x panel in the formats the assembly and the journal need.

    Unlike save_panel() above this passes no bbox_inches: the 4x scripts set
    their own margins with subplots_adjust, and cropping to the ink would undo
    that and change every panel's proportions.
    """
    for ext in formats:
        fig.savefig(f'{stem}.{ext}', dpi=dpi, facecolor='white')
    return f'{stem}.png'
