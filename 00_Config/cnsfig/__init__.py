"""
Shared figure helpers, imported as one package.

    style           panel_style_cns: type sizes, colours, save_panel
    slots           printed panel rectangles in millimetres
    spatial_scale   Visium pixels to microns
    legend          legend and dot-size key layout
    layout          axes placed in millimetres
    rich            labels with italic gene symbols
    boxes           box plots, brackets and their x-limits
    corr_stats      correlation statistics written beside a scatter panel
    cache           a panel's plotted numbers cached in data/*.csv
    panel           the canvas, slot and save for one panel
"""

import sys
from pathlib import Path

_CONFIG = Path(__file__).resolve().parent.parent
if str(_CONFIG) not in sys.path:
    sys.path.insert(0, str(_CONFIG))

import panel_style_cns as style        # noqa: E402,F401
import slots                           # noqa: E402,F401
import spatial_scale                   # noqa: E402,F401

from . import legend                   # noqa: E402,F401
from . import layout                   # noqa: E402,F401
from . import rich                     # noqa: E402,F401
from . import boxes                    # noqa: E402,F401
from . import corr_stats               # noqa: E402,F401
from . import cache                    # noqa: E402,F401
from . import panel                    # noqa: E402,F401
from .legend import (dot_size_key, require_size_key, circle_layout,  # noqa: E402,F401
                     scanpy_dot_areas, axis_size_pt, compact_key,
                     scanpy_compact_key, group_key)
from .layout import umap_grid_mm, scatter_pair_mm, corr_annotate, pin_frame_mm  # noqa: E402,F401
from .rich import rich_text, rich_title, rich_xlabel, rich_ylabel  # noqa: E402,F401
from .boxes import box_xlim, finish_two_group, frame_fixed, title_band_mm  # noqa: E402,F401
from .panel import Panel            # noqa: E402,F401

__all__ = ["style", "slots", "spatial_scale", "legend", "layout", "cache",
           "panel", "Panel",
           "group_key",
           "dot_size_key", "require_size_key", "circle_layout",
           "scanpy_dot_areas", "axis_size_pt", "compact_key",
           "scanpy_compact_key", "umap_grid_mm", "scatter_pair_mm",
           "corr_annotate", "pin_frame_mm", "corr_stats", "rich", "boxes",
           "rich_text", "rich_title", "rich_xlabel", "rich_ylabel",
           "box_xlim", "finish_two_group", "frame_fixed", "title_band_mm"]
