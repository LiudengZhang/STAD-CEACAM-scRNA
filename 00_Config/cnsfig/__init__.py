"""
One import for the figure code, over modules that already exist.

Nothing moved to build this. `panel_style_cns` has 145 live importers and
`paths` has 204; relocating either would be a diff across the whole tree with
no gain. So this package re-exports what is already the authority and owns only
what nothing owned before:

    style           00_Config/panel_style_cns.py   type sizes, save_panel
    slots           00_Config/slots.py             the printed millimetre rects
    spatial_scale   00_Config/spatial_scale.py     pixels to microns
    legend          NEW - the only place a legend is laid out
    layout          NEW (2026-09-14) - axes placed at millimetres, so two
                    panels asking for the same grid print the same grid
    rich            NEW (2026-09-14) - a label with a gene symbol in it, the
                    symbol italic digits and all, the rest upright
    boxes           NEW (2026-09-14) - the two-group box frame, and the one
                    x-limit rule that keeps boxes off the spines
    cache           NEW (2026-09-15) - the seam between computing and drawing:
                    a panel's numbers live in data/*.csv and the drawing
                    reads them, so redrawing never recomputes

Read `legend.py`'s docstring for why the new one had to exist: across 58 panel
scripts there was no shared legend code at all, and two panels shipped wrong
because of it.
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
from .legend import (dot_size_key, require_size_key, circle_layout,  # noqa: E402,F401
                     scanpy_dot_areas, axis_size_pt, compact_key,
                     scanpy_compact_key, group_key)
from .layout import umap_grid_mm, scatter_pair_mm, corr_annotate, pin_frame_mm  # noqa: E402,F401
from .rich import rich_text, rich_title, rich_xlabel, rich_ylabel  # noqa: E402,F401
from .boxes import box_xlim, finish_two_group, frame_fixed, title_band_mm  # noqa: E402,F401

__all__ = ["style", "slots", "spatial_scale", "legend", "layout", "cache",
           "group_key",
           "dot_size_key", "require_size_key", "circle_layout",
           "scanpy_dot_areas", "axis_size_pt", "compact_key",
           "scanpy_compact_key", "umap_grid_mm", "scatter_pair_mm",
           "corr_annotate", "pin_frame_mm", "corr_stats", "rich", "boxes",
           "rich_text", "rich_title", "rich_xlabel", "rich_ylabel",
           "box_xlim", "finish_two_group", "frame_fixed", "title_band_mm"]
