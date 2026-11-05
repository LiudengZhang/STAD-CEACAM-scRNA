"""Printed panel rectangles and panel-letter sizes, in millimetres.

Panels are drawn at the exact rectangle they print in, leaving the panel
letter's corner free. The geometry is read from tables in 03_Final_Panels:

    panel_rects.csv     one rect per printed panel, mm from the page corner
    slot_subrects.csv   sub-boxes of panels that print several plots under
                        one letter
    letter_spec.csv     each letter's size, weight and box, relative to its
                        panel's rect

Where a *_v2.csv table exists for a figure it replaces that figure's rows
entirely (a partial override raises).

    import slots
    W, H = slots.size_mm(5, "C")
    fig, ax = style.subplots_mm(W, H)
    ...
    style.fit_margins(fig, cell_mm=slots.letter_cell_mm(5, "C"))

    W, H = slots.size_mm(5, "J", sub=1)      # one box of a grouped panel

Every lookup raises on a miss.
"""

from __future__ import annotations

import csv
from pathlib import Path

PANELS = Path(__file__).resolve().parent.parent / "03_Final_Panels"
RECTS = PANELS / "panel_rects.csv"
SUBRECTS = PANELS / "slot_subrects.csv"
RECTS_V2 = PANELS / "panel_rects_v2.csv"
SUBRECTS_V2 = PANELS / "slot_subrects_v2.csv"
LETTERS = PANELS / "letter_spec.csv"
# One letter size for every panel.
LETTERS_V2 = PANELS / "letter_spec_v2.csv"
# Printed page size of each figure; the v2 table overrides where present.
PAGES = PANELS / "page_size.csv"
PAGES_V2 = PANELS / "page_size_v2.csv"

__all__ = ["rect_mm", "size_mm", "letter_cell_mm", "letter_spec",
           "n_subrects", "page_mm"]

_cache: dict[Path, list[dict]] = {}


def _rows(path: Path) -> list[dict]:
    if path not in _cache:
        if not path.exists():
            raise FileNotFoundError(
                f"{path.name} is missing. It is written "
                f"by 03_Final_Panels/build_slot_geometry.py.")
        with path.open() as fh:
            _cache[path] = list(csv.DictReader(fh))
    return _cache[path]


def _v2(base: Path, over: Path) -> list[dict]:
    """The measured rows, with a figure's rows replaced wholesale where an
    override exists."""
    rows = _rows(base)
    if not over.exists():
        return rows
    new = _rows(over)
    by_fig: dict[str, list[dict]] = {}
    for r in new:
        by_fig.setdefault(r["figure"], []).append(r)
    for fig, rs in by_fig.items():
        # Tables without a printed_panel column (one row per figure) always pass.
        want = {r.get("printed_panel", "") for r in rows if r["figure"] == fig}
        have = {r.get("printed_panel", "") for r in rs}
        missing = want - have
        if missing:
            raise ValueError(
                f"{over.name} covers {fig} only partly: {sorted(missing)} are "
                f"still on the base grid. A figure is all-or-nothing; "
                f"add them or drop {fig} from the override.")
    kept = [r for r in rows if r["figure"] not in by_fig]
    return kept + new


def _figure_key(figure) -> str:
    """Accept 5, "5" or "Figure 5"; the tables spell it out."""
    text = str(figure)
    return text if text.startswith("Figure") else f"Figure {text}"


def _parse(rect: str) -> tuple[float, float, float, float]:
    x0, y0, x1, y1 = (float(v) for v in rect.split(","))
    return x0, y0, x1, y1


def rect_mm(figure, panel, sub=None):
    """(x0, y0, x1, y1) in mm from the page's top-left corner.

    `sub` selects one box of a panel that prints as a group; leave it None for
    the whole panel. Sub-boxes are numbered from 1 in printed reading order.
    """
    fig = _figure_key(figure)
    if sub is None:
        for r in _v2(RECTS, RECTS_V2):
            if r["figure"] == fig and r["printed_panel"] == panel:
                return _parse(r["rect_mm"])
        raise KeyError(f"no rect for {fig} panel {panel}")
    for r in _v2(SUBRECTS, SUBRECTS_V2):
        if (r["figure"] == fig and r["printed_panel"] == panel
                and int(r["sub_index"]) == sub):
            return _parse(r["rect_mm"])
    raise KeyError(f"no sub-box {sub} for {fig} panel {panel}; it has "
                   f"{n_subrects(figure, panel)}")


def n_subrects(figure, panel) -> int:
    """How many boxes this panel prints under its one letter. 0 if it is one."""
    fig = _figure_key(figure)
    return sum(1 for r in _v2(SUBRECTS, SUBRECTS_V2)
               if r["figure"] == fig and r["printed_panel"] == panel)


def size_mm(figure, panel, sub=None) -> tuple[float, float]:
    """(width, height) in mm of the rectangle this panel prints in."""
    x0, y0, x1, y1 = rect_mm(figure, panel, sub)
    return round(x1 - x0, 3), round(y1 - y0, 3)


def letter_spec(figure, panel) -> dict:
    """The panel letter's size, weight and position, relative to its rect."""
    fig = _figure_key(figure)
    for r in _v2(LETTERS, LETTERS_V2):
        if r["figure"] == fig and r["printed_panel"] == panel:
            return {"size_pt": float(r["size_pt"]),
                    "weight": r["weight"],
                    "font": r["font"],
                    "x_mm": float(r["x_mm"]),
                    "y_mm": float(r["y_mm"]),
                    "cell_w_mm": float(r["cell_w_mm"]),
                    "cell_h_mm": float(r["cell_h_mm"])}
    raise KeyError(f"no letter spec for {fig} panel {panel}")


def letter_cell_mm(figure, panel) -> tuple[float, float]:
    """The keep-out box, in mm, measured from the panel's top-left corner.

    Pass it to `panel_style_cns.fit_margins(cell_mm=...)`. The box includes
    the indent between the rect edge and the glyph.
    """
    spec = letter_spec(figure, panel)
    return spec["cell_w_mm"], spec["cell_h_mm"]


def page_mm(figure) -> tuple[float, float]:
    """(width, height) of the printed page, in mm."""
    fig = _figure_key(figure)
    for r in _v2(PAGES, PAGES_V2):
        if r["figure"] == fig:
            return float(r["width_mm"]), float(r["height_mm"])
    raise KeyError(f"no page size for {fig}")
