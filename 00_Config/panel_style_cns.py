#!/usr/bin/env python3
"""
Figure style for all panels, built on cnsplots
(https://github.com/faridrashidi/cnsplots): a Helvetica-metric font stack,
millimetre geometry (a panel is drawn at the size it prints at), an exact-size
save, P-value labels and the two line weights RULE_PT and EDGE_PT.

    import panel_style_cns as style
    style.apply()
    fig, ax = style.subplots_mm(100, 26)
    ...
    style.save_panel(fig, out_dir / "panel_name")

Run directly for a self-test; --dump-snapshot writes panel_style_rc.json for
environments without cnsplots.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

import json
import types

try:
    import cnsplots as cns
    FROM_SNAPSHOT = False
except ImportError:
    cns = None
    FROM_SNAPSHOT = True

# Where the settings snapshot is written and read.
SNAPSHOT = Path(__file__).resolve().parent / "panel_style_rc.json"

def dump_snapshot(path=None, **setting_overrides) -> Path:
    """Write cnsplots' settings and the resulting rcParams to a JSON file,
    for drawing in an environment without cnsplots."""
    if FROM_SNAPSHOT:
        raise RuntimeError("dump_snapshot needs cnsplots; run it in the "
                           "environment where cnsplots is installed")
    apply(**setting_overrides)
    settings = {k: getattr(cns.settings, k) for k in dir(cns.settings)
                if not k.startswith("_")
                and _jsonable(getattr(cns.settings, k))}
    snap = {
        "cnsplots_version": cns.__version__,
        "overrides": {k: v for k, v in setting_overrides.items()},
        "resolved_family": resolved_family(),
        "settings": settings,
        "rcparams": {k: _jsonable_rc(v)
                     for k, v in matplotlib.rcParams.items()
                     if _jsonable(v)},
    }
    path = Path(path) if path else SNAPSHOT
    path.write_text(json.dumps(snap, indent=1, sort_keys=True))
    return path


def _jsonable(v):
    return isinstance(v, (str, int, float, bool, type(None))) or (
        isinstance(v, (list, tuple))
        and all(isinstance(x, (str, int, float, bool, type(None))) for x in v))


def _jsonable_rc(v):
    return list(v) if isinstance(v, tuple) else v


class _SnapshotCnsplots:
    """Stands in for cnsplots from a recorded snapshot; raises if the
    requested settings differ from the recorded ones."""

    # font_sans_serif depends on installed fonts; compared via resolved_family().
    _NOT_COMPARED = ("font_sans_serif",)

    def __init__(self, path):
        if not Path(path).exists():
            raise RuntimeError(
                f"cnsplots is not importable here and no settings snapshot is "
                f"at {path}. Write one from the environment that has cnsplots: "
                f"python panel_style_cns.py --dump-snapshot")
        self._snap = json.loads(Path(path).read_text())
        self.__version__ = self._snap["cnsplots_version"] + " (snapshot)"
        self.settings = types.SimpleNamespace(**self._snap["settings"])
        self._recorded = dict(self._snap["settings"])

    def setup_matplotlib(self):
        differs = [k for k, v in self._recorded.items()
                   if k not in self._NOT_COMPARED
                   and getattr(self.settings, k, None) != v]
        if differs:
            raise RuntimeError(
                "the settings snapshot was taken with different values for "
                + ", ".join(sorted(differs))
                + "; re-take it with these overrides in the environment that "
                  "has cnsplots")
        matplotlib.rcParams.update(self._snap["rcparams"])
        got = resolved_family()
        want = self._snap["resolved_family"]
        if got != want:
            raise RuntimeError(f"the snapshot was taken with {want} but this "
                               f"machine resolves the stack to {got}")


if cns is None:
    cns = _SnapshotCnsplots(SNAPSHOT)


__all__ = [
    "apply", "describe", "resolved_family",
    "PAGE_W_MM", "PAGE_W_PT", "MM_TO_INCH",
    "figsize_mm", "figure_mm", "subplots_mm", "margins_mm", "overflow_mm",
    "fit_margins", "letter_clear", "LETTER_CELL_MM", "dump_snapshot",
    "save_panel", "body_pt", "tick_pt", "letter_pt",
    "RULE_PT", "EDGE_PT", "p_label", "p_is_star", "STAR_PT", "p_text_kw",
]

# Page geometry. cnsplots' multipanel_max_width is in points (540 pt = 190.5 mm).
MM_TO_INCH = 1.0 / 25.4
PT_PER_MM = 72.0 / 25.4

# Line weights.
#   RULE_PT  lines: spines, ticks, box/whisker/median, brackets, rules,
#            colour-bar outlines, error bars (cnsplots' tick width).
#   EDGE_PT  outlines of filled things: marker edges, heatmap cell grids,
#            hexagon outlines, violin bodies.
RULE_PT = 0.6
EDGE_PT = 0.25

# rcParams cnsplots leaves at matplotlib's defaults; apply() sets them to RULE_PT.
_RULE_RC = (
    "axes.linewidth", "lines.linewidth", "patch.linewidth", "grid.linewidth",
    "hatch.linewidth",
    "xtick.major.width", "ytick.major.width",
    "xtick.minor.width", "ytick.minor.width",
    "boxplot.boxprops.linewidth", "boxplot.whiskerprops.linewidth",
    "boxplot.capprops.linewidth", "boxplot.medianprops.linewidth",
    "boxplot.meanprops.linewidth", "boxplot.flierprops.markeredgewidth",
)


def p_exact(p: float) -> str:
    """The exact form of a P value, never a star.

    "P < 0.001" below 0.001; three decimals below 0.01; otherwise two
    decimals, or three when two would print "0.05". At and above 0.05 this
    equals p_label's string.
    """
    p = float(p)
    if p != p or p < 0 or p > 1:
        raise ValueError(f"not a P value: {p!r}")
    if p < 0.001:
        return "P < 0.001"
    if p < 0.01:
        return f"P = {p:.3f}"
    two = f"{p:.2f}"
    if two == "0.05":
        return f"P = {p:.3f}"
    return f"P = {two}"


def p_label(p: float, exact: bool = False) -> str:
    """The text a panel prints for a P value.

    Below 0.05 a star: "*" (< 0.05), "**" (< 0.01), "***" (< 0.001). At or
    above 0.05, two decimals, or three when two would print "0.05" (so 0.052
    prints as "P = 0.052"). exact=True prints p_exact() instead, never a star.
    """
    if exact:
        return p_exact(p)
    p = float(p)
    if p != p or p < 0 or p > 1:
        raise ValueError(f"not a P value: {p!r}")
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    two = f"{p:.2f}"
    if two == "0.05":
        return f"P = {p:.3f}"
    return f"P = {two}"


def p_is_star(p: float) -> bool:
    """Whether p_label prints a star for this value."""
    return p_label(p).startswith("*")


# Star size, pt; larger than the digits because an asterisk's ink is small.
STAR_PT = 9.0


def p_text_kw(p: float, exact: bool = False) -> tuple[str, dict]:
    """(the label p_label prints, the text kwargs that size it).

    A star is set at STAR_PT, a printed value at tick_pt():
        s, kw = style.p_text_kw(p); ax.text(x, y, s, **kw)
    """
    label = p_label(p, exact=exact)
    return label, {"fontsize": STAR_PT if label.startswith("*") else tick_pt()}


# On-plot cluster labels of the categorical UMAPs (Figures 1B, 2A, 3C, 4D).
def umap_label_kw(**overrides) -> dict:
    """Text kwargs for an on-plot UMAP cluster label:

        ax.text(cx, cy, name, **style.umap_label_kw())

    `overrides` is for placement (transform, zorder), not size, weight or
    colour.
    """
    kw = {"fontsize": body_pt(), "fontweight": "normal", "color": "black",
          "ha": "center", "va": "center", "zorder": 10, "gid": "umap_label"}
    kw.update(overrides)
    return kw


def umap_label_underlays(ax, stroke_pt=1.5):
    """White-stroked copies under each UMAP label of `ax`, so the black
    label on top stays plain text in the PDF. Returns the copies."""
    import matplotlib.patheffects as pe
    out = []
    for t in [t for t in ax.texts if t.get_gid() == "umap_label"]:
        u = ax.text(*t.get_position(), t.get_text(), transform=t.get_transform(),
                    ha=t.get_ha(), va=t.get_va(), fontsize=t.get_fontsize(),
                    fontweight=t.get_fontweight(), color="white",
                    zorder=t.get_zorder() - 0.5, gid="umap_label_underlay")
        u.set_path_effects([pe.withStroke(linewidth=stroke_pt,
                                          foreground="white")])
        out.append(u)
    return out


def _page_width_pt() -> float:
    return float(cns.settings.multipanel_max_width)


PAGE_W_PT = _page_width_pt()
PAGE_W_MM = PAGE_W_PT / PT_PER_MM          # 190.5 mm

# Font stack: cnsplots' order with the Helvetica-metric clones ahead of DejaVu
# Sans. DejaVu Sans stays last for its glyph coverage (mathtext, symbols).
_LOCAL_HELVETICA_CLONES = ("Nimbus Sans", "Liberation Sans")


def _extended_font_stack() -> tuple[str, ...]:
    base = tuple(cns.CNSSettings._defaults["font_sans_serif"]) \
        if hasattr(cns, "CNSSettings") else ("Helvetica", "Helvetica Neue",
                                             "Arial", "DejaVu Sans")
    head = [f for f in base if f != "DejaVu Sans"]
    tail = [f for f in base if f == "DejaVu Sans"] or ["DejaVu Sans"]
    for clone in _LOCAL_HELVETICA_CLONES:
        if clone not in head:
            head.append(clone)
    return tuple(head + tail)


def resolved_family(stack: tuple[str, ...] | None = None) -> str:
    """The first family in the stack that this machine can actually draw."""
    stack = stack or tuple(cns.settings.font_sans_serif)
    installed = {f.name for f in fm.fontManager.ttflist}
    for name in stack:
        if name in installed:
            return name
    return "DejaVu Sans"


# Apply

def _stack_led_by_resolved() -> tuple[str, ...]:
    """The font stack with the resolved family first, so other SVG
    renderers pick the face matplotlib drew."""
    stack = list(_extended_font_stack())
    first = resolved_family(tuple(stack))
    if first in stack:
        stack.remove(first)
    return tuple([first] + stack)


def apply(palette: bool = False, **setting_overrides) -> str:
    """Put cnsplots' type system into matplotlib's rcParams.

    Parameters
    ----------
    palette
        False (default) keeps the caller's existing `axes.prop_cycle` and
        `image.cmap`. True lets cnsplots' palettes through.
    **setting_overrides
        Written onto `cns.settings` before `setup_matplotlib()` runs, and
        validated by cnsplots. An unknown name raises.

    Returns
    -------
    str
        The font family that will actually be drawn on this machine.
    """
    cns.settings.font_sans_serif = _extended_font_stack()
    cns.settings.font_sans_serif = _stack_led_by_resolved()

    # Save the full canvas; cnsplots' bbox='tight' would crop it and change the
    # printed size.
    cns.settings.savefig_bbox = "standard"

    cns.settings.savefig_transparent = False

    # Regular-weight titles (bold only for panel letters); spines and ticks at
    # RULE_PT.
    setting_overrides.setdefault("title_fontweight", "normal")
    setting_overrides.setdefault("axes_linewidth", RULE_PT)
    setting_overrides.setdefault("xtick_major_width", RULE_PT)
    setting_overrides.setdefault("ytick_major_width", RULE_PT)

    for key, value in setting_overrides.items():
        setattr(cns.settings, key, value)

    # Legend text at the tick-label size unless the caller set legend_fontsize.
    if "legend_fontsize" not in setting_overrides:
        cns.settings.legend_fontsize = cns.settings.fontsize_legend

    keep_cycle = matplotlib.rcParams["axes.prop_cycle"]
    keep_cmap = matplotlib.rcParams["image.cmap"]

    cns.setup_matplotlib()

    if not palette:
        matplotlib.rcParams["axes.prop_cycle"] = keep_cycle
        matplotlib.rcParams["image.cmap"] = keep_cmap

    # TrueType fonts in EPS output too.
    matplotlib.rcParams["ps.fonttype"] = 42

    # Set after setup_matplotlib() so the snapshot path ends in the same state.
    for key in _RULE_RC:
        matplotlib.rcParams[key] = RULE_PT
    # A line marker's outline is an edge (matplotlib's default is 1.0 pt).
    matplotlib.rcParams["lines.markeredgewidth"] = EDGE_PT

    return resolved_family()


def setup_scanpy() -> None:
    """Configure scanpy through cnsplots. Raises without cnsplots: the
    snapshot records rcParams, not scanpy's settings.
    """
    # Test FROM_SNAPSHOT: `cns` is never None here (the snapshot stand-in).
    if FROM_SNAPSHOT or not hasattr(cns, "setup_scanpy"):
        raise RuntimeError(
            "setup_scanpy() needs cnsplots, which is not importable here. "
            "The style snapshot (panel_style_rc.json) carries rcParams and not "
            "scanpy's own settings, so a scanpy-built figure cannot be drawn to "
            "print size from it. Run this panel in the environment where "
            "cnsplots is installed (see update_release.fix_environment), or use "
            "the approved asset that install_dotplot() copies.")
    cns.setup_scanpy()


def body_pt() -> float:
    """Axis titles and axis labels. cnsplots `title_fontsize`."""
    return float(cns.settings.title_fontsize)


def tick_pt() -> float:
    """Tick labels and legend text. cnsplots `fontsize_legend`."""
    return float(cns.settings.fontsize_legend)


def letter_pt() -> float:
    """Panel letters. cnsplots sets them at `title_fontsize`, bold."""
    return float(cns.settings.title_fontsize)


def letter_font() -> tuple[str, str]:
    """(family, weight) for a panel letter, from cnsplots' settings."""
    return (resolved_family((cns.settings.panel_label_fontname,)
                            + tuple(cns.settings.font_sans_serif)),
            str(cns.settings.panel_label_fontweight))


def describe() -> str:
    """One block of text recording exactly what apply() did. For the report."""
    rc = matplotlib.rcParams
    fam = resolved_family()
    return "\n".join([
        f"cnsplots {cns.__version__}",
        f"  page width          {PAGE_W_PT:.0f} pt = {PAGE_W_MM:.1f} mm",
        f"  font stack          {', '.join(cns.settings.font_sans_serif)}",
        f"  resolves to         {fam}",
        f"  panel letter        {letter_pt():g} pt "
        f"{cns.settings.panel_label_fontweight}",
        f"  axis title/label    {rc['axes.titlesize']:g} / "
        f"{rc['axes.labelsize']:g} pt, weight {rc['axes.titleweight']}",
        f"  tick labels         {rc['xtick.labelsize']:g} / "
        f"{rc['ytick.labelsize']:g} pt",
        f"  legend text         {rc['legend.fontsize']:g} pt, "
        f"frameon={rc['legend.frameon']}, markerscale={rc['legend.markerscale']}",
        f"  axes.linewidth      {rc['axes.linewidth']:g}",
        f"  rule width          {RULE_PT:g} pt (spines, ticks, boxes, "
        f"brackets, rules; lines.linewidth {rc['lines.linewidth']:g}, "
        f"patch {rc['patch.linewidth']:g}, boxplot "
        f"{rc['boxplot.boxprops.linewidth']:g})",
        f"  edge width          {EDGE_PT:g} pt (marker edges, cell grids, "
        f"hex outlines, violin bodies)",
        f"  P values            {p_label(0.057143)}, {p_label(0.054681)}, "
        f"{p_label(0.005)}",
        f"  ticks               size {rc['xtick.major.size']:g}, "
        f"width {rc['xtick.major.width']:g}, pad {rc['xtick.major.pad']:g}",
        f"  spines top/right    {rc['axes.spines.top']} / {rc['axes.spines.right']}",
        f"  svg.fonttype        {rc['svg.fonttype']}",
        f"  pdf.fonttype        {rc['pdf.fonttype']}",
        f"  savefig             dpi {rc['savefig.dpi']:g}, bbox "
        f"{rc['savefig.bbox']}, pad {rc['savefig.pad_inches']:g}, "
        f"transparent {rc['savefig.transparent']}",
        f"  colour cycle        {'cnsplots' if _cycle_is_cns() else 'unchanged (stage 4 pending)'}",
    ])


def _cycle_is_cns() -> bool:
    try:
        first = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"][0]
        return first.lower() != "#1f77b4"
    except Exception:
        return False


# Millimetre geometry

def figsize_mm(w_mm: float, h_mm: float) -> tuple[float, float]:
    """Figure size in inches for a panel that prints w_mm x h_mm. 1:1."""
    return (w_mm * MM_TO_INCH, h_mm * MM_TO_INCH)


def figure_mm(w_mm: float, h_mm: float, **kwargs):
    return plt.figure(figsize=figsize_mm(w_mm, h_mm), **kwargs)


def subplots_mm(w_mm: float, h_mm: float, *args, **kwargs):
    kwargs.setdefault("figsize", figsize_mm(w_mm, h_mm))
    return plt.subplots(*args, **kwargs)


def margins_mm(fig, left=None, right=None, top=None, bottom=None, **kw):
    """subplots_adjust, but the margins are millimetres of paper.

    Any margin left as None keeps its current value.
    """
    w_in, h_in = fig.get_size_inches()
    w_mm, h_mm = w_in / MM_TO_INCH, h_in / MM_TO_INCH
    adj = {}
    if left is not None:
        adj["left"] = left / w_mm
    if right is not None:
        adj["right"] = 1.0 - right / w_mm
    if bottom is not None:
        adj["bottom"] = bottom / h_mm
    if top is not None:
        adj["top"] = 1.0 - top / h_mm
    adj.update(kw)
    fig.subplots_adjust(**adj)
    return adj


def _ink_artists(fig):
    """Every artist that puts marks outside its own axes box."""
    out = list(fig.texts)                       # suptitle, supxlabel, supylabel
    for ax in fig.axes:
        out.append(ax)
        # Left and right titles are separate Text artists from ax.title.
        out.extend([ax.title,
                    getattr(ax, "_left_title", None),
                    getattr(ax, "_right_title", None),
                    ax.xaxis.label, ax.yaxis.label])
        out.extend(ax.texts)
        leg = ax.get_legend()
        if leg is not None:
            out.append(leg)
        for axis, lo_hi in ((ax.xaxis, ax.get_xlim()),
                            (ax.yaxis, ax.get_ylim())):
            lo, hi = sorted(lo_hi)
            for tick, loc in zip(axis.get_major_ticks(),
                                 axis.get_majorticklocs()):
                if lo - 1e-9 <= loc <= hi + 1e-9:
                    out.append(tick.label1)
                    out.append(tick.label2)
    if fig.legends:
        out.extend(fig.legends)
    return [a for a in out if a is not None and getattr(a, "get_visible", bool)()]


def overflow_mm(fig):
    """How far the drawn ink falls outside the canvas, in mm, per side."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi = fig.dpi
    w_in, h_in = fig.get_size_inches()
    left = right = bottom = top = 0.0
    for artist in _ink_artists(fig):
        try:
            bb = artist.get_window_extent(renderer=r)
        except Exception:
            continue
        if bb is None or bb.width <= 0 or bb.height <= 0:
            continue
        left = max(left, -bb.x0 / dpi)
        right = max(right, bb.x1 / dpi - w_in)
        bottom = max(bottom, -bb.y0 / dpi)
        top = max(top, bb.y1 / dpi - h_in)
    return tuple(round(max(0.0, v) * 25.4, 3)
                 for v in (left, right, bottom, top))


def _ink_bbox_in(fig):
    """Union of every ink artist's extent, in inches from the canvas corner."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi = fig.dpi
    x0 = y0 = float("inf")
    x1 = y1 = float("-inf")
    for artist in _ink_artists(fig):
        try:
            bb = artist.get_window_extent(renderer=r)
        except Exception:
            continue
        if bb is None or bb.width <= 0 or bb.height <= 0:
            continue
        x0, y0 = min(x0, bb.x0), min(y0, bb.y0)
        x1, y1 = max(x1, bb.x1), max(y1, bb.y1)
    if x0 == float("inf"):
        raise RuntimeError("the figure draws no measurable ink")
    return x0 / dpi, y0 / dpi, x1 / dpi, y1 / dpi


def fit_margins(fig, pad_mm=0.6, max_iter=24, tol_mm=0.05,
                reserve_letter=True, cell_mm=None):
    """Set the margins from the rendered ink, leaving `pad_mm` of paper.

    Renders the figure, measures the union of all labels, titles, annotations
    and legends, and moves the axes box until the ink sits `pad_mm` inside
    every edge. Works only for axes placed by subplots_adjust (not add_axes).
    Raises if the ink has not converged after `max_iter` passes.

    With `reserve_letter`, the top-left panel-letter cell is also kept clear by
    moving the axes down or right, whichever costs less space.

    Returns the final overflow, which is (0, 0, 0, 0) on success.
    """
    cell_mm = cell_mm or LETTER_CELL_MM
    pad = pad_mm * MM_TO_INCH
    w_in, h_in = fig.get_size_inches()
    dpi = fig.dpi
    cx1 = cell_mm[0] * MM_TO_INCH * dpi
    cy0 = h_in * dpi - cell_mm[1] * MM_TO_INCH * dpi
    move = None                      # "down" or "right", decided once and kept
    tol = tol_mm * MM_TO_INCH
    for _ in range(max_iter):
        fig.canvas.draw()
        x0, y0, x1, y1 = _ink_bbox_in(fig)
        dl, dr = pad - x0, (w_in - pad) - x1
        db, dt = pad - y0, (h_in - pad) - y1
        # Judge the letter corner only once the ink is inside the canvas.
        inside = (dl <= tol and dr >= -tol and db <= tol and dt >= -tol)
        if reserve_letter and (inside or move is not None):
            down, right = _letter_cell_demand(fig, cx1, cy0)
            if down > 0:
                # Vertically centred artists move half as far as the top edge,
                # so the demand is doubled; the slack covers the tolerance.
                slack = 2 * tol_mm * MM_TO_INCH * dpi
                down, right = 2 * down + slack, right + slack
                if move is None:
                    move = "down" if down <= right else "right"
            if move == "down":
                # Never raise the top again once lowered, or the loop oscillates.
                dt = min(dt, -down / dpi if down > 0 else 0.0)
            elif move == "right":
                dl = max(dl, right / dpi if down > 0 else 0.0)
        if max(abs(dl), abs(dr), abs(db), abs(dt)) / MM_TO_INCH < tol_mm:
            break
        sp = fig.subplotpars
        left, right_ = sp.left * w_in + dl, sp.right * w_in + dr
        bottom, top = sp.bottom * h_in + db, sp.top * h_in + dt
        if right_ - left < 0.15 * w_in or top - bottom < 0.15 * h_in:
            raise RuntimeError(
                "the ink cannot be fitted into this canvas: the axes box would "
                "have to shrink below 15% of it. The panel needs shorter "
                "strings or a larger slot.")
        fig.subplots_adjust(left=left / w_in, right=right_ / w_in,
                            bottom=bottom / h_in, top=top / h_in)
    over = overflow_mm(fig)
    if max(over) > tol_mm:
        raise RuntimeError(f"ink still outside the canvas after {max_iter} "
                           f"passes: {over} mm (left, right, bottom, top)")
    if reserve_letter and letter_clear(fig, cell_mm):
        raise RuntimeError("the panel-letter cell is still not clear after "
                           f"{max_iter} passes: {letter_clear(fig, cell_mm)}")
    return over


def _letter_cell_demand(fig, cx1, cy0):
    """How far ink must move down, or right, to leave the letter cell."""
    r = fig.canvas.get_renderer()
    down = right = 0.0
    for artist in _ink_artists(fig):
        try:
            bb = artist.get_window_extent(renderer=r)
        except Exception:
            continue
        if bb is None or bb.width <= 0 or bb.height <= 0:
            continue
        if min(bb.x1, cx1) - bb.x0 > 0 and bb.y1 - max(bb.y0, cy0) > 0:
            down = max(down, bb.y1 - cy0)
            right = max(right, cx1 - bb.x0)
    return down, right


# Cell reserved at a panel's top-left for its letter, mm (w, h).
LETTER_CELL_MM = (3.5, 3.6)


def letter_clear(fig, cell_mm=LETTER_CELL_MM):
    """Which artists intrude into the panel-letter cell, and by how much."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi = fig.dpi
    w_in, h_in = fig.get_size_inches()
    cw, ch = (v * MM_TO_INCH * dpi for v in cell_mm)
    cx0, cx1 = 0.0, cw
    cy0, cy1 = h_in * dpi - ch, h_in * dpi
    hits = []
    for artist in _ink_artists(fig):
        try:
            bb = artist.get_window_extent(renderer=r)
        except Exception:
            continue
        if bb is None or bb.width <= 0 or bb.height <= 0:
            continue
        ox = min(bb.x1, cx1) - max(bb.x0, cx0)
        oy = min(bb.y1, cy1) - max(bb.y0, cy0)
        if ox > 0 and oy > 0:
            area = (ox / dpi / MM_TO_INCH) * (oy / dpi / MM_TO_INCH)
            label = getattr(artist, "get_text", lambda: "")() or type(
                artist).__name__
            hits.append((str(label)[:40], round(area, 3)))
    return sorted(hits, key=lambda h: -h[1])


# Save

def save_panel(fig, stem, formats=("svg", "pdf", "png"), close=True):
    """Write a panel at exactly its figsize, in the formats the assembler needs."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    written = []
    # Marker "+" and powers of ten print as superscripts for the save only
    # (cnsfig.superscript).
    from cnsfig.superscript import superscripts
    with superscripts(fig, family=resolved_family()):
        for ext in formats:
            out = stem.with_suffix(f".{ext}")
            fig.savefig(out, dpi=matplotlib.rcParams["savefig.dpi"],
                        facecolor="white"
                        if not matplotlib.rcParams["savefig.transparent"] else "none")
            written.append(out)
    if close:
        plt.close(fig)
    return written


if __name__ == "__main__":
    if "--dump-snapshot" in sys.argv:
        apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
        print(dump_snapshot(title_fontsize=7, fontsize_legend=6,
                            legend_fontsize=6))
        sys.exit(0)
    from panel_style_cns_selftest import self_test
    sys.exit(self_test())
