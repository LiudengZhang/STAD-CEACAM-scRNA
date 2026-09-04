#!/usr/bin/env python3
"""
The single entry point for the cnsplots restyle (Version B).

Authorised by the author on 2026-09-01: the main figures are redrawn rather
than patched, and every figure is set in one type system, taken from
https://github.com/faridrashidi/cnsplots.

WHAT THIS MODULE IS
-------------------
A thin adapter over `cnsplots`. It does not hold a copy of cnsplots' numbers.
Every size, weight, spine and tick setting is read out of `cns.settings` at
call time and applied by `cns.setup_matplotlib()`, which is cnsplots' own
mechanism. If cnsplots changes its defaults, this module changes with it.

What the module adds on top of cnsplots is only what cnsplots has no opinion
about, because it is specific to this project:

  1. a font stack that resolves on this machine (see FONTS below),
  2. millimetre geometry, so a panel is drawn at the size it prints at,
  3. an exact-size save, because cnsplots' `savefig.bbox='tight'` crops the
     canvas to the ink and destroys the 1:1 relationship,
  4. a staged rollout switch for colour (see COLOUR below).

WHY THE PANELS ARE DRAWN AT PRINT SIZE
--------------------------------------
Measured on the shipped supplementary figures on 2026-09-01: every panel SVG
is drawn 3.5x to 6.3x larger than the box the assembler puts it in, then fitted
into that box with `preserveAspectRatio='xMidYMid meet'`. The fit scale runs
from 0.1600 (S8H) to 0.4538 (S7A), a 2.84x spread. A panel that sets 20 pt type
- which is what the project's 4x convention produces from a nominal 5 pt -
prints at 3.2 pt in S8H and 9.1 pt in S7A. That is the whole of the typography
problem: the number that reaches the reader is the product of the script's font
size, the script's SCALE, and the assembler's fit.

So a restyled panel is drawn at exactly the millimetre box it will occupy. The
fit scale becomes 1.0, and the size set here is the size printed. Nothing else
makes "8 pt" mean 8 pt.

    print_pt = set_pt x SCALE x assembler_fit          <- Version A
    print_pt = set_pt                                  <- Version B

FONTS
-----
cnsplots asks for Helvetica > Helvetica Neue > Arial > DejaVu Sans. Neither
Helvetica nor Arial is installed on this machine, so cnsplots' own stack falls
through to DejaVu Sans, which is neither metrically nor visually Helvetica and
is markedly wider. Two Helvetica-metric faces *are* installed:

    Nimbus Sans      URW's Helvetica clone, metrically identical to Helvetica
    Liberation Sans  metrically identical to Arial, itself a Helvetica metric
                     clone; already 2072 of the 2954 shipped text spans

Both carry regular, bold, italic and bold-italic and both embed cleanly under
`pdf.fonttype=42` (tested 2026-09-01). So the stack is cnsplots' own order with
those two inserted ahead of the DejaVu fallback, set through
`cns.settings.font_sans_serif`. On a machine with real Helvetica installed,
Helvetica still wins and nothing here changes.

COLOUR
------
`cns.setup_matplotlib()` also sets `axes.prop_cycle` and `image.cmap` to
cnsplots' palettes. Palette unification is stage 4 of the restyle and is
deliberately not part of the type work, so by default this module restores the
incoming colour cycle and colormap after cnsplots has run, and says so in
`describe()`. Pass `palette=True` in stage 4 to let cnsplots' palettes through.
This is a staged rollout, not a workaround: one switch, defaulted off, removed
when stage 4 lands.

USAGE
-----
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[N] / "00_Config"))
    import panel_style_cns as style

    style.apply()
    fig, ax = style.subplots_mm(100, 26)
    ...  # the panel's own drawing code, unchanged
    style.save_panel(fig, out_dir / "S8_F_leave_one_out")

Run this file directly for a self-test:

    conda run -n Liudeng_Python_310 python 00_Config/panel_style_cns.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

import cnsplots as cns

__all__ = [
    "apply", "describe", "resolved_family",
    "PAGE_W_MM", "PAGE_W_PT", "MM_TO_INCH",
    "figsize_mm", "figure_mm", "subplots_mm", "margins_mm", "overflow_mm",
    "save_panel", "body_pt", "tick_pt", "letter_pt", "self_test",
]

# ---------------------------------------------------------------------------
# Page geometry
# ---------------------------------------------------------------------------
# cnsplots expresses a full-width figure as `multipanel_max_width` "pixels",
# and its pixel is a 72-dpi point, so 540 px = 540 pt = 190.5 mm. Read it from
# the setting rather than writing 540 here, so the two cannot drift apart.
MM_TO_INCH = 1.0 / 25.4
PT_PER_MM = 72.0 / 25.4


def _page_width_pt() -> float:
    return float(cns.settings.multipanel_max_width)


PAGE_W_PT = _page_width_pt()
PAGE_W_MM = PAGE_W_PT / PT_PER_MM          # 190.5 mm

# ---------------------------------------------------------------------------
# Font stack
# ---------------------------------------------------------------------------
# cnsplots' own preference order, with the locally installed Helvetica-metric
# faces inserted ahead of the DejaVu fallback. Keep DejaVu Sans last: it is the
# only face here with full coverage of the mathtext and symbol glyphs some
# panels use, and dropping it would turn a missing glyph into a silent box.
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
    """The first family in the stack that this machine can actually draw.

    matplotlib resolves the stack silently, so a figure can be set in DejaVu
    Sans while the rcParam still says Helvetica. Thirteen font families ship in
    the current figures partly for this reason. Ask the font manager instead.
    """
    stack = stack or tuple(cns.settings.font_sans_serif)
    installed = {f.name for f in fm.fontManager.ttflist}
    for name in stack:
        if name in installed:
            return name
    return "DejaVu Sans"


# ---------------------------------------------------------------------------
# Apply
# ---------------------------------------------------------------------------

def _stack_led_by_resolved() -> tuple[str, ...]:
    """The font stack with the family that will actually be drawn put first.

    `svg.fonttype='none'` keeps SVG text editable, which the project requires -
    but it means matplotlib writes the *whole* stack into every <text> element
    and leaves the choice to whatever renders the SVG next. The assembler
    renders it with cairo, and cairo does not resolve a CSS font list the way
    matplotlib's font manager does: given the same stack, matplotlib drew Nimbus
    Sans and cairo drew Liberation Sans, so a panel and the page it was placed
    on ended up in different faces. That is one of the mechanisms behind the
    thirteen families in the shipped figures.

    Naming the resolved family first removes the ambiguity without changing what
    is asked for: `resolved_family()` returns the first *installed* entry, so
    promoting it is a no-op on any machine except in what it tells a downstream
    renderer. On a machine with real Helvetica, Helvetica is resolved and stays
    first.
    """
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
        `image.cmap`, so this is a type-only change. True lets cnsplots'
        palettes through; that is stage 4 of the restyle.
    **setting_overrides
        Written onto `cns.settings` before `setup_matplotlib()` runs, and
        validated by cnsplots. An unknown name raises, which is the point.

    Returns
    -------
    str
        The font family that will actually be drawn on this machine.
    """
    cns.settings.font_sans_serif = _extended_font_stack()
    # Put the family that will actually be drawn at the head of the stack, so
    # every renderer downstream of the SVG picks the same one. See
    # _stack_led_by_resolved().
    cns.settings.font_sans_serif = _stack_led_by_resolved()

    # Exact-size output. cnsplots defaults to bbox='tight', which crops the
    # canvas to the ink: the saved SVG's viewBox then no longer matches the
    # figsize, the assembler's fit scale is no longer 1.0, and the type stops
    # printing at the size it was set at. 'standard' is matplotlib's name for
    # "save the canvas I asked for".
    cns.settings.savefig_bbox = "standard"

    # The panels have always been written on white and composited onto a white
    # page; keep that rather than introduce transparency as a side effect of a
    # type change.
    cns.settings.savefig_transparent = False

    # The author's spec puts tick labels and legend text together at 7 pt.
    # cnsplots drives tick labels from fontsize_legend (7) but lets legend text
    # inherit title_fontsize (8) unless legend_fontsize is set. Set it.
    if cns.settings.legend_fontsize is None:
        cns.settings.legend_fontsize = cns.settings.fontsize_legend

    for key, value in setting_overrides.items():
        setattr(cns.settings, key, value)

    keep_cycle = matplotlib.rcParams["axes.prop_cycle"]
    keep_cmap = matplotlib.rcParams["image.cmap"]

    cns.setup_matplotlib()

    if not palette:
        matplotlib.rcParams["axes.prop_cycle"] = keep_cycle
        matplotlib.rcParams["image.cmap"] = keep_cmap

    # ps.fonttype is not in cnsplots' contract but the project's panels write
    # EPS in a few places and Type 3 there would undo pdf.fonttype=42.
    matplotlib.rcParams["ps.fonttype"] = 42

    return resolved_family()


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
    """(family, weight) for a panel letter, from cnsplots' settings.

    Call `apply()` first. Before it runs, `cns.settings.font_sans_serif` is
    still cnsplots' own stack, which on this machine resolves to DejaVu Sans -
    and an assembler that asked for the letter font before applying the style
    duly drew its panel letters in DejaVu Sans while the panels were in Nimbus
    Sans.
    """
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


# ---------------------------------------------------------------------------
# Millimetre geometry
# ---------------------------------------------------------------------------

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

    Version A's panels set their margins as *fractions* of a canvas four times
    the printed size, so `left=0.12` reserved 48 mm on a 400 mm canvas. Carried
    onto a 100 mm canvas unchanged, the same fraction reserves 12 mm and the
    tick labels run off the page - which is exactly what happened the first
    time S8F was redrawn. A margin is a physical distance; express it as one.

    Any margin left as None keeps whatever the figure already has.
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
    """Every artist that puts marks outside its own axes box.

    Explicit, because neither of the obvious shortcuts works:

      `fig.get_children()` - an axes reports only its own rectangle, so a tick
      label, an axis label, a title or a legend hanging off the edge is
      invisible. S10C's x axis label was clipped 3.1 mm off the right edge and
      nothing was reported.

      `fig.get_tightbbox()` - clamps the result. Measured on a deliberately
      over-long axis label that ran 12.4 mm past the right edge, it reported
      0.4 mm. It is built for layout, not for detecting a clip.

    So the artists are listed. Tick labels whose tick sits outside the current
    view interval are skipped: the locator makes them, the canvas never draws
    them, and counting them reports an overflow that no reader will ever see.
    """
    out = list(fig.texts)                       # suptitle, supxlabel, supylabel
    for ax in fig.axes:
        out.append(ax)
        out.extend([ax.title, ax.xaxis.label, ax.yaxis.label])
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
    """How far the drawn ink falls outside the canvas, in mm, per side.

    Returns (left, right, bottom, top); zeros mean nothing is clipped.

    A panel drawn 1:1 has no assembler shrink left to hide an overflow, so this
    is the check that replaces "it looked fine at 4x". Run it before saving.
    """
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


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

def save_panel(fig, stem, formats=("svg", "pdf", "png"), close=True):
    """Write a panel at exactly its figsize, in the formats the assembler needs.

    No `bbox_inches`. Cropping to the ink is what breaks 1:1 - it is also what
    the 4x scripts already avoid, for the same reason: they set their own
    margins with subplots_adjust and a tight box would undo them.
    """
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for ext in formats:
        out = stem.with_suffix(f".{ext}")
        fig.savefig(out, dpi=matplotlib.rcParams["savefig.dpi"],
                    facecolor="white"
                    if not matplotlib.rcParams["savefig.transparent"] else "none")
        written.append(out)
    if close:
        plt.close(fig)
    return written


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def self_test(tmpdir=None) -> int:
    """Assert the style took, and that a saved panel is exactly its own size.

    The second half is the one that matters. Every previous attempt at
    consistent type in this project failed not because the rcParams were wrong
    but because something downstream rescaled the canvas afterwards.
    """
    import re
    import tempfile

    matplotlib.use("Agg")
    fam = apply()
    rc = matplotlib.rcParams
    fails = []

    def check(name, got, want):
        if got != want:
            fails.append(f"{name}: got {got!r}, want {want!r}")

    check("axes.titlesize", rc["axes.titlesize"], cns.settings.title_fontsize)
    check("axes.labelsize", rc["axes.labelsize"], cns.settings.title_fontsize)
    check("axes.titleweight", rc["axes.titleweight"],
          cns.settings.title_fontweight)
    check("xtick.labelsize", rc["xtick.labelsize"], cns.settings.fontsize_legend)
    check("ytick.labelsize", rc["ytick.labelsize"], cns.settings.fontsize_legend)
    check("legend.fontsize", rc["legend.fontsize"], cns.settings.fontsize_legend)
    check("axes.linewidth", rc["axes.linewidth"], cns.settings.axes_linewidth)
    check("xtick.major.size", rc["xtick.major.size"],
          cns.settings.xtick_major_size)
    check("xtick.major.width", rc["xtick.major.width"],
          cns.settings.xtick_major_width)
    check("xtick.major.pad", rc["xtick.major.pad"], cns.settings.xtick_major_pad)
    check("axes.spines.top", rc["axes.spines.top"], False)
    check("axes.spines.right", rc["axes.spines.right"], False)
    check("legend.frameon", rc["legend.frameon"], False)
    check("legend.markerscale", rc["legend.markerscale"],
          cns.settings.legend_markerscale)
    check("svg.fonttype", rc["svg.fonttype"], "none")
    check("pdf.fonttype", rc["pdf.fonttype"], 42)
    # matplotlib normalises savefig.bbox='standard' to None, which is its way
    # of saying "no bbox adjustment". 'tight' here would mean the canvas gets
    # cropped to the ink and 1:1 is gone.
    if rc["savefig.bbox"] not in (None, "standard"):
        fails.append(f"savefig.bbox: got {rc['savefig.bbox']!r}, want no "
                     f"bbox adjustment")

    if fam == "DejaVu Sans":
        fails.append("font stack resolved to DejaVu Sans - no Helvetica-metric "
                     "face found; type will be wider than the target")

    # 1:1 geometry, through a real save.
    tmp = Path(tmpdir or tempfile.mkdtemp())
    w_mm, h_mm = 100.0, 26.0
    fig, ax = subplots_mm(w_mm, h_mm)
    ax.plot([0, 1], [0, 1])
    ax.set_xlabel("x label")
    save_panel(fig, tmp / "selftest")

    svg = (tmp / "selftest.svg").read_text()[:1500]
    m = re.search(r'viewBox="[\d.]+ [\d.]+ ([\d.]+) ([\d.]+)"', svg)
    if not m:
        fails.append("saved SVG has no viewBox")
    else:
        vw, vh = float(m.group(1)), float(m.group(2))
        want_w, want_h = w_mm * PT_PER_MM, h_mm * PT_PER_MM
        if abs(vw - want_w) > 0.5 or abs(vh - want_h) > 0.5:
            fails.append(f"saved SVG is {vw:.1f} x {vh:.1f} pt, asked for "
                         f"{want_w:.1f} x {want_h:.1f} pt - the 1:1 "
                         f"relationship is broken")

    try:
        import fitz
        doc = fitz.open(tmp / "selftest.pdf")
        page = doc[0]
        got = (page.rect.width, page.rect.height)
        want = (w_mm * PT_PER_MM, h_mm * PT_PER_MM)
        if abs(got[0] - want[0]) > 0.5 or abs(got[1] - want[1]) > 0.5:
            fails.append(f"saved PDF page is {got[0]:.1f} x {got[1]:.1f} pt, "
                         f"asked for {want[0]:.1f} x {want[1]:.1f} pt")
        sizes = {round(s["size"], 2)
                 for b in page.get_text("dict")["blocks"]
                 for l in b.get("lines", []) for s in l["spans"]}
        fonts = {s["font"].split("+")[-1]
                 for b in page.get_text("dict")["blocks"]
                 for l in b.get("lines", []) for s in l["spans"]}
        if not sizes <= {body_pt(), tick_pt()}:
            fails.append(f"printed sizes {sorted(sizes)} are not "
                         f"{{{tick_pt():g}, {body_pt():g}}}")
        if len(fonts) > 1:
            fails.append(f"more than one font family printed: {sorted(fonts)}")
        doc.close()
    except ImportError:
        print("  (PyMuPDF not available - PDF checks skipped)")

    print(describe())
    print()
    if fails:
        for f in fails:
            print(f"  FAIL  {f}")
        print(f"\n{len(fails)} check(s) failed")
        return 1
    print("  all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(self_test())
