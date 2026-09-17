"""
The one place a legend is drawn, because everywhere else got it wrong.

Across the 58 panel scripts of Figures 2-5 there was no shared legend code at
all: every legend was hand-rolled, in four different idioms. Two of them were
wrong in ways nothing could see.

**Figure 5F printed its size key as a solid grey cone.** scanpy's
`_plot_size_legend` puts the circles at `np.arange(n) + 0.5` - a pitch of one
data unit - and sizes them from `largest_dot`, which defaults to 200 pt^2, a
circle 15.96 pt across. In a legend column 8.9 mm wide holding six of them, one
data unit is about 4.1 pt. A 16 pt circle every 4.1 pt is not a key; it is a
wedge. The pitch has to come from the radii, which is what `dot_size_key` does.

**Figure 4E printed no size key at all.** `legend(show_size_legend=False)`, on
the reasoning that the published page carried none and adding one would be new
content. The consequence was recorded in PROVENANCE.csv and left standing: "the
dot area encodes the fraction of expressing cells in each group, and with no
size key nothing on the page decodes it." The author ruled on 2026-09-11 that a
key decoding dots already drawn adds no data, so it is not new content.

Hence `require_size_key`, which raises. A dotplot that encodes a fraction in
dot area and reaches the page without a key is now a crash, not an observation
somebody writes down.

    python legend.py --self-test
"""

import sys

import numpy as np

_CONFIG = __import__("pathlib").Path(__file__).resolve().parent.parent
if str(_CONFIG) not in sys.path:
    sys.path.insert(0, str(_CONFIG))
from panel_style_cns import RULE_PT, EDGE_PT, PT_PER_MM  # noqa: E402

__all__ = ["dot_size_key", "require_size_key", "circle_layout",
           "compact_key", "scanpy_compact_key", "group_key", "GROUP_KEY_MM",
           "MIN_GAP_PT", "self_test"]

#: Clear paper between two circles of a size key. Below about a point they
#: read as one blob, which is the shape of the Figure 5F defect; the gate in
#: 10_Reproduction/check_restyled_panel.py holds text to 0.5 pt for the same
#: reason. Circles are given more room than text because a circle's ink runs
#: to the edge of its box while a glyph's does not.
MIN_GAP_PT = 1.2


def _radii_pt(areas_pt2):
    """Circle radii in points, from scatter `s` values.

    matplotlib's `s` is the marker's diameter squared in points, not pi r^2:
    a dot drawn with s=200 is 14.14 pt across. Until 2026-09-14 this took
    sqrt(s / pi) and every key circle was drawn 13% larger than the dot it
    named; the self-test control on scanpy's 200 pt2 dot pinned the wrong
    number (15.96) to it.
    """
    return [float(np.sqrt(max(a, 0.0)) / 2.0) for a in areas_pt2]


def circle_layout(areas_pt2, labels=None, label_pt=6.0, min_gap_pt=MIN_GAP_PT,
                  min_pitch_pt=0.0):
    """Centres, in points, for circles that must not touch one another.

    The pitch between two neighbours is whatever the two of them need - their
    radii plus the gap - and never a constant, because a constant is what
    produces a cone when the radii grow. If labels are given, the pitch also
    has to hold the two half-labels, so that fixing the circles does not simply
    move the collision into the text.

    `min_pitch_pt` is a floor on the distance between two neighbours whatever
    their radii. A vertical key needs one: its labels sit beside the circles,
    so their widths do not set the run, but two stacked labels still need a
    line of height between them. Without it the two smallest steps of Figure
    5F's key - circles about 2 pt across - were laid out 2 pt apart and their
    labels printed on top of one another. The horizontal branch passes 0,
    because there the label widths already do this job.

    Returns (centres, total_width_pt).
    """
    radii = _radii_pt(areas_pt2)
    if not radii:
        return [], 0.0
    widths = [0.0] * len(radii)
    if labels is not None:
        # 0.6 em per digit is the usual advance for lining figures in a
        # sans face; the gate measures the real thing afterwards.
        widths = [0.6 * label_pt * len(str(t)) for t in labels]

    centres = [radii[0] + widths[0] / 2.0]
    for i in range(1, len(radii)):
        need_circle = radii[i - 1] + radii[i] + min_gap_pt
        need_label = (widths[i - 1] + widths[i]) / 2.0 + min_gap_pt
        centres.append(centres[-1] + max(need_circle, need_label,
                                         min_pitch_pt))
    total = centres[-1] + max(radii[-1], widths[-1] / 2.0)
    return centres, total


def axis_size_pt(ax):
    """The axis's drawn width and height, in points."""
    fig = ax.get_figure()
    box = ax.get_position()
    return (box.width * fig.get_figwidth() * 72.0,
            box.height * fig.get_figheight() * 72.0)


def dot_size_key(ax, areas_pt2, labels, *, title=None, label_pt=6.0,
                 title_pt=7.0, colour="gray", edge_colour="black",
                 edge_lw=0.25, min_gap_pt=MIN_GAP_PT, orientation="auto"):
    """Draw a dot-size key into `ax`, laid out so the circles cannot merge.

    The subtlety that produced the shipped cone: `scatter(s=...)` sizes a
    marker in absolute points squared, independent of the data limits, so the
    circles do NOT shrink when the axis is narrow. Spacing them correctly in
    data coordinates is therefore not enough - the axis has to be physically
    big enough to hold them, and when it is not, no choice of xlim will help.

    Figure 5F's legend column is 8.9 mm wide and four circles up to 15.96 pt
    across need about 60 pt, which is 21 mm. It will never fit across. It fits
    easily down the column, which is 60 mm tall. So `orientation="auto"`
    measures the axis and stacks the key vertically when across will not do,
    and raises if neither direction has the room - because the alternative is
    to shrink the circles, and a key drawn at a different scale from the dots
    it decodes is worse than no key at all.
    """
    if len(areas_pt2) != len(labels):
        raise ValueError(f"{len(areas_pt2)} sizes but {len(labels)} labels")

    radii = _radii_pt(areas_pt2)
    biggest = max(radii) if radii else 1.0
    _, total = circle_layout(areas_pt2, labels, label_pt, min_gap_pt)
    # Down the column the labels sit beside the circles, so their widths do
    # not set the run - but their heights do, one line per step.
    line_pt = label_pt * 1.25
    _, total_v = circle_layout(areas_pt2, None, label_pt, min_gap_pt,
                               min_pitch_pt=line_pt)
    w_pt, h_pt = axis_size_pt(ax)

    if orientation == "auto":
        orientation = "h" if total <= w_pt else "v"
    if orientation == "h" and total > w_pt:
        raise RuntimeError(
            f"a horizontal size key needs {total:.1f} pt and the axis is "
            f"{w_pt:.1f} pt wide; the circles would merge")
    if orientation == "v" and total_v > h_pt:
        raise RuntimeError(
            f"a vertical size key needs {total_v:.1f} pt and the axis is "
            f"{h_pt:.1f} pt tall; the circles would merge")

    ax.clear()
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.grid(False)
    ax.tick_params(axis="both", left=False, bottom=False, right=False,
                   top=False, labelleft=False, labelbottom=False,
                   labelright=False, labeltop=False)

    # ONE DATA UNIT IS ONE POINT. This is the whole mechanism, and getting it
    # wrong is what produced the second wrong key as well as the first:
    # scatter() sizes its markers in absolute points, so if the data range is
    # set to anything other than the axis's own size in points, the computed
    # separations are scaled while the circles are not, and they close up
    # again. Setting the limits to the axis's measured extent makes the layout
    # computed above the layout drawn.
    ax.set_xlim(0, w_pt)
    ax.set_ylim(0, h_pt)

    if orientation == "h":
        centres, _ = circle_layout(areas_pt2, labels, label_pt, min_gap_pt)
        off = max((w_pt - total) / 2.0, 0.0)
        xs = [c + off for c in centres]
        y = h_pt / 2.0 + label_pt * 0.6      # room for the labels beneath
        ax.scatter(xs, [y] * len(xs), s=list(areas_pt2), color=colour,
                   edgecolor=edge_colour, linewidth=edge_lw, zorder=100,
                   clip_on=False)
        for x, r, t in zip(xs, radii, labels):
            ax.text(x, y - r - 1.0, str(t), ha="center", va="top",
                    fontsize=label_pt, clip_on=False)
    else:
        centres, _ = circle_layout(areas_pt2, None, label_pt, min_gap_pt,
                                   min_pitch_pt=line_pt)
        off = max((h_pt - total_v) / 2.0, 0.0)
        # Largest at the top, which is how a size key reads.
        ys = [h_pt - off - c for c in centres]
        x = biggest + 0.5
        ax.scatter([x] * len(ys), ys, s=list(areas_pt2), color=colour,
                   edgecolor=edge_colour, linewidth=edge_lw, zorder=100,
                   clip_on=False)
        for yy, r, t in zip(ys, radii, labels):
            # clip_on=False so a label may sit in the panel margin; the panel
            # reserves that room through style.fit_margins.
            ax.text(x + biggest + 1.5, yy, str(t), ha="left", va="center",
                    fontsize=label_pt, clip_on=False)

    if title:
        ax.set_title(title, fontsize=title_pt)
    return ax


def scanpy_dot_areas(dotplot, n=4):
    """The areas and per-cent labels scanpy would use for `dotplot`'s key.

    Read off the DotPlot rather than recomputed, so the key describes the dots
    actually drawn: same dot_min/dot_max, same size_exponent, same
    largest_dot/smallest_dot.
    """
    lo, hi = float(dotplot.dot_min), float(dotplot.dot_max)
    # scanpy's scale starts at dot_min, which is usually 0, and a 0% step
    # draws a circle of zero area - a step of the key that shows the reader
    # nothing. When the floor is zero the steps start one interval up, so
    # every circle printed is a circle that can be seen and compared.
    if lo <= 0:
        fracs = np.linspace(hi / n, hi, n)
    else:
        fracs = np.linspace(lo, hi, n)
    span = hi - lo
    vals = (fracs - lo) / span if span > 0 else np.ones_like(fracs)
    areas = (vals ** dotplot.size_exponent
             * (dotplot.largest_dot - dotplot.smallest_dot)
             + dotplot.smallest_dot)
    return list(areas), [int(round(f * 100)) for f in fracs]


def grow_axis(ax, want_pt, *, avoid=(), pad_pt=3.0):
    """Give `ax` `want_pt` of height by taking free space below it.

    A size key is as tall as its circles need, and the axis a plotting library
    hands back is whatever that library felt like. Figure 5F's legend column
    is 61.6 mm tall with the size key and the colour bar at opposite ends and
    daylight between them, but scanpy gave the key 31.6 pt and a key needing
    37.1 pt then had nowhere to go. Rather than shrink the circles - which
    would make the key lie about the dots - the axis takes the empty space.

    The top edge stays put; the bottom extends, and stops `pad_pt` short of
    whatever is below it. Returns the height actually obtained, which may be
    less than asked for; the caller still has to check it was enough.
    """
    fig = ax.get_figure()
    h_fig_pt = fig.get_figheight() * 72.0
    box = ax.get_position()
    top = box.y1
    floor = 0.0
    for other in avoid:
        if other is ax:
            continue
        obox = other.get_position()
        if obox.y1 <= box.y0 + 1e-9:            # strictly below
            floor = max(floor, obox.y1)
    floor += pad_pt / h_fig_pt
    want_frac = want_pt / h_fig_pt
    new_y0 = max(floor, top - want_frac)
    ax.set_position([box.x0, new_y0, box.width, top - new_y0])
    return (top - new_y0) * h_fig_pt


def fit_size_key(ax, dotplot, *, max_steps=5, min_steps=3, **kw):
    """Draw the largest dot-size key that fits the axis, and say which it drew.

    How many steps a key can hold is a property of the axis, not a number to
    guess: the circles are sized in absolute points from the dotplot's own
    scale and may not be shrunk, so on a short axis the only thing that can
    give is how many of them there are. Four steps need 38.3 pt down Figure
    5F's legend column and the column is 31.6 pt, which is why this exists.

    Fewer than `min_steps` is not a key worth printing - two circles name a
    range without showing the scale between them - so it raises rather than
    quietly degrading to something useless.
    """
    tried = []
    for n in range(max_steps, min_steps - 1, -1):
        areas, labels = scanpy_dot_areas(dotplot, n=n)
        # A step whose circle is under about a point across draws nothing a
        # reader can see, and a key step that shows no circle is not a step.
        # scanpy's scale starts at dot_min, which is often 0.
        keep = [(a, t) for a, t in zip(areas, labels)
                if 2 * np.sqrt(max(a, 0.0) / np.pi) >= 1.0]
        if len(keep) < min_steps:
            tried.append(f"{n}: only {len(keep)} steps draw a visible circle")
            continue
        areas = [a for a, _ in keep]
        labels = [t for _, t in keep]
        try:
            dot_size_key(ax, areas, labels, **kw)
            return areas, labels, len(areas)
        except RuntimeError as exc:
            tried.append(f"{n}: {exc}")
    raise RuntimeError(
        "no dot-size key fits this axis, down to "
        f"{min_steps} steps. Give the legend more room rather than shrinking "
        "the circles - a key drawn at a different scale from the dots it "
        "decodes is worse than none.\n  " + "\n  ".join(tried))


def require_size_key(axes, dot_areas=None, *, panel=""):
    """Raise unless a panel that encodes a quantity in dot area decodes it.

    `axes` is scanpy's DotPlot.get_axes() mapping, or anything with the same
    keys. `dot_areas` is the set of marker areas drawn on the main plot; when
    they are all equal the dot area carries nothing and no key is owed.
    """
    if dot_areas is not None:
        distinct = {round(float(a), 6) for a in np.ravel(dot_areas)}
        if len(distinct) <= 1:
            return                      # constant dots encode nothing
    ax = axes.get("size_legend_ax") if hasattr(axes, "get") else None
    if ax is None:
        raise RuntimeError(
            f"{panel or 'this panel'} varies dot area but has no size legend "
            f"axis. The dot area encodes a fraction and nothing on the page "
            f"decodes it; see 00_Config/cnsfig/legend.py.")
    if not ax.collections and not ax.get_xticklabels():
        raise RuntimeError(
            f"{panel or 'this panel'} has an empty size legend axis. A key "
            f"that draws nothing is the same as no key.")


def self_test():
    """A legend rule nobody has watched fail is not a rule."""
    ok = True

    def control(sign, name, got, want):
        nonlocal ok
        good = got == want
        ok &= good
        print(f"  control {sign}  {name}: "
              f"{'as required' if good else f'WRONG got {got} want {want}'}")

    # The Figure 5F geometry: scanpy's own numbers, six steps up to
    # largest_dot=200, which is a 15.96 pt circle.
    areas = [12.5, 50.0, 112.5, 200.0]
    radii = _radii_pt(areas)
    centres, total = circle_layout(areas, labels=None)
    gaps = [centres[i + 1] - centres[i] - radii[i] - radii[i + 1]
            for i in range(len(centres) - 1)]
    control("0", "every neighbouring pair clears MIN_GAP_PT",
            all(g >= MIN_GAP_PT - 1e-9 for g in gaps), True)
    control("0", "the largest circle is 14.14 pt across, as scatter(s=200) "
            "draws it", round(2 * radii[-1], 2), 14.14)

    # A VERTICAL KEY IS PACED BY ITS LABELS WHEN ITS CIRCLES ARE SMALL
    #   The run down the column used to be computed from the radii alone, on
    #   the reasoning that the labels sit beside the circles rather than under
    #   them. True of their widths, false of their heights: once Figure 5F's
    #   largest dot came down to the published page's 2.63 mm, the two smallest
    #   steps were circles about 2 pt across, laid out 2 pt apart, and their
    #   labels '12' and '24' printed on top of each other. The control is the
    #   small-circle case, because the large-circle case passes either way.
    tiny = [1.5, 4.0, 9.0, 20.0]
    line_pt = 6.0 * 1.25
    paced, _ = circle_layout(tiny, None, 6.0, min_pitch_pt=line_pt)
    pitches = [paced[i + 1] - paced[i] for i in range(len(paced) - 1)]
    control("0", "a vertical key paces small circles by the label line",
            all(q >= line_pt - 1e-9 for q in pitches), True)
    unpaced, _ = circle_layout(tiny, None, 6.0)
    control("X", "without the floor the same key would stack labels",
            any(unpaced[i + 1] - unpaced[i] < line_pt
                for i in range(len(unpaced) - 1)), True)

    # The mutant is scanpy's own layout: a constant pitch of one data unit,
    # which is what produced the cone.
    unit = [i * 4.1 for i in range(len(areas))]
    unit_gaps = [unit[i + 1] - unit[i] - radii[i] - radii[i + 1]
                 for i in range(len(unit) - 1)]
    control("+", "a constant pitch DOES overlap, which is the shipped defect",
            any(g < 0 for g in unit_gaps), True)

    control("0", "labels widen the pitch when they are wider than the circles",
            circle_layout([1.0, 1.0], labels=["100", "100"], label_pt=7.0)[1]
            > circle_layout([1.0, 1.0], labels=None)[1], True)

    class FakeAx:
        collections, _ticks = [], []

        def get_xticklabels(self):
            return self._ticks

    missing = {}
    try:
        require_size_key(missing, dot_areas=[1.0, 9.0, 25.0], panel="probe")
        caught = False
    except RuntimeError:
        caught = True
    control("+", "a varying dotplot with no size axis raises", caught, True)

    try:
        require_size_key({"size_legend_ax": FakeAx()},
                         dot_areas=[1.0, 9.0], panel="probe")
        caught = False
    except RuntimeError:
        caught = True
    control("+", "an empty size legend axis raises", caught, True)

    try:
        require_size_key({}, dot_areas=[7.0, 7.0, 7.0], panel="probe")
        quiet = True
    except RuntimeError:
        quiet = False
    control("-", "constant dot areas owe no key and do NOT raise", quiet, True)

    # Everything above is arithmetic. The two keys that shipped wrong were
    # both drawn by code whose arithmetic was fine and whose RENDERING was
    # not - the first because scanpy spaced the circles in data units, the
    # second because this module set data limits that did not match the
    # axis, so the separations scaled and the circles did not. So the last
    # controls draw a real figure and measure the ink.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def render(orientation, w_mm, h_mm, areas, labels):
        fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4), dpi=600)
        ax = fig.add_axes([0, 0, 1, 1])
        dot_size_key(ax, areas, labels, orientation=orientation, label_pt=6.0)
        fig.canvas.draw()
        pts = ax.collections[0].get_offsets()
        trans = ax.transData.transform
        px = trans(np.asarray(pts))                     # device pixels
        scale = 72.0 / fig.dpi                          # pixels -> points
        texts = [t.get_text() for t in ax.texts]
        plt.close(fig)
        return px * scale, texts

    areas4 = [12.5, 50.0, 112.5, 200.0]
    labels4 = [10, 30, 50, 70]
    rad = _radii_pt(areas4)

    xy, texts = render("v", 16.0, 60.0, areas4, labels4)
    gaps = [abs(xy[i + 1][1] - xy[i][1]) - rad[i] - rad[i + 1]
            for i in range(len(rad) - 1)]
    control("0", "RENDERED vertical circles clear MIN_GAP_PT",
            all(g >= MIN_GAP_PT - 0.05 for g in gaps), True)
    control("0", "every vertical step is labelled", len(texts), len(labels4))

    xy, texts = render("h", 90.0, 18.0, areas4, labels4)
    gaps = [abs(xy[i + 1][0] - xy[i][0]) - rad[i] - rad[i + 1]
            for i in range(len(rad) - 1)]
    control("0", "RENDERED horizontal circles clear MIN_GAP_PT",
            all(g >= MIN_GAP_PT - 0.05 for g in gaps), True)
    control("0", "every horizontal step is labelled", len(texts), len(labels4))

    # The regression this control exists for: the same key drawn into an axis
    # too small for it must raise, never silently compress.
    try:
        render("v", 16.0, 12.0, areas4, labels4)
        caught = False
    except RuntimeError:
        caught = True
    control("+", "an axis too small to hold the key raises", caught, True)

    _group_key_controls(control)

    print("\nevery control behaved as required" if ok
          else "\n*** A CONTROL FAILED ***")
    return 0 if ok else 1



# ---------------------------------------------------------------------------
# One compact key for a dotplot: size circles and colour bar in one narrow
# column, since 2026-09-14.
# ---------------------------------------------------------------------------
# scanpy's DotPlot.legend() reserves 1.5 inches - 38.1 mm - at the right of
# the figure and sets the size key and the colour bar at opposite ends of it,
# with paper between. On Figure 4E that column was a third of a 114.9 mm
# panel, on 2G a framed ax.legend did the same at the right; the author read
# both as waste. compact_key stacks the two keys at the top of a column of
# the width the caller names, flush with the matrix, and hands back the same
# axes dict require_size_key reads.
#
# The circle radius is matplotlib's: scatter's `s` is the marker's diameter
# squared in points, so r = sqrt(s) / 2. _radii_pt used pi r^2 until
# 2026-09-14, which drew every key circle 13% larger than the dot it named.

def compact_key(fig, *, size_areas, size_labels, column_mm, top_mm=0.0,
                right_mm=0.0, size_title=None, cmap=None, vmin=0.0, vmax=1.0,
                cbar_title=None, cbar_ticks=(0.0, 0.5, 1.0), label_pt=6.0,
                title_pt=6.0, colour="gray", edge_colour="black",
                edge_lw=EDGE_PT, cbar_h_mm=1.5, gap_mm=0.8, linespacing=1.15):
    """Size key over colour bar, in a column `column_mm` wide at the right.

    Returns {"size_legend_ax", "color_legend_ax", "height_mm"}. Raises when
    the column cannot hold the circles in either direction: a key drawn at a
    different scale from the dots it decodes is worse than none.
    """
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    W = fig.get_figwidth() * 25.4
    H = fig.get_figheight() * 25.4
    x0 = W - right_mm - column_mm
    y = float(top_mm)

    def add(h_mm):
        return fig.add_axes([x0 / W, 1.0 - (y + h_mm) / H, column_mm / W,
                             h_mm / H])

    def heading(text):
        n = text.count("\n") + 1
        h = (n * title_pt * linespacing + 1.0) / PT_PER_MM
        ax = add(h)
        ax.set_axis_off()
        ax.text(0.5, 1.0, text, ha="center", va="top", fontsize=title_pt,
                linespacing=linespacing, transform=ax.transAxes)
        return h

    if size_title:
        y += heading(size_title) + 0.4 * gap_mm

    radii = _radii_pt(size_areas)
    biggest = max(radii)
    col_pt = column_mm * PT_PER_MM
    _, total_h = circle_layout(size_areas, size_labels, label_pt, MIN_GAP_PT)
    if total_h <= col_pt:
        orientation = "h"
        h_mm = (2 * biggest + 1.6 * label_pt + 2.0) / PT_PER_MM
    else:
        orientation = "v"
        widest = max(0.6 * label_pt * len(str(t)) for t in size_labels)
        if 2 * biggest + 1.5 + widest > col_pt:
            raise RuntimeError(
                f"a {column_mm} mm column holds the key neither across "
                f"({total_h:.1f} pt) nor down ({2 * biggest + 1.5 + widest:.1f}"
                f" pt wide); give it room rather than shrinking the circles")
        _, total_v = circle_layout(size_areas, None, label_pt, MIN_GAP_PT,
                                   min_pitch_pt=label_pt * 1.25)
        h_mm = (total_v + 2.0) / PT_PER_MM
    sax = add(h_mm)
    dot_size_key(sax, size_areas, size_labels, label_pt=label_pt,
                 colour=colour, edge_colour=edge_colour, edge_lw=edge_lw,
                 orientation=orientation)
    y += h_mm + gap_mm

    cax = None
    if cmap is not None:
        if cbar_title:
            y += heading(cbar_title) + 0.4 * gap_mm
        cax = add(cbar_h_mm)
        cb = matplotlib.colorbar.ColorbarBase(
            cax, cmap=plt.get_cmap(cmap), norm=Normalize(vmin, vmax),
            orientation="horizontal", ticks=list(cbar_ticks))
        # Shortest form of each tick ("0", "0.5", "1"): on a 10 mm bar the
        # three two-decimal labels touched.
        cb.set_ticks(list(cbar_ticks))
        cb.set_ticklabels([f"{t:g}" for t in cbar_ticks])
        cax.tick_params(labelsize=label_pt, width=RULE_PT, length=1.5,
                        pad=1.0)
        cb.outline.set_linewidth(RULE_PT)
        # The end labels are aligned inward, so the bar's width is the key's
        # width and nothing overhangs the column.
        labs = cax.get_xticklabels()
        if len(labs) >= 2:
            labs[0].set_ha("left")
            labs[-1].set_ha("right")
        y += cbar_h_mm + (label_pt * 1.3 + 1.5) / PT_PER_MM

    return {"size_legend_ax": sax, "color_legend_ax": cax,
            "height_mm": y - float(top_mm)}


def scanpy_compact_key(dotplot, axes, column_mm, *, size_title, cbar_title,
                       label_pt=6.0, title_pt=6.0, n=4, cbar_ticks=None):
    """compact_key for a scanpy DotPlot whose legend column is already
    reserved: scanpy's own key axes are removed and replaced in place.

    Returns (axes_dict, areas) - pass both to require_size_key.
    """
    fig = dotplot.fig
    main = axes["mainplot_ax"]
    for key in ("size_legend_ax", "color_legend_ax"):
        ax = axes.get(key)
        if ax is not None:
            ax.remove()
    # scanpy colours the dots itself (the collection carries facecolours, not
    # an array and a norm), and draws its colour bar from dot_color_df with
    # its own vmin/vmax when set. The same numbers are read here.
    df = dotplot.dot_color_df
    vb = dotplot.vboundnorm
    vmin = float(vb.vmin) if vb.vmin is not None else float(df.min().min())
    vmax = float(vb.vmax) if vb.vmax is not None else float(df.max().max())
    cmap = dotplot.cmap
    areas, labels = scanpy_dot_areas(dotplot, n=n)
    keep = [(a, t) for a, t in zip(areas, labels)
            if np.sqrt(max(a, 0.0)) >= 1.0]
    areas = [a for a, _ in keep]
    labels = [t for _, t in keep]
    H = fig.get_figheight() * 25.4
    top_mm = (1.0 - main.get_position().y1) * H
    if cbar_ticks is None:
        cbar_ticks = (vmin, (vmin + vmax) / 2, vmax)
    # One millimetre of paper at the canvas edge; the end tick labels are
    # aligned inward, so nothing else overhangs the column.
    out = compact_key(fig, size_areas=areas, size_labels=labels,
                      column_mm=column_mm - 1.0, right_mm=1.0, top_mm=top_mm,
                      size_title=size_title, cmap=cmap,
                      vmin=vmin, vmax=vmax, cbar_title=cbar_title,
                      cbar_ticks=cbar_ticks, label_pt=label_pt,
                      title_pt=title_pt)
    return out, areas


# ---------------------------------------------------------------------------
# A group key whose circles print at a fixed size, since 2026-09-15.
# ---------------------------------------------------------------------------
# The Figure 3 scatter panels drew their keys from the plotted handles, so a
# key circle was the size of a data marker: 0.49 mm on the page shipped on
# 2026-09-15, which the author could not see, beside 6 pt words. A key names
# a colour, not a size, so its circle is set here at GROUP_KEY_MM whatever
# the data are drawn at.
GROUP_KEY_MM = 1.3


def group_key(fig, entries, *, x_mm, y_mm, key_mm=GROUP_KEY_MM, label_pt=None,
              linespacing=0.35, edge_pt=EDGE_PT, ncol=1, columnspacing=1.0):
    """(label, facecolor[, edgecolor]) keys in `ncol` columns, top-left at
    (x, y) mm from the canvas's top-left; an entry whose facecolor is None
    prints as an open circle. Returns the Legend."""
    import matplotlib.lines as mlines
    import matplotlib.pyplot as plt
    if label_pt is None:
        label_pt = float(plt.rcParams["legend.fontsize"]) \
            if isinstance(plt.rcParams["legend.fontsize"], (int, float)) \
            else float(plt.rcParams["font.size"])
    handles = []
    for e in entries:
        label, face = e[0], e[1]
        edge = e[2] if len(e) > 2 else (face if face is not None else "0.4")
        handles.append(mlines.Line2D(
            [], [], linestyle="none", marker="o",
            markersize=key_mm * PT_PER_MM,       # diameter in points
            markerfacecolor="none" if face is None else face,
            markeredgecolor=edge, markeredgewidth=edge_pt, label=label))
    w_mm = fig.get_figwidth() * 25.4
    h_mm = fig.get_figheight() * 25.4
    leg = fig.legend(handles=handles, loc="upper left",
                     bbox_to_anchor=(x_mm / w_mm, 1.0 - y_mm / h_mm),
                     frameon=False, borderpad=0.0, borderaxespad=0.0,
                     handlelength=0.7, handleheight=1.0, handletextpad=0.4,
                     labelspacing=linespacing, fontsize=label_pt,
                     ncol=ncol, columnspacing=columnspacing,
                     markerscale=1.0)          # never cnsplots' 0.x scale
    return leg


def _group_key_render(path, key_mm, markerscale=None):
    """Draw one key on a 30 x 20 mm canvas; return the circle diameters
    PyMuPDF reads back, in mm."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import fitz
    fig = plt.figure(figsize=(30 / 25.4, 20 / 25.4))
    leg = group_key(fig, [("Pre-R", "#74b9ff"), ("Alive", None, "0.5")],
                    x_mm=2.0, y_mm=2.0, key_mm=key_mm, label_pt=6.0)
    if markerscale is not None:            # the mutation: cnsplots' scale
        for h in leg.legend_handles:
            h.set_markersize(h.get_markersize() * markerscale)
    fig.savefig(path, format="pdf")
    plt.close(fig)
    page = fitz.open(path)[0]
    dias = []
    for d in page.get_drawings():
        r = d["rect"]
        if 0.2 < r.width / PT_PER_MM < 5 and abs(r.width - r.height) < 0.3:
            dias.append(round(r.width / PT_PER_MM, 2))
    return sorted(set(dias))


def _group_key_controls(control):
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as tmp:
        pdf = Path(tmp) / "k.pdf"
        dias = _group_key_render(pdf, GROUP_KEY_MM)
        control("0", f"group_key circles print at {GROUP_KEY_MM} mm",
                all(abs(d - GROUP_KEY_MM) < 0.15 for d in dias) and len(dias) >= 1,
                True)
        # MUTATION: the shipped defect - key circles at the data-marker size.
        small = _group_key_render(pdf, 0.49)
        control("+", "a 0.49 mm key is measurably smaller (the defect)",
                all(d < 0.7 for d in small) and len(small) >= 1, True)
        # MUTATION: cnsplots' legend.markerscale applied after the fact.
        scaled = _group_key_render(pdf, GROUP_KEY_MM, markerscale=0.5)
        control("+", "markerscale 0.5 would halve the key; group_key pins 1.0",
                all(d < GROUP_KEY_MM * 0.6 for d in scaled) and len(scaled) >= 1, True)


if __name__ == "__main__":
    sys.exit(self_test() if "--self-test" in sys.argv else 0)
