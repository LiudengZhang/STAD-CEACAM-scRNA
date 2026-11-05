"""
Legends and dot-size keys.

dot_size_key lays out size-key circles from their radii so they never touch
(scanpy's own key spaces them one data unit apart, which merges large
circles in a narrow column). require_size_key raises if a dot plot varies dot
area without a key. compact_key and scanpy_compact_key stack the size key and
colour bar in one narrow column; group_key draws a colour key with
fixed-size circles.

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
           "MIN_GAP_PT"]

# Minimum gap between two size-key circles, pt. Below about 1 pt they read as
# one blob.
MIN_GAP_PT = 1.2


def _radii_pt(areas_pt2):
    """Circle radii in points, from scatter `s` values."""
    return [float(np.sqrt(max(a, 0.0)) / 2.0) for a in areas_pt2]


def circle_layout(areas_pt2, labels=None, label_pt=6.0, min_gap_pt=MIN_GAP_PT,
                  min_pitch_pt=0.0):
    """Centres, in points, for circles that must not touch one another."""
    radii = _radii_pt(areas_pt2)
    if not radii:
        return [], 0.0
    widths = [0.0] * len(radii)
    if labels is not None:
        # Approximate label width: 0.6 em per character.
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

    orientation="auto" lays the key across when it fits and down the axis
    otherwise, and raises if neither fits at the dots' own scale.
    """
    if len(areas_pt2) != len(labels):
        raise ValueError(f"{len(areas_pt2)} sizes but {len(labels)} labels")

    radii = _radii_pt(areas_pt2)
    biggest = max(radii) if radii else 1.0
    _, total = circle_layout(areas_pt2, labels, label_pt, min_gap_pt)
    # Vertically, labels sit beside the circles: one text line per step.
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

    # One data unit = one point, so separations computed in points match the
    # marker sizes, which scatter() also gives in points.
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
    """The areas and per-cent labels scanpy would use for `dotplot`'s key."""
    lo, hi = float(dotplot.dot_min), float(dotplot.dot_max)
    # A 0% step would draw an invisible circle, so with a zero floor the steps
    # start one interval up.
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


def require_size_key(axes, dot_areas=None, *, panel=""):
    """Raise unless a panel that encodes a quantity in dot area decodes it."""
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


# Compact dot-plot key: size circles over the colour bar in one narrow column
# at the right, instead of scanpy's 38.1 mm legend column.

def compact_key(fig, *, size_areas, size_labels, column_mm, top_mm=0.0,
                right_mm=0.0, size_title=None, cmap=None, vmin=0.0, vmax=1.0,
                cbar_title=None, cbar_ticks=(0.0, 0.5, 1.0), label_pt=6.0,
                title_pt=6.0, colour="gray", edge_colour="black",
                edge_lw=EDGE_PT, cbar_h_mm=1.5, gap_mm=0.8, linespacing=1.15):
    """Size key over colour bar, in a column `column_mm` wide at the right."""
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
        # Shortest tick form ("0", "0.5", "1") so labels fit a short bar.
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
    reserved: scanpy's own key axes are removed and replaced in place."""
    fig = dotplot.fig
    main = axes["mainplot_ax"]
    for key in ("size_legend_ax", "color_legend_ax"):
        ax = axes.get(key)
        if ax is not None:
            ax.remove()
    # Colour range as scanpy uses it: dot_color_df, or vmin/vmax when set.
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


# Group (colour) key circles print at a fixed diameter, mm, whatever size the
# data markers are.
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


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        from legend_selftest import self_test
        sys.exit(self_test())
    sys.exit(0)
