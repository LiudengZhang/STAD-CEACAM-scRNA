"""
Box plots, significance brackets and their frames.

    box_xlim(positions, width, clear)   x limits leaving `clear` box widths
                                        between the outer boxes and the spines
    finish_two_group(...)               frame an R-vs-NR box panel: ticks with
                                        counts, rich y label, title on the canvas
    draw_boxes(ax, data, positions, colors, width)
                                        draw the boxes; data goes to ax.boxplot
    bracket(fig, ax, x1, x2, y_max, y_rng, p, kind)
                                        a bracket and its P label; "omnibus"
                                        prints "ns" when P >= 0.05
    ylim_above(ax, text, pad_mm)        raise ylim until a label fits
    assert_no_points(ax, bxp)           raise if individual points are drawn

Bracket labels are placed by the bottom of their ink, not their baseline, so
a star and "P = 0.30" sit the same distance above the line.

    python boxes.py --self-test
"""

import sys

_CONFIG = __import__("pathlib").Path(__file__).resolve().parent.parent
if str(_CONFIG) not in sys.path:
    sys.path.insert(0, str(_CONFIG))
from panel_style_cns import PT_PER_MM  # noqa: E402

__all__ = ["box_xlim", "finish_two_group", "frame_fixed", "title_band_mm", "CLEAR",
           "draw_boxes", "bracket", "ylim_above", "assert_no_points", "FLIER_MM",
           "INK_GAP_MM"]

# Open-circle flier diameter, mm.
FLIER_MM = 0.79
# Gap between a bracket's top line and the bottom of its label's ink, mm.
INK_GAP_MM = 0.4

# Gap between the outermost box and the spine, in box widths (about 1.2 mm
# for 0.5-wide boxes at unit pitch on an 8 mm plotting box).
CLEAR = 0.6


def box_xlim(positions, width, clear=CLEAR):
    """(lo, hi) so the first and last boxes sit `clear` widths off the spines."""
    positions = list(positions)
    if not positions or width <= 0:
        raise ValueError("box_xlim needs positions and a positive width")
    return (min(positions) - width / 2 - clear * width,
            max(positions) + width / 2 + clear * width)


def title_band_mm(lines, pt, linespacing=1.15, pad_mm=0.8):
    """Height of the canvas band a title of `lines` lines at `pt` occupies."""
    return lines * pt * linespacing / PT_PER_MM + pad_mm


def finish_two_group(fig, ax, *, title, tick_labels, positions, width,
                     ylabel_markup, letter_cell, panel_w_mm, panel_h_mm,
                     title_lines=2, ylabel_pt=None, title_pt=None,
                     pad_mm=0.3, clear=CLEAR, ylocator=None,
                     top_extra_mm=0.0, bottom_mm=None, tick_rotation=0):
    """Frame an R/NR box axes and draw its title on the canvas."""
    import panel_style_cns as style
    from cnsfig.rich import rich_ylabel
    tick_pt = style.tick_pt()
    ylabel_pt = tick_pt if ylabel_pt is None else ylabel_pt
    title_pt = style.body_pt() if title_pt is None else title_pt

    ax.set_xticks(list(positions))
    if tick_rotation:
        ax.set_xticklabels(list(tick_labels), rotation=tick_rotation,
                           ha="right", rotation_mode="anchor")
    else:
        ax.set_xticklabels(list(tick_labels))
    ax.set_xlim(*box_xlim(positions, width, clear))
    if ylocator is not None:
        ax.yaxis.set_major_locator(ylocator)
    ax.set_title("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", pad=0.6, length=1.5)
    ax.tick_params(axis="x", pad=1.5)
    rich_ylabel(ax, ylabel_markup, fontsize=ylabel_pt, labelpad_pt=1.0)

    band = title_band_mm(title_lines, title_pt) + top_extra_mm
    style.fit_margins(fig, pad_mm=pad_mm, cell_mm=letter_cell,
                      reserve_letter=False)
    pos = ax.get_position()
    y0 = pos.y0 if bottom_mm is None else bottom_mm / panel_h_mm
    ax.set_position([pos.x0, y0, pos.width,
                     1.0 - band / panel_h_mm - y0])
    if bottom_mm is not None:
        from cnsfig.layout import recover_x
        recover_x(fig, ax, pad_mm=pad_mm)
        pos = ax.get_position()
    # Centre the title on the plotting box; clamp it inside the canvas. The
    # letter cell bounds it only when the title's top is inside the cell.
    cx = pos.x0 + pos.width / 2
    t = fig.text(cx, 1.0 - (0.4 + top_extra_mm) / panel_h_mm, title,
                 ha="center", va="top", fontsize=title_pt, linespacing=1.15)
    fig.canvas.draw()
    bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
    W = fig.get_figwidth() * fig.dpi
    half = bb.width / 2 / W
    in_cell = (0.4 + top_extra_mm) < letter_cell[1]
    lo = max(letter_cell[0] / panel_w_mm if in_cell else 0.0,
             0.3 / panel_w_mm) + half
    hi = 1.0 - 0.3 / panel_w_mm - half
    if lo > hi:
        raise RuntimeError(f"title {title!r} is wider than the canvas")
    t.set_x(min(max(cx, lo), hi))
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"ink outside the {panel_w_mm} x {panel_h_mm} mm "
                           f"canvas (l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, letter_cell)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    return t


def frame_fixed(fig, ax, *, title, left_mm, right_mm, bottom_mm,
                panel_w_mm, panel_h_mm, ylabel=None, ylabel_pt=None,
                title_pt=None, title_lines=1, letter_cell=None,
                title_top_mm=0.4):
    """Margins in millimetres, the title as a rich run on the canvas centred
    over the plotting frame, the y label as a rich run.

    Fixed margins align the frames of Figure 2 K and L whatever the tick
    labels are.
    """
    import panel_style_cns as style
    from cnsfig.rich import rich_text, rich_ylabel
    title_pt = style.body_pt() if title_pt is None else title_pt
    band = title_band_mm(title_lines, title_pt)
    style.margins_mm(fig, left=left_mm, right=right_mm, top=band,
                     bottom=bottom_mm)
    ax.set_title("")
    if ylabel:
        rich_ylabel(ax, ylabel, fontsize=ylabel_pt, labelpad_pt=1.5)
    pos = ax.get_position()
    cx = pos.x0 + pos.width / 2
    texts = rich_text(fig, cx, 1.0 - title_top_mm / panel_h_mm, title,
                      ha="center", va="top", fontsize=title_pt)
    # A title wider than its frame is centred on it, then kept inside the canvas.
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W = fig.get_figwidth() * fig.dpi
    x0 = min(t.get_window_extent(renderer=r).x0 for t in texts) / W
    x1 = max(t.get_window_extent(renderer=r).x1 for t in texts) / W
    edge = 0.3 / panel_w_mm
    lo = (letter_cell[0] / panel_w_mm if letter_cell else 0.0) + edge
    shift = 0.0
    if x1 > 1.0 - edge:
        shift = (1.0 - edge) - x1
    elif x0 < lo:
        shift = lo - x0
    if shift:
        for t in texts:
            t.remove()
        texts = rich_text(fig, cx + shift, 1.0 - title_top_mm / panel_h_mm,
                          title, ha="center", va="top", fontsize=title_pt)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"ink outside the {panel_w_mm} x {panel_h_mm} mm "
                           f"canvas (l,r,b,t mm): {over}")
    if letter_cell is not None:
        intruders = style.letter_clear(fig, letter_cell)
        if intruders:
            raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    return texts


def draw_boxes(ax, data, positions, colors, *, width=0.6):
    """Draw the standard box plot on `ax`."""
    import panel_style_cns as style
    positions = list(positions)
    if len(positions) != len(data) or len(colors) != len(data):
        raise ValueError("draw_boxes: data, positions and colors must match")
    line = dict(color="black", linewidth=style.RULE_PT)
    bxp = ax.boxplot(
        data, positions=positions, widths=width, patch_artist=True,
        boxprops=dict(linewidth=style.RULE_PT, edgecolor="black"),
        whiskerprops=dict(line), capprops=dict(line), medianprops=dict(line),
        flierprops=dict(marker="o", markerfacecolor="white",
                        markersize=FLIER_MM * PT_PER_MM,
                        markeredgecolor="black",
                        markeredgewidth=style.EDGE_PT))
    for patch, color in zip(bxp["boxes"], colors):
        patch.set_facecolor(color)
    return bxp


def _ink_bottom_mm(text, fontsize):
    """How far above its baseline the string's ink starts, in mm."""
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    ext = TextPath((0, 0), text, size=fontsize,
                   prop=FontProperties(size=fontsize)).get_extents()
    return ext.y0 / PT_PER_MM


# ---------------------------------------------------- the printed P value
# A panel prints the one-tailed test; the legend carries the two-sided value
# (Supplementary Table 6).
# Every panel tests the group it is about against its comparison group, in
# that argument order.
SUBMITTED_ALTERNATIVE = "greater"


def printed_and_two_sided(hi, lo, kind="mannwhitney"):
    """Return (the P the panel prints, the P its legend reports).

    `hi` is the group the panel's directional hypothesis is about and `lo` the
    comparison group; the one-tailed test is hi > lo. Both values come from
    the same data in one call.
    """
    from scipy import stats
    if kind == "mannwhitney":
        test = stats.mannwhitneyu
    elif kind == "wilcoxon":
        test = stats.wilcoxon
    else:
        raise ValueError(f"unknown test {kind!r}")
    printed = float(test(hi, lo, alternative=SUBMITTED_ALTERNATIVE).pvalue)
    two_sided = float(test(hi, lo, alternative="two-sided").pvalue)
    print(f"  figure (one-tailed, {kind}) P = {printed:.4g}; "
          f"legend (two-sided) P = {two_sided:.4g}")
    return printed, two_sided


def bracket(fig, ax, x1, x2, y_max, y_rng, p, *, kind="pair", lift=0.10,
            arm=0.03, ink_gap_mm=INK_GAP_MM, label=None, exact=False):
    """A black bracket over [x1, x2] and its P label.

    The top line sits at y_max + (lift + arm) * y_rng. kind="pair" prints
    p_label's star or value; kind="omnibus" prints "ns" when P >= 0.05.
    exact=True prints p_exact for both. `label` overrides the string and must
    equal f"P = {p:.2f}". The label stays ink_gap_mm above the line whatever
    ylim becomes. Returns (lines, text, top_y), top_y in data units.
    """
    import panel_style_cns as style
    from matplotlib import transforms
    if kind not in ("pair", "omnibus"):
        raise ValueError(f"bracket kind {kind!r}: pair or omnibus")
    y0 = y_max + lift * y_rng
    y_top = y0 + arm * y_rng
    lines = ax.plot([x1, x1, x2, x2], [y0, y_top, y_top, y0],
                    color="black", lw=style.RULE_PT, solid_capstyle="butt")
    if label is not None:
        want = f"P = {p:.2f}"
        if label != want:
            raise ValueError(f"bracket(label={label!r}) does not state its own "
                             f"P value, which prints as {want!r}")
        kw = {"fontsize": style.tick_pt()}
    elif exact:
        label, kw = style.p_exact(p), {"fontsize": style.tick_pt()}
    elif kind == "pair":
        label, kw = style.p_text_kw(p)
    else:
        label = "ns" if p >= 0.05 else style.p_label(p)
        kw = {"fontsize": style.tick_pt()}
    offset_mm = ink_gap_mm - _ink_bottom_mm(label, kw["fontsize"])
    tr = transforms.offset_copy(ax.transData, fig=fig, x=0.0,
                                y=offset_mm * PT_PER_MM, units="points")
    text = ax.text((x1 + x2) / 2, y_top, label, ha="center", va="baseline",
                   color="black", transform=tr, **kw)
    fig.canvas.draw()
    top_y = _ink_top_data(fig, ax, text)
    return lines, text, top_y


def one_vs_others(fig, ax, focal_x, other_xs, y_max, y_rng, p, kw_p, *,
                  comb_lift=0.30, tooth=0.05, rise=0.22, gap_mm=0.4,
                  kw_gap_mm=1.4):
    """One-vs-others annotation: a comb over the other groups, a bracket from
    the focal box to the comb's middle labelled "vs. others, P = x", and the
    Kruskal-Wallis P among the others ("KW P = x") under the comb. Both P
    values are passed in. Returns (label, kw_label) for ylim_above()."""
    import panel_style_cns as style
    from matplotlib import transforms
    lw = style.RULE_PT
    y_comb = y_max + comb_lift * y_rng
    t = tooth * y_rng
    lo, hi = min(other_xs), max(other_xs)
    ax.plot([lo, hi], [y_comb, y_comb], color="black", lw=lw,
            solid_capstyle="butt")
    for x in other_xs:
        ax.plot([x, x], [y_comb, y_comb - t], color="black", lw=lw,
                solid_capstyle="butt")
    mid = (lo + hi) / 2
    y_br = y_comb + rise * y_rng
    y_foot = y_max + 0.04 * y_rng
    # Drawn left to right whichever side the focal box is on.
    xs, ys = [focal_x, focal_x, mid, mid], [y_foot, y_br, y_br, y_comb]
    if focal_x > mid:
        xs, ys = xs[::-1], ys[::-1]
    ax.plot(xs, ys, color="black", lw=lw, solid_capstyle="butt")
    up = transforms.offset_copy(ax.transData, fig=fig, x=0,
                                y=gap_mm * PT_PER_MM, units="points")
    label = ax.text((focal_x + mid) / 2, y_br,
                    f"vs. others, {style.p_exact(p)}", ha="center",
                    va="bottom", color="black", fontsize=style.tick_pt(),
                    transform=up)
    # The KW label sits kw_gap_mm under the tooth ends, clear of the middle tooth.
    down = transforms.offset_copy(ax.transData, fig=fig, x=0,
                                  y=-kw_gap_mm * PT_PER_MM, units="points")
    kw = ax.text(mid, y_comb - t, f"KW P = {kw_p:.2f}", ha="center", va="top",
                 color="black", fontsize=6.0, transform=down)
    return label, kw


def _ink_top_data(fig, ax, text):
    r = fig.canvas.get_renderer()
    bb = text.get_window_extent(renderer=r)
    return ax.transData.inverted().transform((bb.x0, bb.y1))[1]


def ylim_above(ax, text, pad_mm=0.5):
    """Raise the top y limit until `text`'s ink is `pad_mm` inside the axes."""
    fig = ax.figure
    for _ in range(8):
        fig.canvas.draw()
        top = _ink_top_data(fig, ax, text)
        lo, hi = ax.get_ylim()
        per_mm = (hi - lo) / (ax.get_window_extent().height / fig.dpi * 25.4)
        need = top + pad_mm * per_mm
        if need <= hi + 1e-6 * (hi - lo):
            return hi
        ax.set_ylim(lo, need)
    raise RuntimeError("ylim_above did not converge")


def assert_no_points(ax, bxp=None):
    """Raise if `ax` draws individual points: any PathCollection (scatter,
    stripplot) or a marker-only Line2D that is not one of `bxp`'s fliers."""
    from matplotlib.collections import PathCollection
    fliers = set(map(id, bxp["fliers"])) if bxp else set()
    bad = [type(c).__name__ for c in ax.collections
           if isinstance(c, PathCollection)]
    for ln in ax.lines:
        if id(ln) in fliers:
            continue
        if ln.get_marker() not in (None, "None", "", " ") and \
                ln.get_linestyle() in ("None", "none", "", " "):
            bad.append(f"Line2D(marker={ln.get_marker()!r})")
    if bad:
        raise RuntimeError(f"individual points drawn on a box axes: {bad} "
                           "(box panels draw boxes only)")


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        from boxes_selftest import self_test
        sys.exit(self_test())
    sys.exit(0)
