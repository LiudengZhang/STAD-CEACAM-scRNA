"""
The box panels that print as one family, drawn by one function.

WHY THIS EXISTS (2026-09-14, evening)
    Figure 5's panel J is four two-group boxes and panel L is a fifth; the
    author read them side by side and asked why five drawings of the same
    thing were five scripts' worth of layout. And on 3H-J, 5J and 5L the
    boxes stood flush against the y axis: `xlim(-0.25, 1.25)` around boxes
    0.5 wide at 0 and 1 leaves no paper between the first box and the spine.

    Two things live here, so they have one owner:

    box_xlim(positions, width, clear)
        The x limits that leave `clear` box-widths of paper between the
        outermost boxes and the spines. Every boxplot and violin in Figures
        2-5 sets its limits through it. The page-level check
        (`check_box_clearance` in 10_Reproduction/check_restyled_panel.py)
        reads the drawn bodies and the drawn spines back off the page and
        holds them a millimetre apart, with a mutation.

    finish_two_group(...)
        The frame of an R-versus-NR box: the tick labels with their counts,
        the y label as a rich run at the tick size, the title drawn on the
        canvas in a band of fixed height above the plotting box, the margins
        fitted with the letter corner clear. The five callers compute their
        own values, draw their own boxes and points (their own RNG, their own
        bracket geometry, so not one vertex moves), then call this.

ONE BOX, ONE BRACKET, NO POINTS (2026-09-16, the author's fifth reading)
    Thirty box panels across Figures 2, 4, 5 and S2-S10 were drawn by
    twenty-eight scripts, each with its own flier marker (0.79 mm open
    circles on 2H/I, 2.5 mm filled discs on 4C/4H, none on 5J), its own
    jitter of the individual samples (seeded 0, 42, 42+i, or not drawn) and
    its own bracket arithmetic. The author read three defects off the pages
    and ruled on all three:

      "Figure 4's outliers are too large"        -> one flier: draw_boxes()
      "the star sits too far above the line"     -> one bracket: bracket()
      "box plots should look alike throughout;
       don't show every point"                   -> no points: assert_no_points()

    The star floated because every script placed its P text by BASELINE at
    `bracket_y + 0.05 * range`; a 9 pt asterisk's ink starts 1.4 mm above its
    baseline (Nimbus Sans, measured with TextPath), so the same offset that
    puts "P = 0.30" 0.4 mm above the line puts "*" 1.8 mm above it, and on a
    30 mm axes 0.05 of the range is another millimetre. bracket() measures
    the ink and puts ITS bottom edge `ink_gap_mm` above the bracket, as a
    display offset from the bracket's data point, so a later change of ylim
    does not move it. Removing the individual points is the author's ruling
    of 2026-09-16 ("接受"), a declared departure from the published panels,
    recorded per row in PROVENANCE; the comparator shows it as removed
    PathCollections and nothing else.

    draw_boxes(ax, data, positions, colors, width)
        The Figure 2 H/I box: RULE_PT black lines, EDGE_PT flier edge, 0.79 mm
        open-circle fliers, faces from `colors`. Data is passed straight to
        ax.boxplot; nothing is computed here.
    bracket(fig, ax, x1, x2, y_max, y_rng, p, kind)
        The black RULE_PT bracket and its label, ink-gapped. kind="pair" sets
        p_text_kw (a STAR_PT star or the two-decimal value); kind="omnibus"
        prints "ns" at the tick size when P >= 0.05 (Kruskal-Wallis over the
        three groups in 2H/I and 5D/E).
    ylim_above(ax, text, pad_mm)
        Raises the top limit until the label's ink is inside the axes.
    assert_no_points(ax, bxp)
        Refuses an axes that carries a scatter or a marker-only line other
        than the boxes' fliers.

    python boxes.py --self-test
"""

import sys

_CONFIG = __import__("pathlib").Path(__file__).resolve().parent.parent
if str(_CONFIG) not in sys.path:
    sys.path.insert(0, str(_CONFIG))
from panel_style_cns import PT_PER_MM  # noqa: E402

__all__ = ["box_xlim", "finish_two_group", "frame_fixed", "title_band_mm", "CLEAR",
           "draw_boxes", "bracket", "ylim_above", "assert_no_points", "FLIER_MM",
           "INK_GAP_MM", "self_test"]

#: The open-circle flier, millimetres across: counted off the published
#: Figure 2 H/I and used by every box since 2026-09-16.
FLIER_MM = 0.79
#: Paper between a bracket's top line and the bottom of its label's ink.
INK_GAP_MM = 0.4

#: Paper between the outermost box and the spine, in box widths. 0.6 of a
#: 0.5-wide box at a 1.0 pitch is 0.3 data units; on an 8 mm plotting box
#: that is about 1.2 mm, above the 1.0 mm the page check demands.
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
                     top_extra_mm=0.0, bottom_mm=None):
    """Frame an R/NR box axes and draw its title on the canvas.

    The title band is `title_lines` lines tall whatever the title's own line
    count, so a row of these panels prints its plotting boxes at one height.
    The title is centred over the plotting box and kept inside the canvas.
    """
    import panel_style_cns as style
    from cnsfig.rich import rich_ylabel
    tick_pt = style.tick_pt()
    ylabel_pt = tick_pt if ylabel_pt is None else ylabel_pt
    title_pt = style.body_pt() if title_pt is None else title_pt

    ax.set_xticks(list(positions))
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

    For the box pairs that print in two columns (Figure 2 K over L): every
    box of the family names the same left and bottom margins, so their
    frames print at the same x and the pairs' frames at the same y, whatever
    their tick labels are.
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
    # A title wider than its frame is centred on the frame and then kept
    # inside the canvas (the right box of Figure 2 L: 16 mm over 13 mm).
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
    """The one box of the paper (Figure 2 H/I's), on `ax`.

    `data` is a list of 1-D arrays handed to ax.boxplot untouched; `colors`
    one face colour per box. Returns the bxp dict so the fliers can be
    excused by assert_no_points().
    """
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
    """How far above its baseline the string's ink starts, in mm (TextPath;
    1.4 mm for a 9 pt asterisk, ~0 for letters and digits)."""
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    ext = TextPath((0, 0), text, size=fontsize,
                   prop=FontProperties(size=fontsize)).get_extents()
    return ext.y0 / PT_PER_MM


def bracket(fig, ax, x1, x2, y_max, y_rng, p, *, kind="pair", lift=0.10,
            arm=0.03, ink_gap_mm=INK_GAP_MM):
    """A black RULE_PT bracket over [x1, x2] and its P label, ink-gapped.

    The top line sits at y_max + (lift + arm) * y_rng, the arms reach down
    `arm * y_rng` (Figure 2 H/I's shape). kind="pair": p_text_kw (star at
    STAR_PT or the two-decimal value at the tick size). kind="omnibus": "ns"
    at the tick size when P >= 0.05, else p_label. The label's baseline is a
    DISPLAY offset from the top line's data point (transforms.offset_copy),
    so the ink stays `ink_gap_mm` above the line whatever ylim becomes.
    Returns (lines, text, top_y): top_y is the top of the label's ink in data
    units, for callers that set ylim by hand; ylim_above() does it for them.
    """
    import panel_style_cns as style
    from matplotlib import transforms
    if kind not in ("pair", "omnibus"):
        raise ValueError(f"bracket kind {kind!r}: pair or omnibus")
    y0 = y_max + lift * y_rng
    y_top = y0 + arm * y_rng
    lines = ax.plot([x1, x1, x2, x2], [y0, y_top, y_top, y0],
                    color="black", lw=style.RULE_PT, solid_capstyle="butt")
    if kind == "pair":
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


def _ink_top_data(fig, ax, text):
    r = fig.canvas.get_renderer()
    bb = text.get_window_extent(renderer=r)
    return ax.transData.inverted().transform((bb.x0, bb.y1))[1]


def ylim_above(ax, text, pad_mm=0.5):
    """Raise the top y limit until `text`'s ink is `pad_mm` inside the axes.

    Iterates, because raising the limit rescales the axes and the label is
    pinned to the bracket by a display offset.
    """
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
                           "(the author's ruling of 2026-09-16: boxes only)")


def measured_ink_gap_mm(fig, ax, text, line):
    """The paper between the bracket's top line and the label's ink, in mm,
    read off a rendered bitmap: the first dark row above the line inside
    the label's x span. What the self-test trusts; TextPath is what it
    checks."""
    import numpy as np
    fig.canvas.draw()
    dpi = fig.dpi
    buf = np.asarray(fig.canvas.buffer_rgba())[..., :3].min(axis=2)
    H = buf.shape[0]
    bb = text.get_window_extent(renderer=fig.canvas.get_renderer())
    xs = np.array(line.get_xdata()[1:3]); ys = np.array(line.get_ydata()[1:3])
    (px0, py), (px1, _) = ax.transData.transform(np.column_stack([xs, ys]))
    row_line = int(round(H - py))
    c0, c1 = int(bb.x0) + 1, int(bb.x1) - 1
    dark = buf[:, c0:c1] < 128
    # walk up from just above the line's own stroke (RULE_PT ~ 2 px at 300 dpi)
    r = row_line - 3
    while r > 0 and not dark[r].any():
        r -= 1
    if r <= 0:
        raise RuntimeError("no ink found above the bracket")
    return (row_line - r) / dpi * 25.4


def self_test():
    """The limits leave paper; the naive limits do not (mutation); the frame
    puts the plotting box below a band of the declared height; the one box
    and its bracket put the label's ink 0.4 mm above the line, and the old
    baseline placement does not (mutation); points are refused (mutation)."""
    import matplotlib
    matplotlib.use("Agg")
    import numpy as np
    import panel_style_cns as style

    ok = True

    def control(sign, name, good):
        nonlocal ok
        ok &= bool(good)
        print(f"  control {sign}  {name}: {'as required' if good else 'WRONG'}")

    lo, hi = box_xlim([0, 1], 0.5)
    control("0", "box_xlim leaves 0.3 units each side", abs(lo + 0.55) < 1e-9
            and abs(hi - 1.55) < 1e-9)
    # MUTATION: the limits the scripts used to set leave 0.0 units.
    naive = (-0.25, 1.25)
    gap_naive = (0 - 0.25) - naive[0]
    control("+", "the old xlim(-0.25, 1.25) leaves no paper (the defect)",
            abs(gap_naive) < 1e-9)
    try:
        box_xlim([], 0.5)
        control("+", "empty positions raise", False)
    except ValueError:
        control("+", "empty positions raise", True)

    style.apply(title_fontsize=7)
    fig, ax = style.subplots_mm(20, 32)
    rng = np.random.default_rng(0)
    ax.boxplot([rng.normal(size=6), rng.normal(size=6)], positions=[0, 1],
               widths=0.5, patch_artist=True)
    finish_two_group(fig, ax, title="Monocytes/\nMacrophages",
                     tick_labels=["R\n(n=5)", "NR\n(n=6)"], positions=[0, 1],
                     width=0.5, ylabel_markup="PD-L1 (*CD274*)\nexpression",
                     letter_cell=(3.5, 3.6), panel_w_mm=20, panel_h_mm=32)
    pos = ax.get_position()
    top_mm = (1 - pos.y1) * 32
    control("0", "plotting box starts below the two-line band",
            top_mm >= title_band_mm(2, style.body_pt()) - 0.05)
    control("0", "x limits are the clearance limits",
            tuple(round(v, 3) for v in ax.get_xlim()) == (-0.55, 1.55))
    fig2, ax2 = style.subplots_mm(20, 32)
    ax2.boxplot([rng.normal(size=6)], positions=[0], widths=0.5)
    finish_two_group(fig2, ax2, title="Epithelial", tick_labels=["R"],
                     positions=[0], width=0.5, ylabel_markup="x",
                     letter_cell=(3.5, 3.6), panel_w_mm=20, panel_h_mm=32)
    control("0", "a one-line title takes the same band as a two-line one",
            abs((1 - ax2.get_position().y1) * 32 - top_mm) < 0.3)

    # ONE BOX, ONE BRACKET, NO POINTS (2026-09-16). Rendered at 300 dpi so a
    # 0.1 mm error is a pixel.
    from matplotlib.lines import Line2D
    data = [rng.normal(size=6), rng.normal(size=6) + 1.0]
    data[0][0] = 4.0                                     # one flier
    y_all = np.concatenate(data)
    y_max, y_rng = y_all.max(), y_all.max() - y_all.min()
    lo_ok, hi_ok = 0.25, 0.55
    for label, p, kind in (("star", 0.01, "pair"), ("value", 0.30, "pair"),
                           ("ns", 0.92, "omnibus")):
        fig3, ax3 = style.subplots_mm(20, 32)
        fig3.set_dpi(300)
        bxp = draw_boxes(ax3, data, [0, 1], ["#bde0fe", "#ffcfd2"], width=0.5)
        lines, text, top_y = bracket(fig3, ax3, 0, 1, y_max, y_rng, p, kind=kind)
        ax3.set_xlim(*box_xlim([0, 1], 0.5))
        ylim_above(ax3, text)
        style.fit_margins(fig3, pad_mm=0.6, cell_mm=(3.5, 3.6))
        gap = measured_ink_gap_mm(fig3, ax3, text, lines[0])
        control("0", f"{label} ink {gap:.2f} mm above the bracket "
                     f"(window {lo_ok}-{hi_ok})", lo_ok <= gap <= hi_ok)
        control("0", f"{label} label inside the axes after ylim_above",
                ax3.get_ylim()[1] >= _ink_top_data(fig3, ax3, text) - 1e-9)
        if label == "star":
            control("0", "the star is set at STAR_PT",
                    text.get_fontsize() == style.STAR_PT)
            fl = bxp["fliers"][0]
            control("0", "flier is 0.79 mm across, white, EDGE_PT edge",
                    abs(fl.get_markersize() - FLIER_MM * PT_PER_MM) < 1e-6
                    and fl.get_markerfacecolor() == "white"
                    and fl.get_markeredgewidth() == style.EDGE_PT)
            control("0", "box, whisker and median lines are RULE_PT",
                    all(a.get_linewidth() == style.RULE_PT for a in
                        bxp["boxes"] + bxp["whiskers"] + bxp["medians"]
                        + bxp["caps"]))
            try:
                assert_no_points(ax3, bxp)
                control("0", "a boxes-only axes passes assert_no_points", True)
            except RuntimeError:
                control("0", "a boxes-only axes passes assert_no_points", False)
            # MUTATION: the star placed by baseline at y + 0.05 * range, as
            # every script did until today, floats well outside the window.
            fig4, ax4 = style.subplots_mm(20, 32)
            fig4.set_dpi(300)
            draw_boxes(ax4, data, [0, 1], ["#bde0fe", "#ffcfd2"], width=0.5)
            yb = y_max + 0.10 * y_rng
            old_line = ax4.plot([0, 0, 1, 1], [yb, yb + 0.03 * y_rng,
                                               yb + 0.03 * y_rng, yb],
                                "k-", lw=style.RULE_PT)[0]
            s, kw = style.p_text_kw(p)
            old_text = ax4.text(0.5, yb + 0.05 * y_rng, s, ha="center", **kw)
            ax4.set_xlim(*box_xlim([0, 1], 0.5))
            ax4.set_ylim(top=y_max + 0.55 * y_rng)
            style.fit_margins(fig4, pad_mm=0.6, cell_mm=(3.5, 3.6))
            old_gap = measured_ink_gap_mm(fig4, ax4, old_text, old_line)
            control("+", f"the old baseline placement floats the star "
                         f"{old_gap:.2f} mm off (convicted)",
                    not (lo_ok <= old_gap <= hi_ok))
            # MUTATION: a strip point on the axes is refused.
            ax3.scatter([0.0], [float(np.median(data[0]))], s=4)
            try:
                assert_no_points(ax3, bxp)
                control("+", "a scatter point on the box axes is refused", False)
            except RuntimeError:
                control("+", "a scatter point on the box axes is refused", True)
            ax3.collections[-1].remove()
            ax3.add_line(Line2D([1.0], [1.0], marker="o", linestyle="None"))
            try:
                assert_no_points(ax3, bxp)
                control("+", "a marker-only line on the box axes is refused", False)
            except RuntimeError:
                control("+", "a marker-only line on the box axes is refused", True)
    print("\nevery control behaved as required" if ok
          else "\n*** A CONTROL FAILED ***")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(self_test() if "--self-test" in sys.argv else 0)
