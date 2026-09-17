"""
Axes placed at millimetres, so two panels that call this with the same numbers
print the same.

WHY THIS EXISTS (2026-09-14)
    Figure 5's panels B and C are the same drawing - four UMAPs in a 2 x 2
    grid, each with its own colour bar - drawn by two scripts that both handed
    the layout to `tight_layout`. B's slot is 48.0 x 43.0 mm and C's is
    43.7 x 41.0 mm, so `tight_layout` gave B maps of 12.88 x 14.32 mm and C
    maps of 13.23 x 13.31 mm, side by side on the page across a 1.9 mm
    gutter. Nothing in either script chose a size; the author read the
    mismatch off the page.

    A map here is a square of `map_mm`, its colour bar `cbar_w_mm` wide and
    `cbar_gap_mm` to its right, the columns `col_gap_mm` apart, the rows
    `row_gap_mm` apart with `title_mm` above each map. Every number is a
    millimetre, so the panel prints what the call says and a second panel
    calling with the same numbers prints the same.

    python layout.py --self-test
"""

import sys

PT_PER_MM = 72.0 / 25.4

__all__ = ["umap_grid_mm", "scatter_pair_mm", "corr_annotate", "pin_frame_mm", "recover_x",
           "self_test"]

#: The one geometry of the eight scatter panels of Figure 3 (A, B, D, E) and
#: the two of Figure 5 (K, M), since 2026-09-14 (evening): a square box of
#: SCATTER_BOX_MM, so that "the four are the same size and the axes are
#: symmetric" is a number rather than a wish.
SCATTER_BOX_MM = 13.5
SCATTER_TOP_MM = 3.4         # one title line at 6 pt and its clearance
SCATTER_GAP_MM = 1.5         # between the two boxes of a pair


def scatter_pair_mm(fig, *, left_mm, box_mm=SCATTER_BOX_MM,
                    top_mm=SCATTER_TOP_MM, gap_mm=SCATTER_GAP_MM, n=2):
    """`n` square axes of `box_mm` in a row, at millimetres. Returns them.

    The right axes drop their y tick labels: a pair plots one quantity on one
    scale, and the gutter is `gap_mm` wide for exactly that reason.
    """
    W = fig.get_figwidth() * 25.4
    H = fig.get_figheight() * 25.4
    need_w = left_mm + n * box_mm + (n - 1) * gap_mm
    if need_w > W + 1e-6 or top_mm + box_mm > H + 1e-6:
        raise RuntimeError(f"{n} boxes of {box_mm} mm need {need_w:.1f} x "
                           f"{top_mm + box_mm:.1f} mm; the figure is "
                           f"{W:.1f} x {H:.1f}")
    axes = []
    for i in range(n):
        x = left_mm + i * (box_mm + gap_mm)
        ax = fig.add_axes([x / W, 1 - (top_mm + box_mm) / H, box_mm / W, box_mm / H])
        if i:
            ax.tick_params(axis="y", labelleft=False)
        axes.append(ax)
    return axes


def pin_frame_mm(fig, ax, *, top_mm, bottom_mm):
    """Set the axes' top and bottom edges at millimetres from the canvas
    edges, keeping its x extent. Used after fit_margins by the panels of one
    row so their frames print on one line (sweep_pages.ROW_ALIGN measures
    it); the overflow check that follows says whether the ink still fits."""
    H = fig.get_figheight() * 25.4
    if top_mm + bottom_mm >= H:
        raise ValueError(f"top {top_mm} + bottom {bottom_mm} mm leave no "
                         f"frame in {H:.1f} mm")
    pos = ax.get_position()
    ax.set_position([pos.x0, bottom_mm / H, pos.width,
                     (H - top_mm - bottom_mm) / H])
    recover_x(fig, ax)
    return ax


def recover_x(fig, ax, pad_mm=0.6, rounds=3):
    """Move and narrow `ax` by the measured horizontal overflow.

    A taller frame can grow a tick label ("0.6" where "0.5" was) and push
    ink off the left edge; the frame is then moved and narrowed by the
    measured overflow plus the pad the fit used. Vertical overflow is left
    for the caller's check: it is the pin being wrong.
    """
    import panel_style_cns as style
    W = fig.get_figwidth() * 25.4
    for _ in range(rounds):
        left, right, _b, _t = style.overflow_mm(fig)
        if left <= 0 and right <= 0:
            break
        pos = ax.get_position()
        dx = (left + (pad_mm if left > 0 else 0.0)) / W
        dr = (right + (pad_mm if right > 0 else 0.0)) / W
        ax.set_position([pos.x0 + dx, pos.y0, pos.width - dx - dr, pos.height])
    return ax


def corr_annotate(ax, rho, *, loc="upper left", fontsize=None):
    """Spearman rho inside the box, where the author put it (2026-09-14):
    the P value goes to the legend, through corr_stats.write."""
    import panel_style_cns as style
    fontsize = style.tick_pt() if fontsize is None else fontsize
    x, ha = (0.04, "left") if loc.endswith("left") else (0.96, "right")
    y, va = (0.96, "top") if loc.startswith("upper") else (0.04, "bottom")
    return ax.text(x, y, f"ρ = {rho:.2f}", transform=ax.transAxes, ha=ha,
                   va=va, fontsize=fontsize, zorder=6)


def umap_grid_mm(fig, *, n_rows, n_cols, map_mm, left_mm, top_mm,
                 cbar_w_mm=1.2, cbar_gap_mm=0.6, tick_mm=3.6,
                 col_gap_mm=1.5, row_gap_mm=1.2, title_mm=3.2):
    """Place n_rows x n_cols square map axes with a colour-bar axes each.

    Returns a list of (map_ax, cbar_ax) in row-major order. Raises if the grid
    does not fit the figure; a grid that silently spills is the fault this
    replaces.
    """
    W = fig.get_figwidth() * 25.4
    H = fig.get_figheight() * 25.4
    cell_w = map_mm + cbar_gap_mm + cbar_w_mm + tick_mm
    cell_h = title_mm + map_mm
    need_w = left_mm + n_cols * cell_w + (n_cols - 1) * col_gap_mm
    need_h = top_mm + n_rows * cell_h + (n_rows - 1) * row_gap_mm
    if need_w > W + 1e-6 or need_h > H + 1e-6:
        raise RuntimeError(
            f"a {n_rows} x {n_cols} grid of {map_mm} mm maps needs "
            f"{need_w:.1f} x {need_h:.1f} mm; the figure is {W:.1f} x {H:.1f}")
    out = []
    for r in range(n_rows):
        for c in range(n_cols):
            x = left_mm + c * (cell_w + col_gap_mm)
            y_top = top_mm + r * (cell_h + row_gap_mm) + title_mm
            ax = fig.add_axes([x / W, 1 - (y_top + map_mm) / H,
                               map_mm / W, map_mm / H])
            cb_h = 0.8 * map_mm
            cax = fig.add_axes([(x + map_mm + cbar_gap_mm) / W,
                                1 - (y_top + map_mm - 0.1 * map_mm) / H,
                                cbar_w_mm / W, cb_h / H])
            out.append((ax, cax))
    return out


def self_test() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ok = True

    def control(sign, name, got, want):
        nonlocal ok
        good = got == want
        ok &= good
        print(f"  control {sign}  {name}: "
              f"{'as required' if good else f'WRONG got {got} want {want}'}")

    def grid(w, h, **kw):
        fig = plt.figure(figsize=(w / 25.4, h / 25.4))
        cells = umap_grid_mm(fig, n_rows=2, n_cols=2, map_mm=12.8,
                             left_mm=4.2, top_mm=1.0, **kw)
        sizes = [(round(ax.get_position().width * w, 3),
                  round(ax.get_position().height * h, 3)) for ax, _ in cells]
        plt.close(fig)
        return sizes

    a = grid(43.7, 41.0)
    b = grid(48.0, 43.0)
    control("0", "the same call on two canvases draws maps of one size",
            a == b and len(set(a)) == 1, True)
    control("0", "the map is the millimetre asked for", a[0], (12.8, 12.8))
    try:
        grid(30.0, 41.0)
        caught = False
    except RuntimeError:
        caught = True
    control("+", "MUTANT: a canvas too narrow for the grid raises", caught, True)
    print("\nevery control behaved as required" if ok
          else "\n*** A CONTROL FAILED ***")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(self_test() if "--self-test" in sys.argv else 0)
