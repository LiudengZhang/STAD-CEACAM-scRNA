"""
Axes placed in millimetres, so panels called with the same numbers print at
the same size.

umap_grid_mm places a grid of square maps, each with a colour bar to its
right; scatter_pair_mm places a row of square scatter boxes; pin_frame_mm
aligns frames within a row; umap_corner_axis draws the UMAP1/UMAP2 arrows.

    python layout.py --self-test
"""

import sys

PT_PER_MM = 72.0 / 25.4

__all__ = ["umap_grid_mm", "scatter_pair_mm", "corr_annotate", "pin_frame_mm",
           "umap_corner_axis",
           "recover_x"]

# Square box size shared by the scatter panels of Figure 3 (A, B, D, E) and
# Figure 5 (K, M).
SCATTER_BOX_MM = 13.5
SCATTER_TOP_MM = 3.4         # one title line at 6 pt and its clearance
SCATTER_GAP_MM = 1.5         # between the two boxes of a pair


def scatter_pair_mm(fig, *, left_mm, box_mm=SCATTER_BOX_MM,
                    top_mm=SCATTER_TOP_MM, gap_mm=SCATTER_GAP_MM, n=2):
    """`n` square axes of `box_mm` in a row, at millimetres. Returns them."""
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
    edges, keeping its x extent. Called after fit_margins so the panels of
    one row print their frames on one line."""
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
    """Move and narrow `ax` by the measured horizontal overflow."""
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


def corr_annotate(ax, rho, *, p=None, loc="upper left", fontsize=None):
    """Spearman rho inside the box. With `p`, the P value under it, in the
    legend's form (corr_stats.p_string), so panel and legend match."""
    import panel_style_cns as style
    from cnsfig import corr_stats
    fontsize = style.tick_pt() if fontsize is None else fontsize
    x, ha = (0.04, "left") if loc.endswith("left") else (0.96, "right")
    y, va = (0.96, "top") if loc.startswith("upper") else (0.04, "bottom")
    # Use the true minus sign (U+2212), as the axes do, not a hyphen.
    r = f"{rho:.2f}".replace("-", "\u2212")
    text = f"ρ = {r}" if p is None else f"ρ = {r}\n{corr_stats.p_string(p)}"
    return ax.text(x, y, text, transform=ax.transAxes, ha=ha,
                   va=va, fontsize=fontsize, zorder=6)


def umap_grid_mm(fig, *, n_rows, n_cols, map_mm, left_mm, top_mm,
                 cbar_w_mm=1.2, cbar_gap_mm=0.6, tick_mm=3.6,
                 col_gap_mm=1.5, row_gap_mm=1.2, title_mm=3.2):
    """Place n_rows x n_cols square map axes with a colour-bar axes each."""
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


def umap_corner_axis(ax, *, inset_mm=(3.2, 3.2), length_mm=6.0,
                     names=("UMAP1", "UMAP2"), gap_mm=0.8):
    """Draw UMAP1/UMAP2 corner arrows from the lower-left of a UMAP axes."""
    import matplotlib.transforms as mtransforms
    from matplotlib.patches import FancyArrowPatch, ArrowStyle
    import panel_style_cns as style

    fig = ax.figure
    # millimetres -> display, anchored at the axes' (0, 0); follows the axes
    # if the margins are refitted after this call.
    mm = (mtransforms.Affine2D().scale(1 / 25.4) + fig.dpi_scale_trans
          + mtransforms.ScaledTranslation(0, 0, ax.transAxes))
    x0, y0 = inset_mm
    head = ArrowStyle("-|>", head_length=0.5, head_width=0.22)
    arrows = []
    for end in ((x0 + length_mm, y0), (x0, y0 + length_mm)):
        a = FancyArrowPatch((x0, y0), end, transform=mm, arrowstyle=head,
                            mutation_scale=style.tick_pt(), shrinkA=0,
                            shrinkB=0, lw=style.RULE_PT, color="black",
                            joinstyle="miter", capstyle="butt", clip_on=False,
                            zorder=6)
        ax.add_patch(a)
        arrows.append(a)

    # The axes' own labels, shown: the axis is switched on for them alone.
    ax.set_axis_on()
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(which="both", bottom=False, left=False, top=False,
                   right=False, labelbottom=False, labelleft=False,
                   labeltop=False, labelright=False)
    ax.set_xlabel(names[0])
    ax.set_ylabel(names[1])
    fs = style.tick_pt()
    ax.xaxis.set_label_coords(x0, y0 - gap_mm, transform=mm)
    ax.xaxis.label.set(visible=True, ha="left", va="top", rotation=0,
                       fontsize=fs, fontweight="normal", color="black")
    ax.yaxis.set_label_coords(x0 - gap_mm, y0, transform=mm)
    # rotation_mode "anchor": the alignment is taken in the unrotated frame,
    # so the name starts at the corner and reads upward, left of the arrow.
    ax.yaxis.label.set(visible=True, ha="left", va="bottom", rotation=90,
                       rotation_mode="anchor", fontsize=fs,
                       fontweight="normal", color="black")
    return arrows


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        from layout_selftest import self_test
        sys.exit(self_test())
    sys.exit(0)
