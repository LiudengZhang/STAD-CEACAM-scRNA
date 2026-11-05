"""
Labels with gene symbols in italic, digits included.

mathtext's \\it leaves digits upright ("CEACAM" italic, "5" roman), so a
label is drawn as a run of Text artists on one baseline, each part in its own
face. Markup: `*CEACAM5*` is italic, everything else upright. rich_text places
a run anywhere; rich_title, rich_xlabel and rich_ylabel place it where
matplotlib would put the title or axis label.

    python rich.py --self-test
"""

import re
import sys

import matplotlib.transforms as mtransforms

PT_PER_MM = 72.0 / 25.4

__all__ = ["rich_text", "rich_title", "rich_xlabel", "rich_ylabel",
           "parse"]

_PART = re.compile(r"(\*[^*]+\*)")


def parse(markup):
    """[(text, italic)] for one line of markup. `*ABC1*` is italic."""
    parts = []
    for piece in _PART.split(markup):
        if not piece:
            continue
        if piece.startswith("*") and piece.endswith("*") and len(piece) > 2:
            parts.append((piece[1:-1], True))
        else:
            if "*" in piece:
                raise ValueError(f"unbalanced '*' in {markup!r}")
            parts.append((piece, False))
    return parts


def _renderer(fig):
    return fig.canvas.get_renderer()


_FONT_CACHE = {}


def _advance_pt(text, fontsize, italic):
    """The advance width of `text` in points, from the face's own metrics."""
    import matplotlib
    from matplotlib import font_manager, ft2font
    key = (italic, matplotlib.rcParams["font.family"][0]
           if isinstance(matplotlib.rcParams["font.family"], list)
           else matplotlib.rcParams["font.family"])
    if key not in _FONT_CACHE:
        prop = font_manager.FontProperties(
            family=matplotlib.rcParams["font.family"],
            style="italic" if italic else "normal")
        _FONT_CACHE[key] = ft2font.FT2Font(font_manager.findfont(prop))
    font = _FONT_CACHE[key]
    font.set_size(fontsize, 72)
    total = 0.0
    for ch in text:
        flags = getattr(ft2font, "LoadFlags", None)
        no_hint = flags.NO_HINTING if flags is not None else ft2font.LOAD_NO_HINTING
        glyph = font.load_char(ord(ch), flags=no_hint)
        total += glyph.linearHoriAdvance / 65536.0
    return total


def rich_text(owner, x, y, markup, *, transform=None, ha="center",
              va="baseline", rotation=0, fontsize=None, linespacing=1.2,
              color="black", zorder=None):
    """Draw `markup` at (x, y) in `transform`, parts in their own faces."""
    import matplotlib.pyplot as plt  # noqa: F401  (backend must be up)
    fig = owner.figure if hasattr(owner, "figure") and owner.figure is not None \
        else owner
    if owner is fig:
        make = fig.text
    else:
        make = owner.text
    if transform is None:
        transform = fig.transFigure if owner is fig else owner.transAxes
    if rotation not in (0, 90):
        raise ValueError("rich_text supports rotation 0 or 90")
    if fontsize is None:
        import matplotlib
        fontsize = matplotlib.rcParams["font.size"]
    r = _renderer(fig)
    dpi = fig.dpi
    px_per_pt = dpi / 72.0

    lines = [parse(ln) for ln in markup.split("\n")]
    pitch_pt = linespacing * fontsize

    # Leading spaces: a Text is drawn from its first glyph, so a part such
    # as " (Epi)" is shifted right by the advance of its leading spaces.
    space_pt = _advance_pt(" ", fontsize, False)

    # Measure every part at the anchor, with va='baseline', so widths and the
    # ascent/descent of each line are known before anything is placed.
    made = []          # (line_index, Text, width_pt)
    line_metrics = []  # (width_pt, ascent_pt, descent_pt)
    for li, parts in enumerate(lines):
        widths, asc, desc = [], 0.0, 0.0
        texts = []
        for text, italic in parts:
            # Parts are drawn without end spaces (Agg and PDF treat a leading
            # space differently); the spaces are added to the cursor instead.
            n_lead = len(text) - len(text.lstrip(" "))
            n_trail = len(text) - len(text.rstrip(" "))
            stripped = text.strip(" ")
            t = make(x, y, stripped, transform=transform, ha="left", va="baseline",
                     rotation=rotation, rotation_mode="anchor",
                     fontsize=fontsize, fontstyle="italic" if italic else "normal",
                     color=color, zorder=zorder)
            bb = t.get_window_extent(renderer=r)
            ax0, ay0 = transform.transform((x, y))
            if rotation == 0:
                asc = max(asc, (bb.y1 - ay0) / px_per_pt)
                desc = max(desc, (ay0 - bb.y0) / px_per_pt)
            else:
                asc = max(asc, (ax0 - bb.x0) / px_per_pt)
                desc = max(desc, (bb.x1 - ax0) / px_per_pt)
            widths.append(_advance_pt(stripped, fontsize, italic)
                          + (n_lead + n_trail) * space_pt)
            t._rich_lead_pt = n_lead * space_pt
            texts.append(t)
        line_metrics.append((sum(widths), asc, desc))
        for t, w in zip(texts, widths):
            made.append((li, t, w))

    n = len(lines)
    block_asc = line_metrics[0][1]
    block_desc = line_metrics[-1][2] + (n - 1) * pitch_pt
    # Baseline of line i relative to the anchor, along the "down" direction.
    if va == "baseline":
        first_baseline = 0.0
    elif va == "top":
        first_baseline = -block_asc
    elif va == "bottom":
        first_baseline = block_desc
    elif va in ("center", "center_baseline"):
        first_baseline = (block_desc - block_asc) / 2.0
    else:
        raise ValueError(f"va {va!r}")

    for li in range(n):
        width, _, _ = line_metrics[li]
        along = {"left": 0.0, "center": -width / 2.0, "right": -width}[ha]
        down = first_baseline - li * pitch_pt      # points, +up
        cursor = along
        for lj, t, w in made:
            if lj != li:
                continue
            lead = t._rich_lead_pt
            if rotation == 0:
                dx, dy = cursor + lead, down
            else:
                # Reading direction is +y; "down" (across) is +x.
                dx, dy = -down, cursor + lead
            t.set_transform(mtransforms.offset_copy(
                transform, fig=fig, x=dx, y=dy, units="points"))
            cursor += w
    return [t for _, t, _ in made]


def _tick_extent_pt(ax, axis):
    """How far the tick marks and labels of `axis` reach outside the axes
    box, in points (0 when the labels are off)."""
    fig = ax.figure
    fig.canvas.draw()          # tick labels are laid out by the draw
    r = _renderer(fig)
    px_per_pt = fig.dpi / 72.0
    box = ax.get_window_extent(renderer=r)
    reach = 0.0
    for tick in axis.get_major_ticks():
        for lab in (tick.label1, tick.label2):
            if not lab.get_visible() or not lab.get_text():
                continue
            bb = lab.get_window_extent(renderer=r)
            if bb.width <= 0:
                continue
            if axis is ax.xaxis:
                reach = max(reach, (box.y0 - bb.y0) / px_per_pt)
            else:
                reach = max(reach, (box.x0 - bb.x0) / px_per_pt)
    # Tick marks alone, when there are no labels.
    length = axis.get_major_ticks()[0]._size if axis.get_major_ticks() else 0.0
    return max(reach, length)


def rich_title(ax, markup, *, fontsize=None, pad_pt=None, linespacing=1.2):
    """The centred axes title, as a rich run; the plain title is cleared."""
    import matplotlib
    if fontsize is None:
        fontsize = matplotlib.rcParams["axes.titlesize"]
    if pad_pt is None:
        pad_pt = matplotlib.rcParams["axes.titlepad"]
    ax.set_title("")
    anchor = mtransforms.offset_copy(ax.transAxes, fig=ax.figure, x=0,
                                     y=pad_pt, units="points")
    return rich_text(ax, 0.5, 1.0, markup, transform=anchor, ha="center",
                     va="bottom", fontsize=fontsize, linespacing=linespacing)


def rich_xlabel(ax, markup, *, fontsize=None, labelpad_pt=None,
                linespacing=1.2):
    """The x axis label, as a rich run, below the tick labels."""
    import matplotlib
    if fontsize is None:
        fontsize = matplotlib.rcParams["axes.labelsize"]
    if labelpad_pt is None:
        labelpad_pt = matplotlib.rcParams["axes.labelpad"]
    ax.set_xlabel("")
    down = _tick_extent_pt(ax, ax.xaxis) + labelpad_pt
    anchor = mtransforms.offset_copy(ax.transAxes, fig=ax.figure, x=0,
                                     y=-down, units="points")
    return rich_text(ax, 0.5, 0.0, markup, transform=anchor, ha="center",
                     va="top", fontsize=fontsize, linespacing=linespacing)


def rich_ylabel(ax, markup, *, fontsize=None, labelpad_pt=None,
                linespacing=1.3, y=0.5):
    """The y axis label, as a rich run rotated 90, left of the tick labels."""
    import matplotlib
    if fontsize is None:
        fontsize = matplotlib.rcParams["axes.labelsize"]
    if labelpad_pt is None:
        labelpad_pt = matplotlib.rcParams["axes.labelpad"]
    ax.set_ylabel("")
    left = _tick_extent_pt(ax, ax.yaxis) + labelpad_pt
    anchor = mtransforms.offset_copy(ax.transAxes, fig=ax.figure, x=-left,
                                     y=0, units="points")
    # Rotated 90, "top" of the block faces left, so va='bottom' puts the
    # block's right edge (its last line's descent) at the anchor.
    return rich_text(ax, 0.0, y, markup, transform=anchor, ha="center",
                     va="bottom", rotation=90, fontsize=fontsize,
                     linespacing=linespacing)


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        from rich_selftest import self_test
        sys.exit(self_test())
    sys.exit(0)
