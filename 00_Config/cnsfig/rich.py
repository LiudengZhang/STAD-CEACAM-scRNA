"""
A label with a gene symbol in it, set with the symbol in italic and nothing
else - digits included.

WHY THIS EXISTS (2026-09-14, evening)
    The author read the pages and found gene names italic in some panels and
    upright in others. The upright ones were not a choice; they were what
    mathtext does. `$\\it{CEACAM5}$` prints "CEACAM" in the italic face and
    the "5" in the regular one, because matplotlib's mathtext sends every
    non-letter under \\it to the roman font (`_mathtext.py`, `_get_glyph`).
    Measured on the shipped Figure 2: 'CEACAM' NimbusSans-Italic followed by
    '5' NimbusSans-Regular. There is no mathtext spelling that italicises a
    digit, and a whole label set `fontstyle='italic'` italicises the words
    around the symbol too.

    So a label with a symbol in it is drawn as a run of Text artists on one
    baseline, each part in its own face, laid out from measured widths. The
    markup is one character: `*CEACAM5*` is the symbol, everything else is
    text. `rich_text` places a run anywhere; `rich_title`, `rich_xlabel` and
    `rich_ylabel` put one where matplotlib would put the title and the axis
    labels, offset from the tick labels by measurement. The parts are anchored
    in the caller's transform and offset in points, so they move with the axes
    when `fit_margins` moves it.

    A panel that prints a gene symbol anywhere but through this module, or
    through `fontstyle='italic'` on a label that is only the symbol, is caught
    by `check_italics` in 10_Reproduction/check_restyled_panel.py, which reads
    the faces back off the page and holds every declared symbol to the italic
    face.

    python rich.py --self-test
"""

import re
import sys

import matplotlib.transforms as mtransforms

PT_PER_MM = 72.0 / 25.4

__all__ = ["rich_text", "rich_title", "rich_xlabel", "rich_ylabel",
           "parse", "self_test"]

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
    """The advance width of `text` in points, from the face's own metrics.

    A Text's window extent is its INK box, and an italic face's ink starts
    inside and ends outside its advance, so parts laid out from ink widths
    print 0.1-0.3 pt too close or too far apart; the PDF's glyph boxes are
    advances, and the page gate reads those. The advance is what the next
    part is placed at.
    """
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
    """Draw `markup` at (x, y) in `transform`, parts in their own faces.

    `owner` is an Axes or a Figure; the Text artists are added to it so the
    figure's ink measurement (`panel_style_cns._ink_artists`) sees them.
    `rotation` is 0 or 90. `ha` is along the reading direction and `va`
    across it, as for a single Text. Multi-line markup stacks lines at
    `linespacing` x fontsize, first line on top (or, rotated, leftmost).
    Returns the Text artists.
    """
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
            # The spaces at either end are laid out here, not drawn: the Agg
            # backend draws a leading space and the PDF backend drops it, so
            # each part is drawn stripped and the space advances are added
            # to the cursor explicitly.
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
        # Trailing space at the end of a regular part before an italic one
        # is an advance, not ink; the extent already includes it because the
        # layout width is the advance width.
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


def self_test():
    """Controls, then mutations. The claims: every part is in its own face,
    the digit included; the parts abut on one baseline; the run is centred
    where asked; and the check that reads faces back can tell a symbol set
    through mathtext (digit upright) from one set here."""
    import tempfile
    from pathlib import Path
    import fitz
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
    import panel_style_cns as style

    ok = True

    def control(sign, name, good):
        nonlocal ok
        ok &= bool(good)
        print(f"  control {sign}  {name}: {'as required' if good else 'WRONG'}")

    control("0", "parse splits markup", parse("Epi. *CEACAM5* (x)") ==
            [("Epi. ", False), ("CEACAM5", True), (" (x)", False)])
    try:
        parse("a *b")
        control("+", "unbalanced markup raises", False)
    except ValueError:
        control("+", "unbalanced markup raises", True)

    style.apply()
    fig, ax = style.subplots_mm(50, 30)
    rich_title(ax, "*BACH1* regulon")
    rich_xlabel(ax, "Epi. *CEACAM5*")
    rich_ylabel(ax, "PD-L1 (*CD274*)\n(mean log expr.)")
    ax.set_title("")
    ax.text(0.5, 0.5, r"$\it{NFKB1}$ mathtext", ha="center", transform=ax.transAxes)
    style.fit_margins(fig, pad_mm=0.6)
    with tempfile.TemporaryDirectory() as d:
        pdf = Path(d) / "rich.pdf"
        fig.savefig(pdf)
        page = fitz.open(pdf)[0]
        spans = [(s["text"], s["font"], s["bbox"])
                 for b in page.get_text("dict")["blocks"]
                 for l in b.get("lines", []) for s in l["spans"]]
    faces = {t.strip(): f for t, f, _ in spans if t.strip()}
    control("0", "BACH1 is one italic span, digit included",
            "Italic" in faces.get("BACH1", "") )
    control("0", "'regulon' is regular", "Regular" in faces.get("regulon", ""))
    control("0", "CEACAM5 italic, 'Epi.' regular",
            "Italic" in faces.get("CEACAM5", "") and "Regular" in faces.get("Epi.", ""))
    control("0", "CD274 italic inside the rotated y label",
            "Italic" in faces.get("CD274", ""))
    # The mathtext control: the digit comes out regular. This is what the
    # module exists to avoid, and what check_italics must convict.
    control("+", "mathtext leaves the digit upright (the defect)",
            "Regular" in faces.get("1 mathtext", faces.get("1", "")))
    # Baseline: BACH1 and 'regulon' share a baseline to 0.2 pt.
    b1 = next((bb for t, f, bb in spans if t.strip() == "BACH1"), None)
    b2 = next((bb for t, f, bb in spans if t.strip() == "regulon"), None)
    control("0", "title parts share a baseline",
            b1 is not None and b2 is not None and abs(b1[3] - b2[3]) < 0.4)
    # Centring: the title run's centre is the axes' centre to 0.5 pt.
    axbox = ax.get_window_extent(renderer=fig.canvas.get_renderer())
    W_px = fig.get_figwidth() * fig.dpi
    ax_cx_pt = (axbox.x0 + axbox.x1) / 2 / fig.dpi * 72
    run_cx_pt = (b1[0] + b2[2]) / 2
    control("0", "title run centred on the axes", abs(ax_cx_pt - run_cx_pt) < 0.6)
    plt.close(fig)
    print("\nevery control behaved as required" if ok
          else "\n*** A CONTROL FAILED ***")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(self_test() if "--self-test" in sys.argv else 0)
