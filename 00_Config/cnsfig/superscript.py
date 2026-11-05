"""
Superscript "+" for cell-type markers (CD4+, CD8+, CEACAM5+, CEACAM6+, DAB+)
and superscript exponents for Unicode powers of ten ("10⁻⁷").

superscripts(fig) is a context manager used by panel_style_cns.save_panel:
while it is active, matching strings are set through mathtext (CD8$^{+}$, at
0.7 of the label size, in the panel's own font), and on exit the original
strings are restored. The glyphs stay real text in the PDF. The arithmetic
"+" of "log2(x+1)" is not matched.

    python superscript.py --self-test
"""

import re
import sys
from contextlib import contextmanager

# A marker name standing as a word, then its "+".
MARKER_RE = re.compile(r"(?<![A-Za-z0-9_])(CD4|CD8|CEACAM5|CEACAM6|DAB)\+")
# mathtext's superscript scale.
SUP_SCALE = 0.7


# A power of ten in Unicode superscripts ("10⁻⁷", from corr_stats.p_string).
# The panel font has no superscript minus, so the exponent goes through mathtext.
POWER_RE = re.compile(r"10([⁻]?)([⁰¹²³⁴⁵⁶⁷⁸⁹]+)")
_DIGITS = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")


def superscripted_text(s):
    """`s` with each marker "+" and each Unicode power of ten as a mathtext
    superscript; None if unchanged or already mathtext (a "$" in it)."""
    if not s or "$" in s or not (MARKER_RE.search(s) or POWER_RE.search(s)):
        return None
    s = MARKER_RE.sub(lambda m: m.group(1) + "$^{+}$", s)
    return POWER_RE.sub(lambda m: "10$^{" + ("-" if m.group(1) else "")
                        + m.group(2).translate(_DIGITS) + "}$", s)


def _segments(line):
    """[(text, is_superscript)] for one plain line: a marker's "+" and a power
    of ten's exponent (its minus as U+2212, which the face has) are the
    superscript pieces; everything else is plain."""
    out, pos = [], 0
    pat = re.compile(MARKER_RE.pattern + "|" + POWER_RE.pattern)
    for m in pat.finditer(line):
        if m.group(1):                      # a marker: name plain, "+" raised
            out.append((line[pos:m.start()] + m.group(1), False))
            out.append(("+", True))
        else:                               # 10 plain, the exponent raised
            out.append((line[pos:m.start()] + "10", False))
            out.append((("\u2212" if m.group(2) else "")
                        + m.group(3).translate(_DIGITS), True))
        pos = m.end()
    out.append((line[pos:], False))
    return [(t, sup) for t, sup in out if t]


@contextmanager
def superscripts(fig, family=None):
    """Within the block, every marker "+" on `fig` is a superscript; on exit
    every string is restored. Yields the number of strings changed."""
    import matplotlib
    from matplotlib.text import Text
    from matplotlib.ticker import FuncFormatter
    changed = []
    for t in fig.findobj(Text):
        new = superscripted_text(t.get_text())
        if new is not None:
            changed.append((t, t.get_text()))
            t.set_text(new)
    # Tick labels are regenerated at draw time, so wrap their formatter instead.
    formatters = []
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            orig = axis.get_major_formatter()
            labels = [orig(v, i) for i, v in enumerate(axis.get_majorticklocs())]
            if not any(superscripted_text(str(lab)) for lab in labels):
                continue

            def wrapped(v, pos, _orig=orig):
                out = _orig(v, pos)
                return superscripted_text(out) or out
            formatters.append((axis, orig))
            axis.set_major_formatter(FuncFormatter(wrapped))
    rc = {}
    if changed or formatters:
        if family is None:
            import panel_style_cns as style
            family = style.resolved_family()
        rc = {"mathtext.fontset": "custom", "mathtext.rm": family}
    split = []          # (hidden original, [per-line stand-ins])
    try:
        with matplotlib.rc_context(rc):
            # mathtext ignores linespacing in multi-line strings, so these are
            # drawn as one Text per line at the plain layout's baselines, with
            # the superscript pieces raised; the original is hidden meanwhile.
            multi = [(t, s) for t, s in changed if "\n" in s]
            if multi:
                import math
                import matplotlib.transforms as mtransforms
                r = fig.canvas.get_renderer()
                k = 72.0 / fig.dpi
                for t, s in multi:
                    t.set_text(s)
                    _, info, _ = t._get_layout(r)
                    prop = t.get_fontproperties()
                    size = prop.get_size_in_points()
                    a = math.radians(t.get_rotation())
                    along = (math.cos(a), math.sin(a))      # reading direction
                    up = (-math.sin(a), math.cos(a))
                    parts = []

                    def put(text, x_pt, y_pt, fp):
                        tr = mtransforms.offset_copy(
                            t.get_transform(), fig=fig, x=x_pt, y=y_pt,
                            units="points")
                        parts.append(fig.text(
                            *t.get_unitless_position(), text, transform=tr,
                            ha="left", va="baseline", rotation=t.get_rotation(),
                            rotation_mode="anchor", fontproperties=fp,
                            color=t.get_color(), zorder=t.get_zorder()))
                    sup_fp = prop.copy()
                    sup_fp.set_size(size * SUP_SCALE)
                    for line, _wh, x, y in info:
                        if not line.strip():
                            continue
                        x0, y0 = x * k, y * k      # this line's baseline, pt
                        cursor = 0.0
                        rise = 0.36 * size
                        for piece, is_sup in _segments(line):
                            fp = sup_fp if is_sup else prop
                            lead = len(piece) - len(piece.lstrip(" "))
                            if lead:
                                cursor += r.get_text_width_height_descent(
                                    " " * lead, fp, False)[0] * k
                            body = piece.strip(" ")
                            if body:
                                dr = rise if is_sup else 0.0
                                put(body, x0 + cursor * along[0] + dr * up[0],
                                    y0 + cursor * along[1] + dr * up[1], fp)
                                cursor += r.get_text_width_height_descent(
                                    body, fp, False)[0] * k
                            trail = len(piece) - len(piece.rstrip(" "))
                            if trail and body:
                                cursor += r.get_text_width_height_descent(
                                    " " * trail, fp, False)[0] * k
                    t.set_visible(False)
                    split.append((t, parts))
            yield len(changed) + len(formatters)
    finally:
        for t, parts in split:
            for part in parts:
                part.remove()
            t.set_visible(True)
        for t, s in changed:
            t.set_text(s)
        for axis, orig in formatters:
            axis.set_major_formatter(orig)


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        from superscript_selftest import self_test
        sys.exit(self_test())
