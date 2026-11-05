"""
The frame every panel script draws inside: slot, canvas, guards, save.

Each panel script builds one `Panel(figure, printed_letter, __file__)` at
module level - MARK is captured at construction, before `apply()` changes the
tick size - and draws on its canvas at the millimetre size it prints at.

USAGE

    from cnsfig import Panel, pin_frame_mm, style

    PANEL = Panel(3, "F", __file__, small_pt=6.0, letter_gutter_mm=1.0)
    MARK, AREA = PANEL.mark, PANEL.area

    def main():
        fig, ax = PANEL.subplots()
        ...                                   # the panel's own drawing
        PANEL.fit(fig)
        pin_frame_mm(fig, ax, top_mm=5.0, bottom_mm=PANEL.h_mm - 29.0)
        PANEL.finish(fig, 'spatial_ceacam_ratio')
"""


import sys
from pathlib import Path

# Put 00_Config on the path when run directly (the self-test).
_CONFIG = Path(__file__).resolve().parent.parent
if str(_CONFIG) not in sys.path:
    sys.path.insert(0, str(_CONFIG))

import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402

__all__ = ["Panel"]

# Canvas multiplier used to convert 4x marker sizes to a 1:1 canvas (MARK).
SCALE = 4

# Type sizes passed to style.apply() by every panel.
TYPE_KW = dict(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)


class Panel:
    """One printed panel's slot, canvas, guards and output path.

    Parameters
    ----------
    figure
        The printed figure number.
    letter
        The printed panel letter (not necessarily the directory's letter).
    script
        `__file__` of the panel script. Used for the directory the drawing is
        written to, and for nothing else.
    small_pt
        Smallest body type at the 4x scale; sets MARK and AREA.
    letter_gutter_mm
        Millimetres added to the panel-letter keep-out cell's width when
        `fit()` runs. `finish()` checks the cell without the gutter.

    Build the Panel at module level: MARK is computed from the tick size,
    which apply() changes.
    """

    SCALE = SCALE

    def __init__(self, figure, letter, script, *, small_pt,
                 letter_gutter_mm=0.0):
        self.figure = figure
        self.letter = str(letter)
        self.dir = Path(script).resolve().parent
        self.small_pt = float(small_pt)
        self.letter_gutter_mm = float(letter_gutter_mm)

        # Printed slot size, from slots.py.
        self.w_mm, self.h_mm = slots.size_mm(figure, self.letter)
        self.letter_cell = slots.letter_cell_mm(figure, self.letter)
        self.letter_fit = (self.letter_cell[0] + self.letter_gutter_mm,
                           self.letter_cell[1])

        # Tick size before apply() (7.0); after apply() it is 6.0.
        self.tick_pt_at_build = style.tick_pt()
        self.mark = self.tick_pt_at_build / (self.small_pt * SCALE)
        self.area = self.mark ** 2

        self._applied = False

    # -------------------------------------------------------------- the style
    def apply(self, **overrides) -> str:
        """Apply the cnsplots style at TYPE_KW sizes. Idempotent.

        Returns the font family actually used.
        """
        family = style.apply(**{**TYPE_KW, **overrides})
        if not self._applied:
            print(f"  type set in {family}; body {style.body_pt():g} pt, "
                  f"ticks/legend {style.tick_pt():g} pt; MARK {self.mark:.3f}")
            self._applied = True
        return family

    # ------------------------------------------------------------- the canvas
    def subplots(self, *args, **kwargs):
        """`style.subplots_mm` at this panel's slot, with the style applied."""
        self.apply()
        return style.subplots_mm(self.w_mm, self.h_mm, *args, **kwargs)

    def canvas(self, **kwargs):
        """`style.figure_mm` at this panel's slot, for a figure whose axes are
        placed by `add_axes` (cnsfig.layout) rather than by subplots."""
        self.apply()
        return style.figure_mm(self.w_mm, self.h_mm, **kwargs)

    # -------------------------------------------------------------- the margin
    def fit(self, fig, pad_mm=0.6, *, reserve_letter=True):
        """`style.fit_margins` with this panel's letter cell and gutter."""
        return style.fit_margins(fig, pad_mm=pad_mm, cell_mm=self.letter_fit,
                                 reserve_letter=reserve_letter)

    # -------------------------------------------------------- guards and save
    def check(self, fig):
        """Raise if any ink left the canvas or entered the letter cell.

        Measured on the rendered figure; called by `finish()` before saving.
        """
        over = style.overflow_mm(fig)
        if max(over) > 0:
            raise RuntimeError(
                f"ink outside the {self.w_mm} x {self.h_mm} mm canvas "
                f"(l,r,b,t mm): {over}")
        intruders = style.letter_clear(fig, self.letter_cell)
        if intruders:
            raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    def finish(self, fig, stem, **save_kw):
        """Check the two guards, then write `<panel dir>/<stem>` 1:1.

        `stem` is a bare name; the directory is the script's own. bbox_inches
        is refused, because cropping would change the printed size.
        """
        if "bbox_inches" in save_kw:
            raise TypeError(
                "bbox_inches would crop the canvas to the ink, and a panel is "
                "drawn at the millimetre size it prints at. Refused here so "
                "that it cannot be passed panel by panel.")
        self.check(fig)
        out = self.dir / str(stem)
        written = style.save_panel(fig, out, **save_kw)
        print(f"  Saved: {out}.[svg|pdf|png] at {self.w_mm:g} x "
              f"{self.h_mm:g} mm")
        return written

    def __repr__(self):
        return (f"Panel(figure={self.figure}, letter={self.letter!r}, "
                f"{self.w_mm:g} x {self.h_mm:g} mm, mark={self.mark:.4f})")


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        from panel_selftest import self_test
        sys.exit(self_test())
    print(__doc__)
    sys.exit(0)
