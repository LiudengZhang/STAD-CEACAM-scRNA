#!/usr/bin/env python3
"""
Append the adaptive-immune composition panel to the submitted S3 page.

S3 is carried over from the submission unaltered. It gains one panel, and it
gains it by composition rather than by re-assembly: the submitted page is drawn
onto a taller canvas exactly as it stands, and a strip carrying the new panel is
drawn below it. Nothing above the join is re-rendered from panel sources, so
nothing above the join can move - which `gate()` proves, pixel row by pixel row,
and which re-assembling the figure could not.

Two things differ from the equivalent patch on the main figures:

  * The strip is placed with an explicit rectangle at 1:1. It is composed at the
    page's own width and dropped at that width, so the type in it prints at the
    size it was set at. Scaling a strip to fit the page width is what puts the
    same 7 pt of type on the paper at some other size.

  * The panel letter is drawn inside the strip's own SVG, by the supplementary
    assembler, so it reaches the PDF as an embedded subset of the same family
    the panels are set in. Written onto the PDF afterwards it would be one of
    the base-14 names, unembedded and rendered in whatever the reader's viewer
    substitutes.

    python patch_S3_add_adaptive.py
"""

import sys
from pathlib import Path

import fitz
import numpy as np

HERE = Path(__file__).resolve().parent               # Supplementary_Fixes
REVISED = HERE.parent                                # 03_Final_Panels
ROOT = REVISED.parent
SUPP = REVISED / "Supplementary_New"

sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(SUPP))

import panel_style_cns as style                      # noqa: E402
import assemble_new_supplementaries as asmb          # noqa: E402

# The submitted page. 00_GROUND_TRUTH/figures is the authority for what was
# submitted; the copy in the reviewer materials is byte-identical to it.
SOURCE = ROOT / "00_GROUND_TRUTH" / "figures" / "S3_CD8_TCells.pdf"
PANEL = HERE / "S3_E" / "S3_E_adaptive_composition.svg"
STRIP_DIR = HERE / "_S3_E_strip"
OUT = HERE / "_patched"
NAME = "S3_CD8_TCells"

MM = 72.0 / 25.4                # mm -> points
TOL = 0.15                      # points, a quarter of a 150 dpi pixel
GATE_DPI = 150


def build_strip():
    """Compose the strip: one panel at its own size, with its letter above it.

    Returns (pdf path, width mm, height mm).
    """
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    _, weight = style.letter_font()
    letter = asmb.printed_letters()["S3_E"]
    w, h = asmb.natural_size_mm(PANEL)
    strip_h = asmb.LETTER_BAND + h
    asm = asmb.PanelAssembler(asmb.PAGE_W, strip_h, title=None)
    asm.place_panel(None, PANEL, asmb.MARGIN, asmb.LETTER_BAND, w, h)
    asm.add_text(letter, asmb.MARGIN + asmb.LETTER_DX, asmb.LETTER_DY,
                 font_size=asmb.LETTER_PT, font_weight=weight)
    STRIP_DIR.mkdir(parents=True, exist_ok=True)
    asm.save(f"{NAME}_strip", output_dir=STRIP_DIR)
    print(f"  strip   {asmb.PAGE_W:.1f} x {strip_h:.1f} mm, panel {letter} "
          f"{w:.1f} x {h:.1f} mm at ({asmb.MARGIN:.1f}, "
          f"{asmb.LETTER_BAND:.1f})  scale 1.000")
    return STRIP_DIR / f"{NAME}_strip.pdf", asmb.PAGE_W, strip_h


def compose(strip_pdf, strip_w_mm, strip_h_mm):
    """The submitted page, unaltered, with the strip below it."""
    src = fitz.open(SOURCE)
    strip = fitz.open(strip_pdf)
    s, r = src[0].rect, strip[0].rect

    if abs(r.width - strip_w_mm * MM) > TOL or \
            abs(r.height - strip_h_mm * MM) > TOL:
        raise AssertionError(
            f"the strip PDF is {r.width / MM:.2f} x {r.height / MM:.2f} mm "
            f"where it was composed at {strip_w_mm:.2f} x {strip_h_mm:.2f} mm")
    if abs(r.width - s.width) > TOL:
        raise AssertionError(
            f"the strip is {r.width / MM:.2f} mm wide and the submitted page "
            f"is {s.width / MM:.2f} mm. The strip is placed at its own size, "
            f"never scaled to the page, so the two widths have to agree - "
            f"compose the strip at the page's width instead.")

    out = fitz.open()
    page = out.new_page(width=s.width, height=s.y1 + r.height)
    # The submitted page goes down at its own rectangle, not at one built from
    # its width and height. Its crop box starts 0.0004 pt below the media box,
    # and placing it at (0, 0) instead moves the whole page by that much: too
    # little to see and enough to change the anti-aliasing of 2060 pixels,
    # which is exactly what the gate below is there to catch.
    page.show_pdf_page(s, src, 0)
    # The strip's rectangle carries the strip's own width and height, so it is
    # placed at the size it was composed at.
    page.show_pdf_page(
        fitz.Rect(0, s.y1, r.width, s.y1 + r.height), strip, 0)

    OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT / f"{NAME}.pdf"
    if dst.exists():
        dst.chmod(0o644)
    out.save(dst, garbage=3, deflate=True)
    out.close()
    dst.chmod(0o644)
    print(f"  page    {s.width / MM:.2f} x {(s.y1 + r.height) / MM:.2f} mm"
          f"   join at y = {s.y1 / MM:.2f} mm")
    join, width = s.y1, s.width
    src.close()
    strip.close()
    return dst, join, width


def gate(dst, join_pt, width_pt):
    """Every whole pixel above the join must be the pixel the submitted page has.

    Rendered rather than compared as objects: the question is what the reader
    sees, and a PDF can be rewritten byte for byte without changing it.

    "Whole pixel" is the raster the two pages share, and it is a definition
    rather than a tolerance. The page is 1080.71 pixels wide at this
    resolution and the join falls 0.71 of a pixel into a row, so the last
    column and the row across the join are each only partly covered by the
    page. A partly covered pixel is not a pixel above the join; it is the edge
    of the paper, and the two renders resolve it one level apart out of 255
    because one arrives through a form and the other does not. The count over
    those edges is reported too, so nothing is hidden.
    """
    a = fitz.open(SOURCE)
    b = fitz.open(dst)
    pa = a[0].get_pixmap(dpi=GATE_DPI)
    pb = b[0].get_pixmap(dpi=GATE_DPI)
    if pa.width != pb.width or pa.n != pb.n:
        raise AssertionError(
            f"the pages render at different widths ({pa.width} vs {pb.width} "
            f"px) and cannot be compared row by row")
    rows = min(int(join_pt / 72.0 * GATE_DPI), pa.height, pb.height)
    cols = min(int(width_pt / 72.0 * GATE_DPI), pa.width)
    A = np.frombuffer(pa.samples, np.uint8).reshape(pa.height, pa.width, pa.n)
    B = np.frombuffer(pb.samples, np.uint8).reshape(pb.height, pb.width, pb.n)
    moved = (A[:rows] != B[:rows]).any(axis=2)
    differing = int(moved[:, :cols].sum())
    edge = int(moved.sum()) - differing
    total = rows * cols
    worst = int(np.abs(A[:rows, :cols].astype(int)
                       - B[:rows, :cols].astype(int)).max())
    a.close()
    b.close()
    print(f"  gate    {differing} of {total} whole pixels differ above the "
          f"join ({rows} rows x {cols} columns at {GATE_DPI} dpi, "
          f"largest difference {worst} of 255)")
    print(f"          {edge} pixel(s) differ in the partly covered edge "
          f"column, which is not above the join")
    if differing:
        raise AssertionError(
            f"{differing} pixel(s) above the join differ from the submitted "
            f"page. The submitted page is carried over unaltered; anything "
            f"above the join that moves is a fault in the composition, not a "
            f"tolerance.")
    return differing


def main():
    print(f"{NAME}: appending the adaptive-immune composition panel")
    strip_pdf, w, h = build_strip()
    dst, join, width = compose(strip_pdf, w, h)
    gate(dst, join, width)
    print(f"\nWritten to {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
