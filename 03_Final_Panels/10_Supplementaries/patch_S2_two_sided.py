"""
Print the two-sided P value on panel D of Supplementary Figure S2.

The submitted S2D - the CEACAM5/6+ epithelial proportion adjusted for tumor
content, pre-treatment non-responders against responders, four samples each -
prints "P = 0.03". That is the one-sided exact Mann-Whitney value for the
configuration the panel plots (U = 1 at 4 vs 4: 2/70 = 0.029). The revision
reports two-sided tests throughout and Table S6 tabulates every converted
comparison, including this one at 4/70 = 0.057. The panel therefore has to say
what the table says.

The figure is patched, not rebuilt: the live panel script draws a different
title from the one submitted, so a rebuild would replace the published panel
rather than correct its one annotation (the project's standing rules rules 1 and 3). The text
span is redacted and the new value drawn on the same origin at the same size;
every dot, box, whisker, axis and bracket is line art and survives untouched.

Which span that is is read out of the page, never assumed: the panel letters
are located by their size, panel D's box is the region from its letter to the
next letter in the row and the next row of letters, and the span replaced is
the one span inside that box whose text is exactly the old value. Zero or
several matches is an error, not a choice. After saving, the page is rendered
against the submitted one at 150 dpi and the run fails if any pixel outside
the old and new text boxes changed.

Run: python patch_S2_two_sided.py
"""

from pathlib import Path
import sys

import fitz
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from paths import REVIEWER_MATERIALS  # noqa: E402

FIGURE = "S2_CEACAM_Metaprogram_Validation.pdf"
SRC = REVIEWER_MATERIALS / "figures_submitted" / FIGURE
OUT_DIR = Path(__file__).parent / "_patched"

PANEL = "D"
OLD = "P = 0.03"                       # one-sided exact value, as submitted
NEW = "P = 0.057"                      # two-sided exact value, Table S6
LETTER_SIZE = 11.0                     # pt; panel letters are 11.5, nothing else is above 10.6
FONT = "helv"                          # metrically the ArialMT the panel uses
DPI = 150
MARGIN_PT = 1.0                        # antialiasing halo allowed around the text boxes


def panel_box(page, letter):
    """The rectangle of a printed panel, from its letter to the next ones."""
    letters = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                text = span["text"].strip()
                if span["size"] >= LETTER_SIZE and len(text) == 1 and text.isupper():
                    letters.append((text, fitz.Rect(span["bbox"])))
    mine = [r for name, r in letters if name == letter]
    if len(mine) != 1:
        raise SystemExit(f"panel letter {letter!r} found {len(mine)} times, expected 1")
    me = mine[0]
    same_row = [r for _, r in letters if abs(r.y0 - me.y0) < 5 and r.x0 > me.x0 + 5]
    below = [r for _, r in letters if r.y0 > me.y0 + 5]
    x1 = min(r.x0 for r in same_row) if same_row else page.rect.x1
    y1 = min(r.y0 for r in below) if below else page.rect.y1
    return fitz.Rect(me.x0, me.y0, x1, y1)


def spans_reading(page, text, box):
    """Every text span reading exactly `text` whose box lies inside `box`."""
    out = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                if span["text"].strip() == text and box.contains(fitz.Rect(span["bbox"])):
                    out.append(span)
    return out


def render(path):
    with fitz.open(path) as doc:
        pix = doc[0].get_pixmap(dpi=DPI, colorspace=fitz.csGRAY)
    return np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width)


def main():
    doc = fitz.open(SRC)
    if doc.page_count != 1:
        raise SystemExit(f"{SRC.name} has {doc.page_count} pages, expected 1")
    page = doc[0]
    box = panel_box(page, PANEL)

    hits = spans_reading(page, OLD, box)
    if len(hits) != 1:
        raise SystemExit(f"{OLD!r} found {len(hits)} times inside panel {PANEL} "
                         f"{tuple(round(v, 2) for v in box)}, expected exactly 1")
    span = hits[0]
    on_page = len(page.search_for(OLD))
    old_rect = fitz.Rect(span["bbox"])
    origin = fitz.Point(span["origin"])
    size = span["size"]
    print(f"panel {PANEL} box      : {tuple(round(v, 2) for v in box)} pt")
    print(f"span '{OLD}'      : bbox {tuple(round(v, 2) for v in old_rect)} pt, "
          f"origin {tuple(round(v, 2) for v in origin)}, font {span['font']} "
          f"{size:.2f} pt; {on_page} occurrence(s) on the whole page")
    old_w = fitz.get_text_length(OLD, fontname=FONT, fontsize=size)
    new_w = fitz.get_text_length(NEW, fontname=FONT, fontsize=size)
    print(f"{FONT} width of '{OLD}' at {size:.2f} pt = {old_w:.2f} pt against "
          f"{old_rect.width:.2f} pt printed; '{NEW}' will be {new_w:.2f} pt wide")
    new_rect = fitz.Rect(origin.x, old_rect.y0, origin.x + new_w, old_rect.y1)
    if not box.contains(new_rect):
        raise SystemExit("the two-sided value would not fit inside the panel")

    page.add_redact_annot(old_rect)
    # Text only. The dots, boxes, whiskers, axes and bracket are line art.
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                          graphics=fitz.PDF_REDACT_LINE_ART_NONE)
    page.insert_text(origin, NEW, fontname=FONT, fontsize=size, color=(0, 0, 0))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dst = OUT_DIR / FIGURE
    if dst.exists():
        dst.chmod(0o644)
    doc.save(dst, garbage=3, deflate=True)
    doc.close()
    dst.chmod(0o644)

    # ---- the patched file says what it should, and nothing else moved
    with fitz.open(dst) as check:
        text = check[0].get_text()
        if text.count(NEW) != 1 or OLD in text:
            raise SystemExit(f"{dst.name}: expected {NEW!r} once and {OLD!r} never; "
                             f"got {text.count(NEW)} and {text.count(OLD)}")
        if len(spans_reading(check[0], NEW, box)) != 1:
            raise SystemExit(f"{NEW!r} is not inside panel {PANEL} in {dst.name}")
    before, after = render(SRC), render(dst)
    if before.shape != after.shape:
        raise SystemExit(f"page size changed: {before.shape} -> {after.shape}")
    changed = before != after
    allowed = (old_rect | new_rect) + (-MARGIN_PT, -MARGIN_PT, MARGIN_PT, MARGIN_PT)
    scale = DPI / 72.0
    ys, xs = np.nonzero(changed)
    if len(ys):
        px_box = (xs.min() / scale, ys.min() / scale, (xs.max() + 1) / scale,
                  (ys.max() + 1) / scale)
        outside = ~((xs / scale >= allowed.x0) & (xs / scale <= allowed.x1)
                    & (ys / scale >= allowed.y0) & (ys / scale <= allowed.y1))
        if outside.any():
            raise SystemExit(f"{int(outside.sum())} changed pixel(s) fall outside "
                             f"the text boxes {tuple(round(v, 2) for v in allowed)} "
                             f"pt; changed region spans {px_box}")
    else:
        px_box = None
    frac = changed.mean()
    print(f"{DPI} dpi diff vs submitted: {int(changed.sum())} of {changed.size} "
          f"pixels changed ({100 * frac:.4f}%), all inside "
          f"{tuple(round(v, 2) for v in allowed)} pt"
          + (f"; changed region {tuple(round(v, 2) for v in px_box)} pt" if px_box else ""))
    print(f"'{OLD}' -> '{NEW}' in panel {PANEL}; written to {dst}")


if __name__ == "__main__":
    main()
