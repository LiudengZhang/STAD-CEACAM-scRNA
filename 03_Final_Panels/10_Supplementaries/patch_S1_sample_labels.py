"""
Put the study sample IDs on the sample axes of Supplementary Figure S1.

The submitted S1 prints internal hospital specimen identifiers - the codes the
samples carry in the lab, in six different naming formats - on the x axis of
the per-sample QC box plots, where the paper names every sample P01-M1 style.
Panel H already uses the study IDs; the box-plot panels never did.

The figure is patched, not rebuilt. Re-running assemble_S1.py does not
reproduce the submitted page - it comes out 1187 pt tall against the submitted
877 pt and loses seven of the eight panel letters - so a rebuild would replace
the published figure rather than de-identify it (CLAUDE.md rules 1 and 3).
Redacting the tick glyphs and redrawing the study ID on the same anchor leaves
every box, whisker, axis, legend and statistic exactly as submitted; the tick
text is the only thing that changes.

Which panels those are is read out of the page, never inferred from a name
(CLAUDE.md rule 2): the panel letters are located by position and a panel is
patched because it contains rotated tick labels that are specimen identifiers,
not because of what it is called. The panels found are printed at the end of
the run.

The crosswalk is read from the dataset at run time and cross-checked against
Supplementary Table 1, so no internal identifier is written into this file or
into anything it produces. A specimen with no study ID is an error, not a
label to fall back on.

Run: python patch_S1_sample_labels.py
"""

from pathlib import Path
import re
import sys

import fitz
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from paths import FULL_DATASET_H5AD, MANUSCRIPT, REVIEWER_MATERIALS  # noqa: E402

FIGURE = "S1_QC_Annotation.pdf"
SRC = REVIEWER_MATERIALS / "figures_submitted" / FIGURE
OUT_DIR = Path(__file__).parent / "_patched"
ST1_CSV = MANUSCRIPT / "04_Tables" / "ST1_patient_sample_characteristics.csv"

STUDY_ID = re.compile(r"P\d+-\w+")     # the Supplementary Table 1 label format
LETTER = re.compile(r"[A-Z]")
LETTER_SIZE = 8.0                      # pt; panel letters are 10, nothing else is
TICK_SIZE = 4.0                        # pt; tick labels are ~2.7, axis titles ~5.3
UPWARD = (0.0, -1.0)                   # text direction of a rotation=90 tick label
FONT = "helv"                          # metrically the ArialMT the figure uses


def crosswalk():
    """Internal specimen identifier -> study sample ID, for every sample."""
    obs = sc.read_h5ad(FULL_DATASET_H5AD, backed="r").obs
    pairs = obs[["sample", "Sample ID"]].drop_duplicates().astype(str)
    st1 = set(pd.read_csv(ST1_CSV)["Sample"].astype(str))
    # Eight lymph-node specimens carry two different study IDs across their own
    # cells, differing in whether the node is numbered or marked positive, so
    # the column alone does not define the label. Supplementary Table 1 is the
    # authority and lists exactly one of each pair; a candidate it does not list
    # is not a study ID whatever it looks like.
    keep = pairs[pairs["Sample ID"].isin(st1)]
    lost = sorted(set(pairs["sample"]) - set(keep["sample"]))
    if lost:
        raise SystemExit(f"{len(lost)} specimen(s) with no {ST1_CSV.name} label")
    clash = keep["sample"].duplicated().sum()
    if clash:
        raise SystemExit(f"{clash} specimen(s) with two {ST1_CSV.name} labels")
    return dict(keep.values)


def panel_bands(page):
    """Printed panel letter -> (y_top, y_bottom), read off the page itself."""
    letters = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                if span["size"] >= LETTER_SIZE and LETTER.fullmatch(span["text"].strip()):
                    letters.append((round(span["bbox"][1], 1), span["text"].strip()))
    tops = sorted({y for y, _ in letters})
    bands = {}
    for y, name in letters:
        below = [t for t in tops if t > y]
        bands.setdefault(name, (y, below[0] if below else page.rect.y1))
    return bands


def tick_spans(page):
    """The rotated, tick-sized text spans, in reading order along the axis."""
    out = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            if line["dir"] != UPWARD:
                continue
            for span in line["spans"]:
                if span["size"] < TICK_SIZE:
                    out.append(span)
    return sorted(out, key=lambda s: (s["bbox"][1], s["bbox"][0]))


def main():
    ids = crosswalk()
    doc = fitz.open(SRC)
    page = doc[0]
    bands = panel_bands(page)
    ticks = tick_spans(page)

    # A tick label that is neither a known specimen nor already a study ID is
    # an unmapped specimen: stop rather than ship a partial de-identification.
    unmapped = [s for s in ticks
                if s["text"] not in ids and not STUDY_ID.fullmatch(s["text"])]
    if unmapped:
        raise SystemExit(f"{len(unmapped)} axis label(s) with no study ID; "
                         "the crosswalk is incomplete")

    todo, panels = [], {}
    for span in ticks:
        if span["text"] not in ids:
            continue
        letter = next(n for n, (top, bot) in bands.items()
                      if top <= span["bbox"][1] < bot)
        panels.setdefault(letter, set()).add(span["text"])
        todo.append((span, ids[span["text"]], letter))

    for span, _, _ in todo:
        page.add_redact_annot(fitz.Rect(span["bbox"]))
    # Text only. The boxes, whiskers and axes are line art and must survive.
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                          graphics=fitz.PDF_REDACT_LINE_ART_NONE)

    for span, new, _ in todo:
        # Rotation 90 reads upward, so the run grows downward from the top edge
        # as it gets longer: the anchor is the top of the span it replaces, and
        # the baseline x is the one the submitted label already sits on.
        length = fitz.get_text_length(new, fontname=FONT, fontsize=span["size"])
        page.insert_text((span["origin"][0], span["bbox"][1] + length), new,
                         fontname=FONT, fontsize=span["size"], rotate=90,
                         color=(0, 0, 0))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dst = OUT_DIR / FIGURE
    if dst.exists():
        dst.chmod(0o644)
    doc.save(dst, garbage=3, deflate=True)
    doc.close()
    dst.chmod(0o644)

    for letter in sorted(panels):
        n = sum(1 for _, _, p in todo if p == letter)
        print(f"panel {letter}: {len(panels[letter])} samples, {n} tick labels")
    print(f"{len(todo)} labels replaced; written to {dst}")


if __name__ == "__main__":
    main()
