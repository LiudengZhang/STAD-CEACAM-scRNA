#!/usr/bin/env python3
"""
Assemble the restyled main Figures 1-5 (Version B).

This is a rewrite of `01_Figure_1/assemble_figure_1.py` ...
`05_Figure_5/assemble_figure_5.py`, not an edit of them. Version A is frozen
and those five files stay exactly as they are.

WHY A REWRITE AND NOT A REUSE
-----------------------------
The five assemblers on disk do not build the figure the paper prints. Measured
against `PROVENANCE.csv`:

    assemble_figure_2.py places  A B C D G H I J K M N1 N2 O1 O2 P Q   (16)
    Figure 2 prints              A B C D E F G H I J K L M N           (14)

    assemble_figure_5.py places  A B C D E F G H I J K L M N O P       (16)
    Figure 5 prints              A B C D E F G H I J K L M N           (14)

    assemble_figure_4.py places  A B C D E F G                          (7)
    Figure 4 prints              A B C D E F G H                        (8)

Figure 5's assembler splits the printed panel J - which is four CD274 boxplots
under one letter - into four letters J K L M, and every later panel then shifts
by three: printed K becomes N, printed L becomes O, printed M becomes P, and
printed N (the four-cell-type GSEA summary) is not placed at all. Figure 2's
splits N and O into N1/N2 and O1/O2, drops printed E, F and L, and adds a
Jaccard heatmap the paper does not print. Figure 4's drops printed H.

Running any of them would therefore *renumber the paper*. This assembler does
not: **it is driven by the printed letters in PROVENANCE.csv**, so Version B
carries the same letters on the same content as the published figures, and no
reference in the main text, the response letter or Table S6 has to move for a
panel that Version B can draw. Where Version B cannot draw a panel that is a
loss, not a renumbering, and it is recorded in `../../RENUMBERING.md`.

CLAUDE.md rule 2 in force: **a panel directory's letter is not the printed
panel letter.** `02_H` holds printed panel G, `02_G` holds H *and* I, `05_I`
holds J, `05_F` holds H, `05_C` holds G. Nothing here infers a letter from a
directory name; `LAYOUT` names printed letters and `resolve()` looks up which
drawing supplies each one. `check_layout_against_provenance()` refuses to
assemble if the two disagree.

WHAT IS DIFFERENT FROM VERSION A, MECHANICALLY
----------------------------------------------
Panels are placed at their natural size and `check_scale()` refuses to assemble
if any placement is not 1:1. Version A declared a millimetre box per panel and
let `preserveAspectRatio='xMidYMid meet'` shrink the panel into it, which is
what put identical source type on the page anywhere between 1.7 and 7.5 pt.
Type set at 7 pt here prints at 7 pt.

    conda run -n Liudeng_Python_310 python assemble_main_figures.py [Figure_2 ...]
"""

import csv
import re
import sys
from collections import OrderedDict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAIN = HERE.parent                       # 11_Version_B_Restyle/Main_Figures
ROOT = MAIN.parents[1]                   # repository root
sys.path.insert(0, str(ROOT / "00_Config"))

import panel_style_cns as style                       # noqa: E402
from shared.svg_assembler import VectorAssembler      # noqa: E402
from lxml import etree                                # noqa: E402

OUT = MAIN / "_restyled"
PROVENANCE = ROOT / "03_Revised_Panels" / "PROVENANCE.csv"
RECTS = ROOT / "03_Revised_Panels" / "panel_rects.csv"

PAGE_W = style.PAGE_W_MM        # 190.5 mm, cnsplots' full-width figure
GUTTER = 4.0                    # between panels in a row
ROW_GAP = 5.0                   # between rows
TOP = 4.0
BOTTOM = 4.0
LETTER_BAND = 4.6               # mm reserved above every row for its letters
LETTER_DX = 0.0
LETTER_DY = 3.4

MM = 25.4 / 72.0                # points -> mm
_SVGNS = "http://www.w3.org/2000/svg"

# What a main figure can occupy once the journal has placed it: A4 less 20 mm
# margins. Reported for every figure, never enforced - the same treatment
# stage 2 gave the supplementary pages, and for the same reason: the decision
# to accept a taller page than the text block belongs to the author and to
# production, not to this script.
PRINTABLE_W, PRINTABLE_H = 170.0, 247.0


# ---------------------------------------------------------------------------
# Layout. Rows of printed panel letters, top to bottom.
# ---------------------------------------------------------------------------
# A row is a list of cells; a cell is a printed letter, or a tuple of printed
# letters that share ONE drawing (Figure 2's H and I are the two halves of a
# single side-by-side boxplot SVG), or a list of letters whose drawings sit
# side by side under one letter (Figure 2's K is 02_J and 02_K).
#
# The letters here are the letters the PAPER prints. They are checked against
# PROVENANCE.csv before anything is placed.

# Letters that share ONE drawing. Figure 2's H and I are the left and right
# halves of a single side-by-side boxplot SVG, so they are placed once and
# lettered twice. Everything else is derived.
SHARED = {("Figure_2", "H"): ("H", "I")}

# The rows are NOT hand-written. They are packed by `pack_rows()` in **printed
# letter order**, greedily, at the sizes the panels were actually drawn at.
#
# Letter order is not a stylistic preference: a figure's panel letters have to
# read left to right and top to bottom, so A may not sit beside C with B on the
# next row. That constraint is what makes the pages tall - two panels that are
# each just over half the page width cannot share a row, and the second one
# takes a row of its own. The per-figure numbers are printed on every run.
#
# A hand-written table was tried first and thrown away: it put panels out of
# reading order to save height, which is a worse fault than a tall page.
LAYOUT = None

# Printed panels Version B CANNOT carry, and why. A panel gets in here only
# when its Version A script does not run against the data on disk today, so the
# redraw cannot be proved CONTENT IDENTICAL against it. CLAUDE.md rule 1: the
# published figure is the ground truth and the code is the suspect, so a script
# that fails is reported, never repaired - and a panel that cannot be gated
# must not ship.
#
# This is NOT a renumbering. The letters of the panels that ARE carried do not
# move; an omitted letter simply is not drawn, and every citation of it in the
# main text, the response letter and Table S6 would be left pointing at nothing.
# That is the cost, and it is spelled out per panel in ../../RENUMBERING.md.
#
# Populated from the runnability sweep, `_audit/runnability.csv` +
# `_audit/runnability_divert.log`: 61 of the 63 Version A main-figure scripts
# run. These are the two that do not.
# Figure 2's printed panel A was here until 2026-09-01. It was STOPPED, not
# failed: every one of its 149,373 points, all nine label strings and every
# colour were identical, but 55 label POSITIONS moved, because `adjust_text`
# repels labels by their RENDERED bounding box and Version B's canvas is
# 66 mm where Version A's was 280 mm. It is now carried, with the nine label
# positions and the six leader lines pinned to the values measured off Version
# A - see 02_A/create_epithelial_umap.py and ../../PANEL_2A_LABEL_PINS.md. The
# gate reports CONTENT IDENTICAL and no label needed a nudge at print size.
OMITTED = {
    ("Figure_1", "C"):
        "01_C/generate_stomach_stacked_bar.py raises KeyError: "
        "'Original Sample ID'. It maps the h5ad's sample IDs onto the printed "
        "ones through that column of "
        "Round_5/04_Manuscript/04_Tables/ST1_patient_sample_characteristics.csv, "
        "and the table on disk today has no such column - its columns are "
        "Patient ID, Sample, Age, Sex, cTNM stage, Differentiation, Stomach "
        "site, Anatomical site, Tumor type, Biopsy method, Treatment phase, "
        "R/NR Grouping, Other Notes. Pre-existing and nothing to do with the "
        "restyle: Version A's own script fails identically. CLAUDE.md rule 1 - "
        "the published panel is the ground truth and the script is the "
        "suspect - so it is reported, not repaired.",
    ("Figure_2", "F"):
        "02_F/create_epithelial_milo_pre_rvsnr.py cannot run in "
        "Liudeng_Python_310 (ModuleNotFoundError: jaxlib.xla_extension, as its "
        "own docstring predicts) and cnsplots cannot run in pertpy_milo, where "
        "the panel does run. 22 of cnsplots' 38 declared requirements are "
        "absent there and 8 are needed at import time. Installing them would "
        "move versions in an environment that produces a shipped panel, which "
        "is exactly what CLAUDE.md's --no-deps rule forbids. Measured, "
        "attempted and fully reverted: see _audit/FIGURE_2F_ENVIRONMENT.md.",
}

# Figure 6 has no panel code. `Round_5/03_Final_Panels/06_Figure_6/canvas.svg`
# is an earlier hand-drawn draft of the same concept, not the published
# artwork; nothing on disk can regenerate the shipped PDF. Version B therefore
# has no Figure 6 and Version A's Figure 6 stands. See ../../RENUMBERING.md.
NO_CODE = {"Figure_6": "hand-drawn schematic; no panel script exists anywhere"}


class RestyledAssembler(VectorAssembler):
    """VectorAssembler with the cnsplots type system for its own text.

    Version A's assembler hardcodes `Liberation Sans` and draws panel letters
    at `font_weight='normal'`. cnsplots sets panel letters in the panel-label
    font at `title_fontsize`, bold. Overridden here rather than edited in
    place, because `svg_assembler.py` is shared with Version A and Version A
    is frozen. Identical to the stage-2 supplementary assembler's override.
    """

    FONT_STACK = None

    def __init__(self, *a, **kw):
        family, _ = style.letter_font()
        # ONE family, named exactly, with NO generic fallback.
        #
        # The stage-2 assembler's comment says "one family, named exactly" but
        # its code writes `f"{family}, sans-serif"`, which is a list. Measured
        # on the first Version B main-figure build: Figures 1-4 came out with
        # their panel letters in NimbusSans-Bold and Figure 5's fourteen came
        # out in LiberationSans-Bold, off the identical string in the identical
        # run. cairo resolves a CSS font list by its own rules and may take the
        # generic; matplotlib's font manager does not. That divergence is one
        # of the mechanisms behind the thirteen font families measured across
        # the shipped figures, and a trailing generic is all it needs.
        self.FONT_STACK = family
        super().__init__(*a, **kw)

    def _embed_text(self, root, data):
        x, y = data["x_mm"], data["y_mm"]
        txt = etree.SubElement(root, "{%s}text" % _SVGNS)
        txt.set("x", f"{x}")
        txt.set("y", f"{y}")
        txt.set("font-family", self.FONT_STACK)
        txt.set("font-size", f"{data['font_size'] * 0.3528}")
        if data["font_weight"] == "bold":
            txt.set("font-weight", "bold")
        txt.set("text-anchor",
                {"start": "start", "middle": "middle",
                 "end": "end"}.get(data.get("anchor", "start"), "start"))
        txt.set("fill", data.get("color") or "black")
        if data.get("rotation"):
            txt.set("transform", f"rotate({-data['rotation']}, {x}, {y})")
        txt.text = data["text"]


# ---------------------------------------------------------------------------
# PROVENANCE - the only authority on which drawing carries which letter
# ---------------------------------------------------------------------------

def provenance_main():
    """{figure number -> {printed letter -> [(source_dir, filename), ...]}}.

    Rows whose `printed_panel` is empty are directories that map to no printed
    panel at all - Figure 2's 02_C1, 02_C2, 02_E2 and 02_P. They are NOT in the
    published figure, so putting them into Version B would change what the
    figure shows, which the redraw is not allowed to do. They are returned
    separately so the run can report them.
    """
    letters, orphans = {}, []
    with open(PROVENANCE, newline="") as fh:
        for row in csv.DictReader(fh):
            fig = (row["figure"] or "").strip()
            if not fig or fig.startswith("S"):
                continue
            panel = (row["printed_panel"] or "").strip()
            dirs = [d.strip() for d in (row["source_dir"] or "").split(";")
                    if d.strip()]
            files = [f.strip() for f in (row["source_file"] or "").split(";")
                     if f.strip()]
            if not panel:
                orphans.append((fig, dirs, (row["note"] or "").strip()))
                continue
            if not files:
                letters.setdefault(fig, {})[panel] = []
                continue
            pairs = []
            for f in files:
                p = Path(f)
                pairs.append((p.parent.name or (dirs[0] if dirs else ""),
                              p.name))
            letters.setdefault(fig, {})[panel] = pairs
    return letters, orphans


def printed_rects():
    """{(figure, printed letter) -> (w_mm, h_mm)} off the published page."""
    out = {}
    with open(RECTS, newline="") as fh:
        for row in csv.DictReader(fh):
            fig = row["figure"].strip()
            if not fig.startswith("Figure "):
                continue
            try:
                x0, y0, x1, y1 = [float(v) for v in row["rect_mm"].split(",")]
            except Exception:
                continue
            out[(fig.split()[1], row["printed_panel"].strip())] = (x1 - x0,
                                                                  y1 - y0)
    return out


def check_layout_against_provenance(letters):
    """Every letter declared OMITTED or SHARED must be one PROVENANCE prints.

    The rows themselves are derived from PROVENANCE by `pack_rows()`, so they
    cannot disagree with it. What CAN disagree are the two hand-written tables
    above, and this is where they are checked. CLAUDE.md rule 2: PROVENANCE is
    the authority and a hand-written letter is a suspect.
    """
    problems = []
    for (figure, L) in OMITTED:
        n = figure.split("_")[1]
        if L not in letters.get(n, {}):
            problems.append(f"{figure}: {L} is declared OMITTED but "
                            f"PROVENANCE does not list it as a printed panel")
    for (figure, L), group in SHARED.items():
        n = figure.split("_")[1]
        for g in group:
            if g not in letters.get(n, {}):
                problems.append(f"{figure}: {g} is declared SHARED but "
                                f"PROVENANCE does not print it")
        files = {tuple(letters.get(n, {}).get(g, [])) for g in group}
        if len(files) != 1:
            problems.append(
                f"{figure}: {'/'.join(group)} are declared to share one "
                f"drawing, but PROVENANCE gives them different drawings: "
                f"{files}")
    if problems:
        raise AssertionError(
            "a hand-written table disagrees with PROVENANCE.csv:\n  "
            + "\n  ".join(problems))


def figures_from_provenance(letters):
    """{Figure_N: [printed letters in order]}, straight out of PROVENANCE."""
    out = OrderedDict()
    for n in sorted(letters, key=lambda s: int(s)):
        ls = sorted(letters[n])
        if ls:
            out[f"Figure_{n}"] = ls
    return out


def natural_size_mm(svg_path):
    """The size the panel was drawn at, in mm, from its own viewBox.

    This is the number the whole restyle rests on: a panel drawn at its printed
    size has a viewBox equal to that size, and `check_scale()` refuses to
    assemble if a placement does not match it.
    """
    head = svg_path.read_text(errors="replace")[:2000]
    m = re.search(r'viewBox="[\d.eE+-]+ [\d.eE+-]+ ([\d.eE+-]+) ([\d.eE+-]+)"',
                  head)
    if m:
        return float(m.group(1)) * MM, float(m.group(2)) * MM
    w = re.search(r'width="([\d.]+)pt"', head)
    h = re.search(r'height="([\d.]+)pt"', head)
    if not (w and h):
        raise ValueError(f"cannot read a size out of {svg_path}")
    return float(w.group(1)) * MM, float(h.group(1)) * MM


def figure_dir(figure):
    n = figure.split("_")[1]
    return MAIN / f"0{n}_Figure_{n}"


def drawing(figure, source_dir, filename):
    """The restyled SVG for one Version A drawing.

    Resolved by the directory and stem PROVENANCE names, never by scanning -
    three different Figure 3 panels write a file called
    `ceacam_spatial_boxplot.png`, in 03_H, 03_I and 03_J, and only the directory
    tells them apart.
    """
    d = figure_dir(figure) / source_dir
    svg = d / (Path(filename).stem + ".svg")
    if svg.exists():
        return svg
    # PROVENANCE and the script can disagree about the FILENAME while agreeing
    # about the panel. Figure 2's printed M is the one case: PROVENANCE records
    # `02_L/ihc_representative_2x3.svg`, the script is called
    # `create_ihc_representative_2x2.py`, it draws a 1x6 layout, and it writes
    # `ihc_representative_1x6.*`. Three older layouts sit beside it in Version
    # A. Version B keeps Version A's own output stem, as PANEL_SPEC requires,
    # so the stem does not match.
    #
    # Falling back is allowed ONLY when the directory holds exactly one SVG -
    # any ambiguity and this must stop, because picking the wrong one is
    # precisely how a panel ends up under the wrong letter. The mismatch is
    # printed on every run rather than silently absorbed; PROVENANCE is not
    # edited here.
    hits = sorted(p for p in d.glob("*.svg") if p.parent == d)
    if len(hits) == 1:
        print(f"    NOTE {figure} {source_dir}: PROVENANCE names "
              f"'{Path(filename).stem}.svg' but the only SVG in the directory "
              f"is '{hits[0].name}'. Using it, and reporting the mismatch - "
              f"PROVENANCE is not edited by this script.")
        return hits[0]
    raise FileNotFoundError(
        f"no restyled SVG for {figure} {source_dir}/{filename} at {svg}"
        + (f"\n    {len(hits)} other SVGs are in that directory "
           f"({[p.name for p in hits]}); refusing to guess which one is the "
           f"panel." if hits else
           "\n    Run the restyled panel script that produces it before "
           "assembling. Nothing is carried over from Version A."))


def drawings_in_order(figure, ordered_letters, letters):
    """Every drawing this figure needs, in printed-letter order.

    Returns a list of dicts: the SVG, its size, and the printed letter(s) it
    carries. A letter supplied by several drawings (Figure 5's J is four CD274
    boxplots, its N four GSEA dotplots) puts the letter on the FIRST of them and
    leaves the rest unlettered - the paper letters the group once. A drawing
    that carries several letters (Figure 2's H and I are the two halves of one
    SVG) carries them both.
    """
    n = figure.split("_")[1]
    by_letter = letters.get(n, {})
    out, consumed, placed_files = [], set(), set()
    for L in ordered_letters:
        if L in consumed or (figure, L) in OMITTED:
            continue
        group = SHARED.get((figure, L), (L,))
        consumed.update(group)
        files, seen = [], set()
        for g in group:
            for pair in by_letter.get(g, []):
                if pair not in seen:
                    seen.add(pair)
                    files.append(pair)
        if not files:
            raise KeyError(f"{figure}: PROVENANCE names no drawing for "
                           f"printed panel {'/'.join(group)}")
        for i, (sd, fn) in enumerate(files):
            if (sd, fn) in placed_files:
                continue
            placed_files.add((sd, fn))
            svg = drawing(figure, sd, fn)
            w, h = natural_size_mm(svg)
            out.append(dict(panel=f"{sd}/{Path(fn).stem}", svg=svg, w=w, h=h,
                            letters=list(group) if i == 0 else []))
    return out


def pack_rows(items):
    """Greedy rows, in printed-letter order, at the panels' natural sizes.

    Reading order always: a figure's letters have to run left to right and top
    to bottom, so a drawing is never moved past a later one to save height.
    That constraint is what makes these pages tall - two panels each just over
    half the page width cannot share a row - and it is the right trade, because
    a figure whose letters do not read in order is simply wrong.

    Drawings that share one printed letter DO wrap across rows. Figure 5's
    panel N is four GSEA dotplots at 88 mm each: 364 mm side by side, which no
    190.5 mm page can hold. They tile two-and-two and the letter N sits on the
    first.
    """
    for it in items:
        if it["w"] > PAGE_W + 0.01:
            raise ValueError(
                f"{it['panel']} is {it['w']:.1f} mm wide, wider than the "
                f"{PAGE_W:.1f} mm page on its own. Redraw it narrower - do not "
                f"scale it down, that is what Version A did.")

    # Drawings that share a printed letter are packed as a UNIT, so the four
    # CD274 boxplots the paper letters J stay together instead of being sliced
    # by whatever happened to precede them. The first attempt did slice them:
    # Figure 2's K came out with 02_J at the end of one row and 02_K at the
    # start of the next, which reads as two unrelated panels.
    groups, cur = [], []
    for it in items:
        if it["letters"] and cur:
            groups.append(cur)
            cur = []
        cur.append(it)
    if cur:
        groups.append(cur)

    def widths(seq):
        return sum(i["w"] for i in seq) + GUTTER * (len(seq) - 1)

    rows, row, width = [], [], 0.0
    for g in groups:
        gw = widths(g)
        if gw <= PAGE_W + 0.01:
            # the whole group fits one row: keep it whole
            if row and width + GUTTER + gw > PAGE_W + 0.01:
                rows.append(row)
                row, width = [], 0.0
            row.extend(g)
            width = widths(row)
            continue
        # too wide for any row: give the group its own rows and tile it
        if row:
            rows.append(row)
            row, width = [], 0.0
        sub, sw = [], 0.0
        for it in g:
            if sub and sw + GUTTER + it["w"] > PAGE_W + 0.01:
                rows.append(sub)
                sub, sw = [], 0.0
            sub.append(it)
            sw = widths(sub)
        if sub:
            rows.append(sub)
    if row:
        rows.append(row)
    return rows


def layout(rows):
    """Positions for every drawing, at natural size. (placements, page_h)."""
    placements, y = [], TOP
    for row in rows:
        row_top = y
        y += LETTER_BAND
        total_w = sum(c["w"] for c in row) + GUTTER * (len(row) - 1)
        x = (PAGE_W - total_w) / 2
        for c in row:
            c.update(x=x, y=y, letter_y=row_top)
            placements.append(c)
            x += c["w"] + GUTTER
        y += max(c["h"] for c in row) + ROW_GAP
    return placements, y - ROW_GAP + BOTTOM


def check_scale(placements):
    """Every panel must be placed at exactly the size it was drawn at."""
    bad = []
    for p in placements:
        w, h = natural_size_mm(p["svg"])
        if abs(w - p["w"]) > 0.05 or abs(h - p["h"]) > 0.05:
            bad.append(f"{p['panel']}: drawn {w:.2f} x {h:.2f} mm, "
                       f"placed {p['w']:.2f} x {p['h']:.2f} mm")
    if bad:
        raise AssertionError(
            "panels would be rescaled at assembly, which is the fault the "
            "restyle exists to remove:\n  " + "\n  ".join(bad))


def fit_report(page_w, page_h):
    scale = min(1.0, PRINTABLE_W / page_w, PRINTABLE_H / page_h)
    return scale, style.tick_pt() * scale


def _selected(argv, figures):
    want = [a for a in argv if not a.startswith("-")]
    if not want:
        return figures
    return OrderedDict((k, v) for k, v in figures.items()
                       if k in want or any(k.startswith(w) for w in want))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    style.apply()                       # before letter_font(); see stage 2
    family, weight = style.letter_font()
    letters, orphans = provenance_main()
    check_layout_against_provenance(letters)
    rects = printed_rects()

    print(f"page {PAGE_W:.1f} mm wide; panel letters {style.letter_pt():g} pt "
          f"{weight} {family}")
    print(f"letters looked up in {PROVENANCE.relative_to(ROOT)} - never "
          f"inferred from a directory name\n")

    written = []
    figures = figures_from_provenance(letters)
    for figure, ordered in _selected(sys.argv[1:], figures).items():
        n = figure.split("_")[1]
        if figure in NO_CODE:
            continue
        items = drawings_in_order(figure, ordered, letters)
        placements, page_h = layout(pack_rows(items))
        check_scale(placements)
        asm = RestyledAssembler(PAGE_W, page_h, title=None)
        print(f"{figure}   {PAGE_W:.1f} x {page_h:.1f} mm")
        for p in placements:
            asm.place_panel(None, p["svg"], p["x"], p["y"], p["w"], p["h"])
            # Two letters on one drawing (Figure 2 H and I) sit over the two
            # halves they label; one letter sits at the drawing's left edge.
            k = len(p["letters"])
            for i, L in enumerate(p["letters"]):
                dx = 0.0 if k <= 1 else i * p["w"] / k
                asm.add_text(L, p["x"] + LETTER_DX + dx,
                             p["letter_y"] + LETTER_DY,
                             font_size=style.letter_pt(), font_weight=weight)
            pr = rects.get((n, p["letters"][0] if p["letters"] else ""), None)
            was = f"  (printed {pr[0]:.1f} x {pr[1]:.1f})" if pr else ""
            print(f"    {''.join(p['letters']) or '-':<3} {p['panel']:<38}"
                  f"{p['w']:6.1f} x {p['h']:5.1f} mm at "
                  f"({p['x']:5.1f}, {p['y']:5.1f})  scale 1.000{was}")
        scale, smallest = fit_report(PAGE_W, page_h)
        print(f"    page fits {PRINTABLE_W:.0f} x {PRINTABLE_H:.0f} mm at "
              f"{scale:.3f}; smallest type would then print at "
              f"{smallest:.2f} pt")
        asm.save(figure, output_dir=OUT)
        written.append((figure, page_h, scale, smallest))

    if orphans:
        print(f"\n{len(orphans)} panel director"
              f"{'y' if len(orphans) == 1 else 'ies'} map to no printed panel "
              f"and are therefore NOT in Version B - drawing them would change "
              f"what the figure shows:")
        for fig, dirs, note in orphans:
            print(f"   Figure {fig}  {', '.join(dirs)}"
                  f"{'  - ' + note if note else ''}")
    if OMITTED:
        print(f"\n{len(OMITTED)} printed panel(s) are NOT in Version B "
              f"because they could not be gated:")
        for (fig, L), why in sorted(OMITTED.items()):
            print(f"   {fig} panel {L}: {why}")
        print("   Every reference to those letters would be left pointing at "
              "nothing. See ../../RENUMBERING.md.")
    for fig, why in NO_CODE.items():
        print(f"\n{fig} is not rebuilt: {why}. Version A's {fig} stands.")

    print(f"\n{'figure':<14}{'page mm':>16}{'journal fit':>13}"
          f"{'smallest pt':>13}")
    for figure, page_h, scale, smallest in written:
        print(f"{figure:<14}{PAGE_W:6.1f} x {page_h:6.1f}{scale:13.3f}"
              f"{smallest:13.2f}")
    print(f"\nAt 1:1 the smallest body type on every page above is "
          f"{style.tick_pt():g} pt. The journal-fit column is reported, not "
          f"enforced; see ../../RESTYLE_REPORT.md.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
