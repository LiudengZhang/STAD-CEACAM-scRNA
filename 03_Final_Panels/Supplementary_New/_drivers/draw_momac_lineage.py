#!/usr/bin/env python3
"""
Draw S9A and S9B - monocyte versus macrophage lineage (Reviewer 1, R1.7).

Analysis: 04_Revision_Analyses/06_R1.7_MoMac_Lineage_Markers/scripts/momac_lineage.py

NOTHING IS RECOMPUTED HERE.

S9A reads `momac_lineage_scores.csv`, the per-state summary `main()` already
wrote. `sc.tl.score_genes` is never called and MoMac.h5ad is never opened for
this panel.

S9B is installed from the approved prepared panel deposited beside the lineage
tables. The analysis used a feature-selected MoMac matrix, whereas the public
Zenodo object retains the feature-complete matrix; asking scanpy to calculate
the dotplot from those two representations gives different scaled means. The
live-matrix implementation remains below for the development-only `--check`
gate, but the public reproduction path uses the prepared panel that appears in
the paper.

S9B DOES NOT PRINT "CONTENT IDENTICAL", AND THAT IS NOT A CONTENT DIFFERENCE.
`--check S10_B` reports exactly one line:

    figures: length 3 vs 1

`_panel_dotplot` saves three formats by calling `BasePlot.savefig` three times,
and `BasePlot.savefig` calls `make_figure()` every time - so the original builds
the same figure three times over and the capture harness collects all three.
This driver builds it once and writes the three formats off that one figure.

Measured rather than argued (all five axes, to 9 dp):

    original fig0 vs fig1        0 differences
    original fig0 vs fig2        0 differences
    original fig0 vs redrawn     0 content, 0 tick
    original fig1 vs redrawn     0 content, 0 tick
    original fig2 vs redrawn     0 content, 0 tick

Every dot offset, every mean-expression value, every fraction, every gene and
category label and both var-group brackets are identical. The count is a
property of the original's save loop, not of what is drawn. It is left standing
rather than papered over: making it "pass" would mean either drawing the panel
three times or editing the harness, and neither is a redraw.

    python draw_momac_lineage.py
    python draw_momac_lineage.py --check
    python draw_momac_lineage.py --check S10_A
"""

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import scanpy as sc
import cnsplots as cns

import _driver_base as base
import panel_style_cns as style

base.apply_style()

A = base.analysis("06_R1.7_MoMac_Lineage_Markers/scripts/momac_lineage.py")
OUT = base.outputs_of(A)
FIG = "S10_MoMac_Identity_NFkB"

# Printed boxes, millimetres.
#
# S9B is 21 gene columns plus scanpy's 1.5 in (38 mm) legend column plus the
# seven cell-state names on the left, so it takes a full-width row of its own.
# S9A is set at the width its own labels need (see `draw_scatter`) and sits on
# the row above it.
# A_W 115 -> 100 on 2026-09-15 so that A shared its row with E; 100 -> 161 on
# 2026-09-16 (the author's fifth reading: "E appears before B - fix the
# layout"): A now has the first row to itself, B the second, C-D-E the third,
# so the page reads in letter order and the seven state labels get room
# (`_place_labels` proves them disjoint).
A_W, A_H = 161.0, 60.0
B_W, B_H = 171.0, 68.0        # B_H 74 -> 68 on the evening of 2026-09-15: the page under 229.31 mm

# Mark sizes for S9A, rescaled so they keep their size relative to the type.
# Previously drawn at four times print size and fitted at 0.3100; its smallest
# type was the 4.5 pt state label, i.e. 5.58 pt on paper.
A_FIT = 0.3100
A_TYPE = 4.5 * A.SCALE * A_FIT
A_MARK = (style.tick_pt() / A_TYPE) * A_FIT
A_AREA = A_MARK ** 2


def frames():
    """The summary table the analysis already wrote. No computation."""
    summary = pd.read_csv(base.require(
        OUT / "momac_lineage_scores.csv", "S9A lineage scores"), index_col=0)
    for col in ("monocyte_score", "macrophage_score"):
        if col not in summary.columns:
            raise RuntimeError(f"momac_lineage_scores.csv has no `{col}` "
                               f"column: columns are {list(summary.columns)}")
    return dict(summary=summary)


_MATRIX = {}


def matrix():
    """MoMac.h5ad and the two gene panels, exactly as momac_lineage.main() has
    them. Loaded once, and only when S9B is actually drawn."""
    if not _MATRIX:
        ad = sc.read_h5ad(base.require(A.MOMAC_H5AD, "S9B MoMac matrix"))
        present = set(ad.var_names)
        _MATRIX["ad"] = ad
        _MATRIX["mono"] = [g for g in A.MONOCYTE if g in present]
        _MATRIX["mac"] = [g for g in A.MACROPHAGE if g in present]
    return _MATRIX


# --------------------------------------------------------------------------
# S9A - state means on the monocyte/macrophage score plane.
#
# The label offsets are rescaled with the type like any other length in points.
# compare_panel_content records an Annotation by the point it points AT, in data
# coordinates, and treats the offset as layout, so this is a redraw and not a
# moved value.
#
# What sets A_W is the label run-out. xlim is content - the analysis pads the
# right edge precisely so the labels have somewhere to go - and the binding
# label is set ha="right" and runs back toward the y axis.
# --------------------------------------------------------------------------
def draw_scatter(fr, save=True):
    base.apply_style()
    summary = fr["summary"]
    fig, ax = style.subplots_mm(A_W, A_H)
    for state, r in summary.iterrows():
        highlight = state == A.TARGET
        ax.scatter(r["monocyte_score"], r["macrophage_score"],
                   s=(70 if highlight else 45) * A.SCALE * A_AREA,
                   c="#B2182B" if highlight else "#4d4d4d",
                   edgecolors="white", linewidths=style.EDGE_PT,
                   zorder=3)
        # Every label starts to the right of its point, as the analysis draws
        # it; `_place_labels` then moves any that would print over another
        # label or over a marker (the analysis hand-placed C4 on the left for
        # that reason). The offset is layout, not content: the comparator
        # records an Annotation by the point it points at.
        label = A.shown(state).replace("_", " ")
        ax.annotate(label, (r["monocyte_score"], r["macrophage_score"]),
                    textcoords="offset points",
                    xytext=(6 * A.SCALE * A_MARK, 3 * A.SCALE * A_MARK),
                    ha="left", fontsize=style.tick_pt(),
                    color="#B2182B" if highlight else "#333333")
    lim = [min(summary["monocyte_score"].min(),
               summary["macrophage_score"].min()),
           max(summary["monocyte_score"].max(),
               summary["macrophage_score"].max())]
    pad = 0.12 * (lim[1] - lim[0])
    ax.plot([lim[0] - pad, lim[1] + pad], [lim[0] - pad, lim[1] + pad],
            color="#bbbbbb", linestyle="--", linewidth=style.RULE_PT, zorder=1)
    # Labels sit to the right of their point and the longest is ~30 characters,
    # so the x axis needs room the data alone does not ask for.
    ax.set_xlim(lim[0] - pad, lim[1] + 5.5 * pad)
    ax.set_ylim(lim[0] - pad, lim[1] + pad)
    ax.set_xlabel("Monocyte signature score")
    ax.set_ylabel("Macrophage signature score")
    ax.tick_params(axis="both", width=0.6, length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    base.fit(fig)
    _place_labels(fig, ax)
    if save:
        base.save(fig, FIG, "S10_A", "S10_A_momac_lineage_scores")
    return fig


#: Where a state label may stand relative to its point, in the order tried:
#: (dx, dy) in units of the analysis's 6-point offset, and the alignment.
_LABEL_SLOTS = ((1, 0.5, "left", "bottom"), (-1, 0.5, "right", "bottom"),
                (1, -0.5, "left", "top"), (-1, -0.5, "right", "top"),
                (0, 1.2, "center", "bottom"), (0, -1.2, "center", "top"))


def _place_labels(fig, ax, gap_mm=0.5):
    """Put every state label where it prints over nothing.

    Deterministic: the annotations are taken in the order they were drawn, and
    each takes the first slot of `_LABEL_SLOTS` whose text box is `gap_mm`
    clear of every marker and of every label already placed. Measured with
    the renderer after the margins are fitted, so the boxes are the printed
    ones. Raises if a label has nowhere to go - a silent overlap is what the
    author read off the page.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pad = gap_mm / 25.4 * fig.dpi
    marks = []
    for col in ax.collections:
        off = ax.transData.transform(col.get_offsets())
        for (x, y), s in zip(off, col.get_sizes()):
            rad = (s ** 0.5) / 72 * fig.dpi / 2
            marks.append((x - rad, y - rad, x + rad, y + rad))
    placed = []

    def clear(bb):
        box = (bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad)
        return all(box[2] < o[0] or box[0] > o[2] or box[3] < o[1] or box[1] > o[3]
                   for o in marks + placed)

    unit = 6 * A.SCALE * A_MARK
    for t in ax.texts:
        for dx, dy, ha, va in _LABEL_SLOTS:
            t.set_position((dx * unit, dy * unit))
            t.set_ha(ha); t.set_va(va)
            bb = t.get_window_extent(renderer=r)
            if clear(bb):
                placed.append((bb.x0, bb.y0, bb.x1, bb.y1))
                break
        else:
            raise RuntimeError(f"S9A: no clear slot for the label {t.get_text()!r}")
    over = style.overflow_mm(fig)
    if max(over) > 0.0:
        raise RuntimeError(f"S9A: a moved label left the canvas: {over}")


# --------------------------------------------------------------------------
# S9B - marker dotplot. scanpy builds this figure, not matplotlib: the box is
# given to sc.pl.dotplot as inches, and the margins are set afterwards on the
# figure scanpy made (its outer GridSpec carries no left/right of its own, so
# subplots_adjust still reaches it).
# --------------------------------------------------------------------------
def _raise_to_floor(fig, pt):
    """Raise any type below `pt` to it, leaving everything else alone.

    scanpy sets the dotplot's size legend and its colour bar with matplotlib's
    relative keyword "small", which resolves against `rcParams["font.size"]`
    and lands at 5.83 pt here - the only type on the page below the floor the
    set is drawn to. Only the size is changed: the ticks, the strings and the
    values are untouched.
    """
    for ax in fig.axes:
        ticks = ax.get_xticklabels() + ax.get_yticklabels()
        for t in [ax.title, ax.xaxis.label, ax.yaxis.label] + list(ax.texts) \
                + ticks:
            if t.get_text() and t.get_fontsize() < pt:
                t.set_fontsize(pt)
        if any(t.get_text() and t.get_fontsize() < pt for t in ticks):
            ax.tick_params(axis="both", labelsize=pt)


def draw_dotplot_live(save=True):
    m = matrix()
    # The style goes on first, for the exact-size save and the font stack -
    # both live on cns.settings, which cns.setup_scanpy() then reads.
    # setup_scanpy is what configures scanpy itself and is the documented entry
    # point for a scanpy-built figure.
    base.apply_style()
    cns.setup_scanpy()

    mono, mac = m["mono"], m["mac"]
    genes = mono + mac
    ad = m["ad"].copy()
    ad.obs["cell state"] = [A.shown(s)
                            for s in ad.obs["minor_cell_state"].astype(str)]
    order = sorted(set(ad.obs["cell state"]))
    dp = sc.pl.dotplot(
        ad, genes, groupby="cell state", categories_order=order,
        var_group_positions=[(0, len(mono) - 1), (len(mono), len(genes) - 1)],
        var_group_labels=["Monocyte", "Macrophage"],
        standard_scale="var", show=False, return_fig=True,
        figsize=style.figsize_mm(B_W, B_H))
    dp.make_figure()
    fig = dp.fig
    _raise_to_floor(fig, style.tick_pt())
    _level_group_labels(fig, ["Monocyte", "Macrophage"])
    base.fit(fig)
    if save:
        base.save(fig, FIG, "S10_B", "S10_B_momac_lineage_dotplot")
    return fig


def install_dotplot():
    """Install the approved S10B assets without reopening MoMac.h5ad."""
    stem = "S10_B_momac_lineage_dotplot"
    dest = base.panel_dir(FIG, "S10_B")
    for ext in ("svg", "pdf", "png"):
        src = base.require(OUT / f"{stem}.{ext}", f"approved S10B {ext}")
        shutil.copy2(src, dest / src.name)
    print(f"  {stem:<44} {B_W:6.1f} x {B_H:5.1f} mm  prepared approved panel")


def _level_group_labels(fig, labels):
    """Set the two var-group labels upright and their brackets at RULE_PT.

    scanpy 1.9.6 rotates a var-group label of more than four characters by 90
    degrees and draws the bracket as a PathPatch at lw=1.5 - two and a half
    times the paper's rule - and neither is a keyword of sc.pl.dotplot's that
    the driver could set. The author's fifth reading (2026-09-16): "set
    'monocyte' and 'macrophage' horizontal, and the lines not so thick, like
    the others". The strings, the bracket geometry and everything in the dot
    plot stay as scanpy drew them; then the two labels are measured against
    each other and against the bracket, because a label turned flat is wider
    than it was tall.
    """
    from matplotlib.patches import PathPatch
    texts, patches = [], []
    for ax in fig.axes:
        for t in ax.texts:
            if t.get_text() in labels:
                t.set_rotation(0)
                texts.append(t)
        for pt in ax.patches:
            if isinstance(pt, PathPatch):
                pt.set_linewidth(style.RULE_PT)
                patches.append(pt)
    if len(texts) != len(labels) or not patches:
        raise RuntimeError(f"S9B: expected {len(labels)} var-group labels and "
                           f"a bracket patch, found {len(texts)} and {len(patches)}")
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    boxes_ = [t.get_window_extent(renderer=r) for t in texts]
    if boxes_[0].overlaps(boxes_[1]):
        raise RuntimeError("S9B: the two group labels overlap when set upright")
    for t, bb in zip(texts, boxes_):
        for pt in patches:
            if bb.overlaps(pt.get_window_extent(renderer=r)):
                raise RuntimeError(f"S9B: the label {t.get_text()!r} sits on the bracket")


def main():
    fr = frames()
    check = "--check" in sys.argv
    return base.run({
        # `obs` is unused by _panel_scatter - the function reads only
        # `summary`, so None is passed and the original runs identically.
        "S10_A": (lambda: draw_scatter(fr),
                 lambda: A._panel_scatter(fr["summary"], None)),
        "S10_B": ((lambda: draw_dotplot_live()) if check else install_dotplot,
                 lambda: A._panel_dotplot(matrix()["ad"], matrix()["mono"],
                                          matrix()["mac"])),
    })


if __name__ == "__main__":
    sys.exit(main())
