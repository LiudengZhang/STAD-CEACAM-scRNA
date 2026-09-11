#!/usr/bin/env python3
"""
Draw S9A and S9B - monocyte versus macrophage lineage (Reviewer 1, R1.7).

Analysis: 04_Revision_Analyses/06_R1.7_MoMac_Lineage_Markers/scripts/momac_lineage.py

NOTHING IS RECOMPUTED HERE.

S9A reads `momac_lineage_scores.csv`, the per-state summary `main()` already
wrote. `sc.tl.score_genes` is never called and MoMac.h5ad is never opened for
this panel.

S9B is the one panel in this driver that needs the matrix, because the original
`_panel_dotplot(ad, mono, mac)` was itself given the AnnData and `sc.pl.dotplot`
computes its own mean expression and fraction expressing from it. The object is
loaded exactly as the analysis loads it and the two gene lists are derived
exactly as `main()` derives them, from `set(ad.var_names)`. `sc.tl.score_genes`
is deliberately NOT run: the dotplot does not read `monocyte_score` or
`macrophage_score`, so running it here would be a recomputation with no
consumer.

S9B DOES NOT PRINT "CONTENT IDENTICAL", AND THAT IS NOT A CONTENT DIFFERENCE.
`--check S9_B` reports exactly one line:

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
    python draw_momac_lineage.py --check S9_A
"""

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
FIG = "S9_MoMac_Identity_NFkB"

# Printed boxes, millimetres.
#
# S9B is 21 gene columns plus scanpy's 1.5 in (38 mm) legend column plus the
# seven cell-state names on the left, so it takes a full-width row of its own.
# S9A is set at the width its own labels need (see `draw_scatter`) and sits on
# the row above it.
A_W, A_H = 115.0, 60.0
B_W, B_H = 171.0, 74.0

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
                   edgecolors="white", linewidths=0.5 * A.SCALE * A_MARK,
                   zorder=3)
        # C4 and C1 sit at almost the same height, and a label to the right of
        # C4 would end under C1's point and read as C1's. C4 is labelled on its
        # left instead, into empty space below the diagonal.
        label = A.shown(state).replace("_", " ")
        left = label.startswith("C4 ")
        ax.annotate(label, (r["monocyte_score"], r["macrophage_score"]),
                    textcoords="offset points",
                    xytext=((-6 if left else 6) * A.SCALE * A_MARK,
                            3 * A.SCALE * A_MARK),
                    ha="right" if left else "left",
                    fontsize=style.tick_pt(),
                    color="#B2182B" if highlight else "#333333")
    lim = [min(summary["monocyte_score"].min(),
               summary["macrophage_score"].min()),
           max(summary["monocyte_score"].max(),
               summary["macrophage_score"].max())]
    pad = 0.12 * (lim[1] - lim[0])
    ax.plot([lim[0] - pad, lim[1] + pad], [lim[0] - pad, lim[1] + pad],
            color="#bbbbbb", linestyle="--", linewidth=0.5, zorder=1)
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
    if save:
        base.save(fig, FIG, "S9_A", "S9_A_momac_lineage_scores")
    return fig


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


def draw_dotplot(save=True):
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
    base.fit(fig)
    if save:
        base.save(fig, FIG, "S9_B", "S9_B_momac_lineage_dotplot")
    return fig


def main():
    fr = frames()
    return base.run({
        # `obs` is unused by _panel_scatter - the function reads only
        # `summary`, so None is passed and the original runs identically.
        "S9_A": (lambda: draw_scatter(fr),
                 lambda: A._panel_scatter(fr["summary"], None)),
        "S9_B": (lambda: draw_dotplot(),
                 lambda: A._panel_dotplot(matrix()["ad"], matrix()["mono"],
                                          matrix()["mac"])),
    })


if __name__ == "__main__":
    sys.exit(main())
