#!/usr/bin/env python3
"""
Restyle S9C and S9D - monocyte versus macrophage lineage (Reviewer 1, R1.7).

Analysis: 02_New_Analyses/06_R1.7_MoMac_Lineage_Markers/scripts/momac_lineage.py

NOTHING IS RECOMPUTED HERE.

S9C reads `momac_lineage_scores.csv`, the per-state summary `main()` already
wrote. `sc.tl.score_genes` is never called and MoMac.h5ad is never opened for
this panel.

S9D is the one panel in this driver that needs the matrix, because the original
`_panel_dotplot(ad, mono, mac)` was itself given the AnnData and `sc.pl.dotplot`
computes its own mean-expression and fraction-expressing from it. The object is
loaded exactly as the analysis loads it (`sc.read_h5ad(A.MOMAC_H5AD)`) and the
two gene lists are derived exactly as `main()` derives them, from
`set(ad.var_names)`. `sc.tl.score_genes` is deliberately NOT run: the dotplot
does not read `monocyte_score` or `macrophage_score`, so running it here would
be a recomputation with no consumer.

S9D DOES NOT PRINT "CONTENT IDENTICAL", AND THAT IS NOT A CONTENT DIFFERENCE.
`--check S9_D` reports exactly one line:

    figures: length 3 vs 1

`_panel_dotplot` saves three formats by calling `BasePlot.savefig` three times,
and `BasePlot.savefig` calls `make_figure()` every time - so the original builds
the same figure three times over and the capture harness collects all three. The
restyle builds it once and writes the three formats off that one figure.

Measured rather than argued (all five axes, to 9 dp):

    original fig0 vs fig1        0 differences
    original fig0 vs fig2        0 differences
    original fig0 vs restyled    0 content, 0 tick
    original fig1 vs restyled    0 content, 0 tick
    original fig2 vs restyled    0 content, 0 tick

Every dot offset, every mean-expression value, every fraction, every gene and
category label and both var-group brackets are identical. The count is a
property of Version A's save loop, not of what is drawn. It is left standing
rather than papered over: making it "pass" would mean either drawing the panel
three times or editing the harness, and neither is a restyle.

    conda run -n Liudeng_Python_310 python restyle_momac_lineage.py
    conda run -n Liudeng_Python_310 python restyle_momac_lineage.py --check
    conda run -n Liudeng_Python_310 python restyle_momac_lineage.py --check S9_C
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import scanpy as sc
import cnsplots as cns

import _driver_base as base
import panel_style_cns as style

A = base.analysis("06_R1.7_MoMac_Lineage_Markers/scripts/momac_lineage.py")
OUT = base.outputs_of(A)
FIG = "S9_Mechanism_Specificity"

# Printed boxes, millimetres. Version A: S9C 84 x 62, S9D 81 x 62, side by side.
#
# S9D cannot stay at 81 mm. It is 21 gene columns plus scanpy's 1.5 in (38 mm)
# legend column plus the seven cell-state names on the left; at 7 pt the names
# alone need 42 mm, which would leave 1 mm for the 21 columns. S9D therefore
# takes a full-width row of its own, and S9C, which no longer has to share it,
# is set at the width its own labels need (see `draw_scatter`).
C_W, C_H = 115.0, 62.0
D_W, D_H = 171.0, 76.0

# Mark sizes for S9C, rescaled so they keep their size relative to the type.
# Version A drew it at 4x and the assembler fitted it at 0.3100; its smallest
# type was the 4.5 pt state label, i.e. 5.58 pt on paper.
C_FIT = 0.3100
C_TYPE = 4.5 * A.SCALE * C_FIT
C_MARK = (style.tick_pt() / C_TYPE) * C_FIT
C_AREA = C_MARK ** 2


def frames():
    """The summary table the analysis already wrote. No computation."""
    summary = pd.read_csv(base.require(
        OUT / "momac_lineage_scores.csv", "S9C lineage scores"), index_col=0)
    for col in ("monocyte_score", "macrophage_score"):
        if col not in summary.columns:
            raise RuntimeError(f"momac_lineage_scores.csv has no `{col}` "
                               f"column: columns are {list(summary.columns)}")
    return dict(summary=summary)


_MATRIX = {}


def matrix():
    """MoMac.h5ad and the two gene panels, exactly as momac_lineage.main() has
    them. Loaded once, and only when S9D is actually drawn."""
    if not _MATRIX:
        ad = sc.read_h5ad(base.require(A.MOMAC_H5AD, "S9D MoMac matrix"))
        present = set(ad.var_names)                       # momac_lineage.py:78
        _MATRIX["ad"] = ad
        _MATRIX["mono"] = [g for g in A.MONOCYTE if g in present]
        _MATRIX["mac"] = [g for g in A.MACROPHAGE if g in present]
    return _MATRIX


# --------------------------------------------------------------------------
# S9C - state means on the monocyte/macrophage score plane. Layout unchanged.
#
# The label offsets are rescaled with the type like any other length in points:
# (24, 12) pt at 4x becomes (9.3, 4.7) pt here. compare_panel_content records an
# Annotation by the point it points AT, in data coordinates, and treats the
# offset as layout, so this is a restyle and not a moved value.
#
# What sets C_W is the label run-out. xlim is content - the analysis pads the
# right edge by 5.5 * pad precisely so the labels have somewhere to go - and
# C1_Mono_Classic_CD14 sits at 63% of that span with a 20-character label to
# its right, so the axes need about 76 mm for it to fit inside the canvas. The
# binding label is in fact C4, set ha="right" and running back toward the y
# axis: at 100 mm its first character touched the canvas edge.
# --------------------------------------------------------------------------
def draw_scatter(fr, save=True):
    style.apply()
    summary = fr["summary"]
    fig, ax = style.subplots_mm(C_W, C_H)
    for state, r in summary.iterrows():
        highlight = state == A.TARGET
        ax.scatter(r["monocyte_score"], r["macrophage_score"],
                   s=(70 if highlight else 45) * A.SCALE * C_AREA,
                   c="#B2182B" if highlight else "#4d4d4d",
                   edgecolors="white", linewidths=0.5 * A.SCALE * C_MARK,
                   zorder=3)
        # C4 and C1 sit at almost the same height, and a label to the right of
        # C4 would end under C1's point and read as C1's. C4 is labelled on its
        # left instead, into empty space below the diagonal.
        label = A.shown(state).replace("_", " ")
        left = label.startswith("C4 ")
        ax.annotate(label, (r["monocyte_score"], r["macrophage_score"]),
                    textcoords="offset points",
                    xytext=((-6 if left else 6) * A.SCALE * C_MARK, 3 * A.SCALE * C_MARK),
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
    style.margins_mm(fig, left=14, right=3, top=4, bottom=12)
    if save:
        base.save(fig, FIG, "S9_C", "S9_C_momac_lineage_scores")
    return fig


# --------------------------------------------------------------------------
# S9D - marker dotplot. scanpy builds this figure, not matplotlib: the box is
# given to sc.pl.dotplot as inches, and the margins are set afterwards on the
# figure scanpy made (its outer GridSpec carries no left/right of its own, so
# subplots_adjust still reaches it).
# --------------------------------------------------------------------------
def draw_dotplot(save=True):
    m = matrix()
    # style.apply() first, for the exact-size save and the font stack - both
    # live on cns.settings, which cns.setup_scanpy() then reads. setup_scanpy
    # is what configures scanpy itself and is the documented entry point for
    # a scanpy-built figure.
    style.apply()
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
        figsize=style.figsize_mm(D_W, D_H))
    dp.make_figure()
    fig = dp.fig
    style.margins_mm(fig, left=42, right=3, top=22, bottom=13)
    if save:
        base.save(fig, FIG, "S9_D", "S9_D_momac_lineage_dotplot")
    return fig


def main():
    fr = frames()
    return base.run({
        # `obs` is unused by _panel_scatter - the function reads only
        # `summary`, so None is passed and the original runs identically.
        "S9_C": (lambda: draw_scatter(fr),
                 lambda: A._panel_scatter(fr["summary"], None)),
        "S9_D": (lambda: draw_dotplot(),
                 lambda: A._panel_dotplot(matrix()["ad"], matrix()["mono"],
                                          matrix()["mac"])),
    })


if __name__ == "__main__":
    sys.exit(main())
