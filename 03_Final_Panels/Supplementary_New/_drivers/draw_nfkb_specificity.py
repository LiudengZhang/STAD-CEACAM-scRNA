#!/usr/bin/env python3
"""
Draw S9C and S9D - NF-kB specificity and the epithelial compartment
(Reviewer 1, point R1.8).

Analysis: 04_Revision_Analyses/07_R1.8_NFkB_Specificity/scripts/nfkb_specificity.py

NOTHING IS RECOMPUTED HERE. In particular the driver never reaches the
Hallmark-table walk in `load_gsea()`, which re-reads the recompute tables,
re-ranks all fifty sets per file, and would silently drop a cell type through
its `if term is None: continue` if a table were renamed. A redraw must not be
able to change which cell types appear in S9C, so the driver reads a table off
disk rather than rebuilding one.

WHICH table is the whole question for S9C, and it is NOT this module's own
`outputs/nfkb_per_celltype.csv`. That file is the LIVE run, computed on the
doubly-normalised .X of 00_Data_Audit/FINDINGS.md sections 1 and 7. `main()`
writes it, but it does not draw the panel from it: the last lines of `main()`
read

    _panel_nfkb(load_adopted())

so the frame the panel is given is the ADOPTED sound-input table,
13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/nfkb_per_celltype_sound13.csv,
which the author ruled on 2026-09-03 is the one every NF-kB panel and every
NF-kB sentence must rest on (07_Archive/2026-09-03_S9E_S10C_drawn_from_live_run).
The two disagree in six quantities, listed at `ADOPTED_GSEA` in the analysis;
against the live table B cells print #37/38 instead of last, and five cell types
print rank 1 instead of three, which contradicts the printed Results sentence
and the Fig. S9C columns of Table S10. So `g` comes from `A.load_adopted()` -
the analysis module's own accessor for that file, which raises if it is missing.

This is what a driver reading `outputs/*.csv` by reflex cannot see, and the
`--check` gate could not catch it while the same wrong frame was handed to both
sides of the comparison. `main()` below now passes `fr["g"]` to
`A._panel_nfkb`, which is the adopted table, so the gate compares what `main()`
would actually draw.

`_panel_cytokines` is given `cyto` in the same way, from
`outputs/epithelial_cytokines.csv`. No h5ad is opened, so `sample_means()` and
its Mann-Whitney tests are not re-run either.

Both panel functions do their own selection inside the panel
(`method == PRIMARY and phase == "post"` for S9C, `phase == "Post"` for S9D);
those selections are copied verbatim below and `--check` proves they select the
same rows.

    python draw_nfkb_specificity.py --check
    python draw_nfkb_specificity.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

import _driver_base as base
import panel_style_cns as style

base.apply_style()

A = base.analysis("07_R1.8_NFkB_Specificity/scripts/nfkb_specificity.py")
OUT = base.outputs_of(A)
FIG = "S9_MoMac_Identity_NFkB"

# Printed boxes, millimetres. The two sit side by side on one row of the
# 183 mm page. S9C is the wider of the two because each bar carries an FDR and
# a Hallmark rank annotated to its right; S9D is four grouped bars.
C_W, C_H = 92.0, 62.0
D_W, D_H = 74.0, 46.0


def frames():
    """The two tables the panels are drawn from. No computation.

    `g` is the adopted sound-input table, read through the analysis module's
    own `load_adopted()` so that one file name is resolved in one place; it is
    the frame `main()` gives `_panel_nfkb`. `cyto` is this module's own written
    output, which is the frame `main()` gives `_panel_cytokines`.
    """
    base.require(A.ADOPTED_GSEA, "S9C NF-kB per cell type (adopted table)")
    return dict(
        g=A.load_adopted(),
        cyto=pd.read_csv(base.require(
            OUT / "epithelial_cytokines.csv", "S9D epithelial cytokines")),
    )


# --------------------------------------------------------------------------
# S9C - post-treatment NES per cell type, with FDR and Hallmark rank annotated.
# One horizontal bar per cell type, sorted by NES.
# --------------------------------------------------------------------------
def draw_C(fr, save=True):
    base.apply_style()
    g = fr["g"]
    post = (g[(g["method"] == A.PRIMARY) & (g["phase"] == "post")]
            .sort_values("nes"))
    fig, ax = style.subplots_mm(C_W, C_H)
    y = np.arange(len(post))
    colors = [A.COLOR_NR if q < 0.25 else "#cccccc" for q in post["fdr_q"]]
    ax.barh(y, post["nes"], color=colors, edgecolor="#444444", linewidth=0.4,
            height=0.65)
    # Annotations always sit to the right of the bar's far end, so the ones on
    # negative bars do not run into the cell-type labels on the axis.
    for i, (_, r) in enumerate(post.iterrows()):
        ax.text(max(r["nes"], 0.0) + 0.05, i,
                f"q={r['fdr_q']:.2f}  #{r['rank']}/{r['n_sets']}",
                va="center", ha="left",
                fontsize=style.tick_pt(), color="#333333")
    ax.axvline(0, color="#666666", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([A.LABELS[c] for c in post["cell_type"]])
    ax.set_xlabel("NES, TNF$\\alpha$ signalling via NF-$\\kappa$B\n"
                  "(post-treatment, non-responders vs responders)")
    ax.set_xlim(min(0, post["nes"].min()) - 0.3, post["nes"].max() + 1.9)
    ax.tick_params(axis="x", width=0.6, length=2)
    ax.tick_params(axis="y", length=0)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="#444444",
                             linewidth=0.4, label=l)
               for c, l in ((A.COLOR_NR, "FDR q < 0.25"),
                            ("#cccccc", "FDR q >= 0.25"))]
    ax.legend(handles=handles, loc="lower right", frameon=False,
              handlelength=1.0, handletextpad=0.4, labelspacing=0.25)
    base.fit(fig)
    if save:
        base.save(fig, FIG, "S9_C", "S9_C_nfkb_per_celltype")
    return fig


# --------------------------------------------------------------------------
# S9D - epithelial versus myeloid cytokine expression after treatment.
# Two grouped bars per gene on a symlog axis.
# --------------------------------------------------------------------------
def draw_D(fr, save=True):
    base.apply_style()
    post = fr["cyto"][fr["cyto"]["phase"] == "Post"]
    fig, ax = style.subplots_mm(D_W, D_H)
    genes = A.CYTOKINES
    x = np.arange(len(genes))
    w = 0.38
    for k, (comp, color) in enumerate((("Epithelial", "#4d4d4d"),
                                       ("Monocytes/Macrophages", "#B2182B"))):
        sub = post[post["compartment"] == comp].set_index("gene")
        vals = [sub["mean_NR"].get(gene, np.nan) for gene in genes]
        ax.bar(x + (k - 0.5) * w, vals, width=w, color=color, alpha=0.85,
               edgecolor="#333333", linewidth=0.4, label=comp)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, style="italic")
    ax.set_ylabel("Mean expression,\npost-treatment non-responders")
    ax.set_yscale("symlog", linthresh=0.01)
    # A symlog axis labels its decades in scientific notation, and mathtext
    # draws a superscript at 70% of its base, so the exponent of 10^-2 would
    # print at 4.2 pt - the only type on the page below the 6 pt floor. The
    # decades are 0.01, 0.1 and 1, which are shorter written out than as powers
    # and put every glyph on this axis at the tick size. The tick positions,
    # the scale and the values are untouched; only the notation changes.
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, _: "0" if v == 0 else f"{v:g}"))
    ax.tick_params(axis="both", width=0.6, length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="upper left",
              handlelength=1.0, handletextpad=0.4, labelspacing=0.25)
    base.fit(fig)
    if save:
        base.save(fig, FIG, "S9_D", "S9_D_epithelial_vs_myeloid_cytokines")
    return fig


def main():
    fr = frames()
    return base.run({
        # `fr["g"]` is the adopted table, so the right-hand side is exactly
        # what the last lines of the analysis's main() draw, `_panel_nfkb(
        # load_adopted())`. Handing this side the live table instead is what
        # let the gate report CONTENT IDENTICAL on 2026-09-08 while the panel
        # printed the live run's numbers; both sides were then wrong together.
        "S9_C": (lambda: draw_C(fr), lambda: A._panel_nfkb(fr["g"])),
        # The symlog decades are printed written out rather than as powers of
        # ten so that no glyph on that axis falls below the type floor; the
        # change is declared in 00_Config/shared/labels.py and the gate names
        # each substitution it accepts.
        # base.aliases, not base.aliases(): resolving the alias set reaches
        # compare_panel_content, and calling it here would do so while this
        # dict is built - i.e. on the drawing path, in the capsule, where that
        # module is deliberately not shipped. run() calls it only under
        # --check.
        "S9_D": (lambda: draw_D(fr), lambda: A._panel_cytokines(fr["cyto"]),
                 base.aliases),
    })


if __name__ == "__main__":
    sys.exit(main())
