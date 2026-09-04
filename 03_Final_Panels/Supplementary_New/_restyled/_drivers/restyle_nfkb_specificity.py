#!/usr/bin/env python3
"""
Restyle S9E and S9F - NF-kB specificity and the epithelial compartment
(Reviewer 1, point R1.8).

Analysis: 02_New_Analyses/07_R1.8_NFkB_Specificity/scripts/nfkb_specificity.py

NOTHING IS RECOMPUTED HERE. In particular the driver never reaches
`nfkb_specificity.py:100`, the `RECOMPUTE_GSEA.glob("*_hallmark.csv")` that
`load_gsea()` walks: that loop re-reads the module-12 Hallmark tables, re-ranks
all 50 sets per file, and would silently drop a cell type through the
`if term is None: continue` on line 105 if a table were renamed. A restyle must
not be able to change which cell types appear in S9E, so the driver reads
`outputs/nfkb_per_celltype.csv` - the table that same loop already wrote, and
which is exactly the `g` frame `main()` hands to `_panel_nfkb`.

`_panel_cytokines` is given `cyto` in the same way, from
`outputs/epithelial_cytokines.csv`. No h5ad is opened, so `sample_means()` and
its Mann-Whitney tests are not re-run either.

Both panel functions do their own selection inside the panel
(`method == PRIMARY and phase == "post"` for S9E, `phase == "Post"` for S9F);
those selections are copied verbatim below and `--check` proves they select the
same rows.

    conda run -n Liudeng_Python_310 python restyle_nfkb_specificity.py --check
    conda run -n Liudeng_Python_310 python restyle_nfkb_specificity.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _driver_base as base
import panel_style_cns as style

A = base.analysis("07_R1.8_NFkB_Specificity/scripts/nfkb_specificity.py")
OUT = base.outputs_of(A)
FIG = "S9_Mechanism_Specificity"

# Printed boxes, millimetres. Version A: S9E 100 x 52, S9F 65 x 46.
E_W, E_H = 100.0, 66.0
F_W, F_H = 78.0, 46.0


def frames():
    """The two tables the analysis already wrote. No computation."""
    return dict(
        g=pd.read_csv(base.require(
            OUT / "nfkb_per_celltype.csv", "S9E NF-kB per cell type")),
        cyto=pd.read_csv(base.require(
            OUT / "epithelial_cytokines.csv", "S9F epithelial cytokines")),
    )


# --------------------------------------------------------------------------
# S9E - post-treatment NES per cell type, with FDR and Hallmark rank annotated.
# Layout unchanged: one horizontal bar per cell type, sorted by NES.
# --------------------------------------------------------------------------
def draw_E(fr, save=True):
    style.apply()
    g = fr["g"]
    post = (g[(g["method"] == A.PRIMARY) & (g["phase"] == "post")]
            .sort_values("nes"))
    fig, ax = style.subplots_mm(E_W, E_H)
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
    style.margins_mm(fig, left=18, right=2, top=3, bottom=13)
    if save:
        base.save(fig, FIG, "S9_E", "S9_E_nfkb_per_celltype")
    return fig


# --------------------------------------------------------------------------
# S9F - epithelial versus myeloid cytokine expression after treatment.
# Layout unchanged: two grouped bars per gene on a symlog axis.
# --------------------------------------------------------------------------
def draw_F(fr, save=True):
    style.apply()
    post = fr["cyto"][fr["cyto"]["phase"] == "Post"]
    fig, ax = style.subplots_mm(F_W, F_H)
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
    ax.tick_params(axis="both", width=0.6, length=2)
    # The one place in the restyled supplementary set where 7 pt is not enough.
    # A symlog axis labels its decades as mathtext powers, and mathtext draws a
    # superscript at 70% of its base, so the exponent of 10^-2 printed at
    # 4.9 pt - the only type below the 5 pt floor on any restyled page. Setting
    # this axis's labels at the body size instead puts the exponent at 5.6 pt.
    # Local, deliberate, and recorded in RESTYLE_REPORT.md; the tick values and
    # the labels themselves are untouched.
    ax.tick_params(axis="y", labelsize=style.body_pt())
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="upper left",
              handlelength=1.0, handletextpad=0.4, labelspacing=0.25)
    style.margins_mm(fig, left=20, right=2, top=3, bottom=9)
    if save:
        base.save(fig, FIG, "S9_F", "S9_F_epithelial_vs_myeloid_cytokines")
    return fig


def main():
    fr = frames()
    return base.run({
        "S9_E": (lambda: draw_E(fr), lambda: A._panel_nfkb(fr["g"])),
        "S9_F": (lambda: draw_F(fr), lambda: A._panel_cytokines(fr["cyto"])),
    })


if __name__ == "__main__":
    sys.exit(main())
