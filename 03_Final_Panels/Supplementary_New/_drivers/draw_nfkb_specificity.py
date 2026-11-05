#!/usr/bin/env python3
"""Draw Supplementary Figure S9 C (post-treatment NF-kB NES per cell type)
and S9 D (epithelial vs myeloid cytokine expression).

Reads the tables written by
04_Revision_Analyses/07_R1.8_NFkB_Specificity/scripts/nfkb_specificity.py:
C from nfkb_per_celltype_sound13.csv (`A.load_adopted()`), D from
outputs/epithelial_cytokines.csv. No h5ad is opened.


    python draw_nfkb_specificity.py [--check]"""

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

# Panel sizes, mm.
C_W, C_H = 66.0, 62.0
D_W, D_H = 48.0, 46.0

def frames():
    """The two tables the panels are drawn from."""
    base.require(A.ADOPTED_GSEA, "S9C NF-kB per cell type (adopted table)")
    return dict(
        g=A.load_adopted(),
        cyto=pd.read_csv(base.require(
            OUT / "epithelial_cytokines.csv", "S9D epithelial cytokines")),
    )

def draw_C(fr, save=True):
    """S9 C - post-treatment NES per cell type, one horizontal bar each sorted
    by NES, with FDR and Hallmark rank annotated."""
    base.apply_style()
    g = fr["g"]
    post = (g[(g["method"] == A.PRIMARY) & (g["phase"] == "post")]
            .sort_values("nes"))
    fig, ax = style.subplots_mm(C_W, C_H)
    y = np.arange(len(post))
    colors = [A.COLOR_NR if q < 0.25 else "#cccccc" for q in post["fdr_q"]]
    ax.barh(y, post["nes"], color=colors, edgecolor="#444444",
            linewidth=style.EDGE_PT, height=0.65)
    # Annotate right of zero for negative bars, clear of the y labels.
    for i, (_, r) in enumerate(post.iterrows()):
        ax.text(max(r["nes"], 0.0) + 0.05, i,
                f"q={r['fdr_q']:.2f}  #{r['rank']}/{r['n_sets']}",
                va="center", ha="left",
                fontsize=style.tick_pt(), color="#333333")
    ax.axvline(0, color="#666666", linewidth=style.RULE_PT)
    ax.set_yticks(y)
    ax.set_yticklabels([A.LABELS[c] for c in post["cell_type"]])
    ax.set_xlabel("NES, TNF$\\alpha$ signalling via NF-$\\kappa$B\n"
                  "(post-treatment, non-responders vs responders)")
    ax.set_xlim(min(0, post["nes"].min()) - 0.3, post["nes"].max() + 1.9)
    ax.tick_params(axis="x", width=style.RULE_PT, length=2)
    ax.tick_params(axis="y", length=0)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="#444444",
                             linewidth=style.EDGE_PT, label=l)
               for c, l in ((A.COLOR_NR, "FDR q < 0.25"),
                            ("#cccccc", "FDR q ≥ 0.25"))]
    # Key above the plotting box, clear of the bar annotations.
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.0),
              borderaxespad=0.0, frameon=False,
              handlelength=1.0, handletextpad=0.4, labelspacing=0.25)
    base.fit(fig)
    if save:
        base.save(fig, FIG, "S9_C", "S9_C_nfkb_per_celltype")
    return fig

def draw_D(fr, save=True):
    """S9 D - epithelial versus myeloid cytokine expression after treatment,
    two grouped bars per gene on a symlog axis."""
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
               edgecolor="#333333", linewidth=style.EDGE_PT, label=comp)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, style="italic")
    ax.set_ylabel("Mean expression,\npost-treatment non-responders")
    ax.set_yscale("symlog", linthresh=0.01)
    # Decades written out (0.01, 0.1, 1) instead of powers of ten.
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, _: "0" if v == 0 else f"{v:g}"))
    base.finish(ax)
    # Key above the plotting box.
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0),
              borderaxespad=0.0, ncol=1,
              handlelength=1.0, handletextpad=0.4, labelspacing=0.25)
    base.fit(fig)
    if save:
        base.save(fig, FIG, "S9_D", "S9_D_epithelial_vs_myeloid_cytokines")
    return fig

def main():
    fr = frames()
    return base.run({
        "S9_C": (lambda: draw_C(fr), lambda: A._panel_nfkb(fr["g"])),
        "S9_D": (lambda: draw_D(fr), lambda: A._panel_cytokines(fr["cyto"]),
                 base.aliases),
    })

if __name__ == "__main__":
    sys.exit(main())
