"""
Supplementary Figure S8E, RESTYLED (Version B) - convergence of the
pre-treatment CEACAM5/6 evidence.

Copy of Supplementary_New/S8_CEACAM_Metaprogram/S8_E/create_S8_E_convergence.py
with the type taken from cnsplots via 00_Config/panel_style_cns.py, the canvas
set to the millimetre box the panel prints in, and the point-specified marks
rescaled so their size relative to the type is unchanged. Every value read,
grouped, sorted and plotted is the same code.

Version A drew 68.4 x 13.6 cm at SCALE 4 into a 171 x 34 mm box at 0.2500, so
its row labels printed at 5.0 pt. Seven rows at 7 pt need a taller panel and a
wider label column.

Reviewer 1 asked what a four-versus-four comparison can support. Read one cohort
at a time each CEACAM measurement sits near P = 0.06; read together, the two
cohorts that share no patients give a combined P of 0.007 for CEACAM6. This
panel shows every measurement on one standardised scale, grouped by whether it
is independent of the discovery cohort, with the combined result stated.

Input : 04_Revision_Analyses/10_R1.3_CrossCohort_Convergence/outputs/
Output: S8_E_convergence_forest.{svg,pdf,png}
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[5] / "00_Config"))
from paths import NEW_ANALYSES  # noqa: E402
import panel_style_cns as style  # noqa: E402

OUT_DIR = Path(__file__).parent
CONV = NEW_ANALYSES / "10_R1.3_CrossCohort_Convergence" / "outputs"

PANEL_W_MM = 171.0
PANEL_H_MM = 52.0

A_FIT = 0.2500
A_TYPE = 5.0 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2

GROUP_COLOR = {
    "discovery cohort": "#2166AC",
    "independent of the discovery cohort": "#1B7837",
    "same patients as the discovery cohort": "#B2182B",
}
GROUP_ORDER = list(GROUP_COLOR)

LABEL = {
    "CEACAM6 expression, pre-treatment epithelium": "CEACAM6, single-cell (4 vs 4)",
    "CEACAM5 expression, pre-treatment epithelium": "CEACAM5, single-cell (4 vs 4)",
    "CEACAM6 expression, PRJEB25780 deconvolved epithelium":
        "CEACAM6, PRJEB25780 (33 vs 12)",
    "CEACAM5 expression, PRJEB25780 deconvolved epithelium":
        "CEACAM5, PRJEB25780 (33 vs 12)",
    "IHC staining, CEACAM5 + CEACAM6 summed": "IHC, summed (4 vs 4)",
    "IHC staining, CEACAM5 only": "IHC, CEACAM5 (4 vs 4)",
    "IHC staining, CEACAM6 only": "IHC, CEACAM6 (4 vs 4)",
}


def main():
    style.apply()
    df = pd.read_csv(CONV / "forest_effect_sizes.csv")
    comb = pd.read_csv(CONV / "crosscohort_combination.csv")
    df["group"] = df["Independent of the scRNA cohort"].map({
        "discovery cohort": "discovery cohort",
        "yes": "independent of the discovery cohort",
        "no - same patients": "same patients as the discovery cohort"})
    df["order"] = df["group"].map({g: i for i, g in enumerate(GROUP_ORDER)})
    df = df.sort_values(["order", "Measurement"]).reset_index(drop=True)
    y = np.arange(len(df))[::-1]

    fig, (ax, axp) = style.subplots_mm(
        PANEL_W_MM, PANEL_H_MM, 1, 2,
        gridspec_kw={"width_ratios": [2.6, 1.0], "wspace": 0.03})

    X_MIN, X_MAX = -0.6, 4.0
    ax.axvline(0, color="#999999", linewidth=0.8 * MARK, linestyle="--",
                zorder=1)
    for i, (_, r) in enumerate(df.iterrows()):
        c = GROUP_COLOR[r["group"]]
        lo, hi = max(r["g 95% CI low"], X_MIN), min(r["g 95% CI high"], X_MAX)
        ax.plot([lo, hi], [y[i]] * 2, color=c, linewidth=1.2 * MARK,
                solid_capstyle="round", zorder=2)
        if r["g 95% CI high"] > X_MAX:
            ax.plot(X_MAX, y[i], marker=">", color=c, markersize=3 * MARK,
                    clip_on=False, zorder=3)
        ax.scatter(min(r["Hedges g"], X_MAX), y[i], s=26 * AREA, c=c,
                   edgecolors="white", linewidths=0.4 * MARK, zorder=4)
    ax.set_xlim(X_MIN - 0.1, X_MAX + 0.1)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[m] for m in df["Measurement"]])
    ax.set_xlabel("Higher in non-responders  (Hedges' $g$, 95% bootstrap CI)")
    ax.tick_params(axis="x", width=0.6, length=2)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.7, len(df) - 0.1)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    axp.set_ylim(ax.get_ylim())
    axp.axis("off")
    axp.text(0.0, len(df) - 0.55, "$P$, two-sided", fontsize=style.tick_pt(),
             ha="left", va="center", style="italic")
    for i, (_, r) in enumerate(df.iterrows()):
        axp.text(0.0, y[i], f"{r['P, two-sided']:.3f}",
                 fontsize=style.tick_pt(), ha="left", va="center",
                 color="#333333")
    stouffer = comb[comb.Method == "Stouffer, unweighted"].set_index("Gene")
    axp.text(0.30, len(df) / 2 - 0.5,
             "Combined across the two cohorts\n"
             "that share no patients (Stouffer):\n"
             f"CEACAM6  $P$ = {stouffer.loc['CEACAM6', 'P, combined two-sided']:.3f}\n"
             f"CEACAM5  $P$ = {stouffer.loc['CEACAM5', 'P, combined two-sided']:.3f}",
             fontsize=style.tick_pt(), ha="left", va="center", linespacing=1.5,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#F2F2F2",
                       edgecolor="#CCCCCC", linewidth=0.5))
    axp.set_xlim(-0.05, 1.05)

    # cnsplots halves legend keys (legend.markerscale 0.5); size each key so it
    # prints at the same diameter as the mark it stands for.
    key = 1.0 / plt.rcParams["legend.markerscale"]
    handles = [plt.Line2D([], [], color=c, linewidth=1.2 * MARK, marker="o",
                          markersize=(26 * AREA) ** 0.5 * key, label=g)
               for g, c in GROUP_COLOR.items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.30),
              ncol=3, handletextpad=0.5, columnspacing=1.4)

    style.margins_mm(fig, left=42, right=2, top=3, bottom=16)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUT_DIR / "S8_E_convergence_forest")
    print(f"Saved S8_E_convergence_forest at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
