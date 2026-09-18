#!/usr/bin/env python3
"""Build the complete reviewer-only two-sided sensitivity figure.

The earlier embedded artwork was a stale 19-row export and shortened three
labels with ellipses.  This version reads the current sweep directly, includes
all 20 comparisons in Table S6, and uses the panel numbering in the submitted
revision.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SWEEP = (ROOT / "04_Revision_Analyses" / "02_R1.3_TwoSided_Stats_Sweep" /
         "outputs" / "twosided_sweep.csv")
OUT = (ROOT / "04_Manuscript_R1" / "05_Response_to_Reviewers" / "_embed" /
       "Figure_R1_complete_twosided_sensitivity.png")

COLORS = {
    "CEACAM": "#2166AC",
    "Metaprogram": "#4393C3",
    "PD-L1": "#B2182B",
    "Spatial": "#7B3294",
    "Cytokine": "#D6604D",
    "Regulon": "#1B7837",
}


def display_label(row):
    analysis = row["analysis"]
    panel = row["panel"]
    if panel == "Fig 2 (N1)": return "Fig. 2K  CEACAM6, pretreatment epithelium"
    if panel == "Fig 2 (N2)": return "Fig. 2K  CEACAM5, pretreatment epithelium"
    if panel == "Fig 2 (O1)": return "Fig. 2L  CEACAM6, external-cohort epithelium"
    if panel == "Fig 2 (O2)": return "Fig. 2L  CEACAM5, external-cohort epithelium"
    if panel == "Fig 2 (Q)":
        marker = analysis.split("IHC staining, ", 1)[1]
        return f"Fig. 2N  IHC, {marker}"
    if panel == "Fig S2 (D)": return "Fig. S3D  Tumor-adjusted CEACAM5/6+ proportion"
    if panel == "Fig 5 (IL6-CD4)": return "Fig. 5L  IL-6/JAK/STAT3, post-treatment CD4+ T"
    if panel == "Fig 2 (K)":
        return ("Fig. 2H  S-MP4, pre-R vs all other groups" if "S-MP4" in analysis
                else "Fig. 2I  S-MP5, pre-R vs all other groups")
    if panel == "Fig 5 (J)": return "Fig. 5J  CD274, monocytes/macrophages"
    if panel == "Fig 5 (K)": return "Fig. 5J  CD274, epithelial cells"
    if panel == "Fig 5 (L)": return "Fig. 5J  CD274, fibroblasts"
    if panel == "Fig 5 (M)": return "Fig. 5J  CD274, dendritic cells"
    if panel == "Fig 5 (regulon-BACH1)": return "Fig. 5D  BACH1 regulon, post-R vs all others"
    if panel == "Fig 5 (regulon-NFKB1)": return "Fig. 5E  NFKB1 regulon, post-R vs all others"
    if panel == "Fig 3 (H)": return "Fig. 3H  Neighborhood epithelial density"
    if panel == "Fig 3 (I)": return "Fig. 3I  Distance to stroma"
    if panel == "Fig 3 (J)": return "Fig. 3J  Distance to immune-rich regions"
    raise ValueError(f"Unmapped sensitivity row: {panel} / {analysis}")


def main():
    df = pd.read_csv(SWEEP).copy()
    flip = df["family"].isin(("Metaprogram", "Regulon"))
    for col in ("hedges_g", "g_ci95_lo", "g_ci95_hi"):
        df.loc[flip, col] = -df.loc[flip, col]
    lo = df.loc[flip, "g_ci95_lo"].copy()
    df.loc[flip, "g_ci95_lo"] = df.loc[flip, "g_ci95_hi"].values
    df.loc[flip, "g_ci95_hi"] = lo.values
    order = ["CEACAM", "Metaprogram", "PD-L1", "Cytokine", "Regulon", "Spatial"]
    df["_family"] = pd.Categorical(df["family"], order, ordered=True)
    df = df.sort_values(["_family", "panel"], kind="stable").reset_index(drop=True)
    labels = [display_label(r) for _, r in df.iterrows()]
    if len(labels) != 20 or any("..." in x for x in labels):
        raise RuntimeError("Sensitivity figure must contain 20 complete labels")

    y = np.arange(len(df))[::-1]
    fig = plt.figure(figsize=(5.6, 7.55), dpi=300)
    gs = fig.add_gridspec(1, 3, width_ratios=[2.85, 1.75, 0.95], wspace=0.03,
                          left=0.015, right=0.985, top=0.935, bottom=0.11)
    axl = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1], sharey=axl)
    axp = fig.add_subplot(gs[0, 2], sharey=axl)

    axl.set_xlim(0, 1); axl.set_ylim(-0.8, len(df)-0.2); axl.axis("off")
    for yy, label in zip(y, labels):
        axl.text(0.99, yy, label, ha="right", va="center", fontsize=6.2)

    xmin, xmax = -1.0, 5.0
    ax.axvline(0, color="#999999", lw=0.55, ls="--", zorder=1)
    for yy, (_, row) in zip(y, df.iterrows()):
        color = COLORS[row["family"]]
        low = max(float(row["g_ci95_lo"]), xmin)
        high = min(float(row["g_ci95_hi"]), xmax)
        ax.plot([low, high], [yy, yy], color=color, lw=0.8, zorder=2)
        if row["g_ci95_lo"] < xmin:
            ax.plot(xmin, yy, marker="<", color=color, ms=2.8, clip_on=False)
        if row["g_ci95_hi"] > xmax:
            ax.plot(xmax, yy, marker=">", color=color, ms=2.8, clip_on=False)
        ax.scatter(np.clip(row["hedges_g"], xmin, xmax), yy, s=9,
                   color=color, edgecolor="white", linewidth=0.3, zorder=3)
    ax.set_xlim(xmin-0.08, xmax+0.08); ax.set_ylim(-0.8, len(df)-0.2)
    ax.set_yticks([]); ax.tick_params(axis="x", labelsize=6.2, length=2, width=0.5)
    ax.set_xlabel("Hedges' g (95% bootstrap CI)", fontsize=7)
    for side in ("top", "right", "left"): ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)

    axp.set_xlim(0, 1); axp.set_ylim(ax.get_ylim()); axp.axis("off")
    axp.text(0.02, len(df)-0.05, "P 1-tail", ha="left", va="bottom",
             fontsize=6.2, fontstyle="italic")
    axp.text(0.55, len(df)-0.05, "P 2-tail", ha="left", va="bottom",
             fontsize=6.2, fontstyle="italic")
    for yy, (_, row) in zip(y, df.iterrows()):
        axp.text(0.02, yy, f"{row['p_one_tailed']:.3f}", ha="left", va="center",
                 fontsize=6.2, color="#666666")
        sig = row["p_two_tailed"] < 0.05
        axp.text(0.55, yy, f"{row['p_two_tailed']:.3f}", ha="left", va="center",
                 fontsize=6.2, color="black" if sig else "#666666",
                 fontweight="bold" if sig else "normal")

    n_ci = int(df["ci_excludes_zero"].sum())
    fig.suptitle(f"All {len(df)} directional comparisons repeated two-sided; "
                 f"{n_ci}/{len(df)} retain a 95% CI excluding zero",
                 fontsize=8.2, y=0.985)
    handles = [plt.Line2D([], [], color=c, marker="o", lw=0.8, ms=2.8, label=k)
               for k, c in COLORS.items()]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=6.2, handlelength=1.7, columnspacing=1.2,
               bbox_to_anchor=(0.55, 0.012))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, facecolor="white")
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
