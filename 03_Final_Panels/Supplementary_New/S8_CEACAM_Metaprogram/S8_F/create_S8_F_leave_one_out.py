"""
Supplementary Figure S8F - leave-one-patient-out stability of the four-versus-
four pre-treatment comparison.

The response label of one patient moves the result substantially (Table S8).
Which patients were sampled does not: dropping any single patient leaves both
the direction and the magnitude of the difference intact. The two questions are
different and this panel separates them.

Input : 04_Revision_Analyses/10_R1.3_CrossCohort_Convergence/outputs/loo_stability.csv
Output: S8_F_leave_one_out.{svg,pdf,png}
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import NEW_ANALYSES  # noqa: E402

OUT_DIR = Path(__file__).parent
LOO = (NEW_ANALYSES / "10_R1.3_CrossCohort_Convergence" / "outputs"
       / "loo_stability.csv")

SCALE = 4
CM = 1 / 2.54
DPI = 300
PANEL_W_CM = 10.0 * SCALE
PANEL_H_CM = 3.0 * SCALE

GENE_COLOR = {"CEACAM6": "#2166AC", "CEACAM5": "#4393C3"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def main():
    df = pd.read_csv(LOO)
    genes = ["CEACAM6", "CEACAM5"]

    fig, ax = plt.subplots(figsize=(PANEL_W_CM * CM, PANEL_H_CM * CM))
    ax.axvline(0, color="#999999", linewidth=0.8, linestyle="--", zorder=1)

    for row, gene in enumerate(genes):
        y = len(genes) - 1 - row
        d = df[df.Gene == gene]
        drops = d[d.Dropped != "none (as published)"]
        published = d[d.Dropped == "none (as published)"].iloc[0]
        c = GENE_COLOR[gene]
        ax.scatter(drops["Hedges g"], np.full(len(drops), y), s=22, c=c,
                   alpha=0.75, edgecolors="white", linewidths=0.3, zorder=3)
        ax.scatter(published["Hedges g"], y, s=70, marker="D", c="white",
                   edgecolors=c, linewidths=1.0, zorder=4)
        ax.text(3.9, y + 0.30,
                f"{int((drops['Hedges g'] > 0).sum())} of {len(drops)} refits keep "
                f"the direction; $g$ {drops['Hedges g'].min():.2f}"
                f"–{drops['Hedges g'].max():.2f}",
                fontsize=5 * SCALE, ha="right", va="center", color="#555555")

    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes[::-1], fontsize=5.5 * SCALE)
    ax.set_ylim(-0.6, len(genes) - 0.4)
    ax.set_xlim(-0.4, 4.0)
    ax.set_xlabel("Effect size after dropping one patient (Hedges' $g$)",
                  fontsize=5.5 * SCALE)
    ax.tick_params(axis="x", labelsize=5 * SCALE, width=0.5, length=2)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)

    handles = [
        plt.Line2D([], [], marker="D", color="white", markeredgecolor="#444444",
                   markersize=5, linestyle="none", label="all eight patients"),
        plt.Line2D([], [], marker="o", color="#444444", markersize=4,
                   linestyle="none", label="one patient dropped"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.34),
              ncol=2, frameon=False, fontsize=5 * SCALE, handletextpad=0.4,
              columnspacing=1.6)

    fig.subplots_adjust(left=0.12, right=0.99, top=0.95, bottom=0.34)
    stem = OUT_DIR / "S8_F_leave_one_out"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
