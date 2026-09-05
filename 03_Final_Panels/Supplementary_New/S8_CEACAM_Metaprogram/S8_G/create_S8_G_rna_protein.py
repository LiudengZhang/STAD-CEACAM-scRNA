"""
Supplementary Figure S8G - transcript against protein in the same eight tumors.

The immunohistochemistry was performed on the same eight pre-treatment patients
as the single-cell cohort, so it is not an independent cohort. It answers a
different question instead: whether the single-cell measurement reflects protein
abundance in the same tissue. It does.

Input : 04_Revision_Analyses/10_R1.3_CrossCohort_Convergence/outputs/
Output: S8_G_rna_protein_concordance.{svg,pdf,png}
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import NEW_ANALYSES  # noqa: E402

OUT_DIR = Path(__file__).parent
CONV = NEW_ANALYSES / "10_R1.3_CrossCohort_Convergence" / "outputs"

SCALE = 4
CM = 1 / 2.54
DPI = 300
PANEL_W_CM = 6.5 * SCALE
PANEL_H_CM = 3.0 * SCALE

GROUP_COLOR = {"NR": "#B2182B", "R": "#2166AC"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def main():
    per = pd.read_csv(CONV / "rna_protein_per_patient.csv")
    conc = pd.read_csv(CONV / "rna_protein_concordance.csv").set_index("Marker")
    per["rna"] = (per["CEACAM5_rna"] + per["CEACAM6_rna"]) * 100
    per["protein"] = per["CEACAM5_protein"] + per["CEACAM6_protein"]

    fig, ax = plt.subplots(figsize=(PANEL_W_CM * CM, PANEL_H_CM * CM))
    for group, d in per.groupby("group_rna"):
        ax.scatter(d["rna"], d["protein"], s=34, c=GROUP_COLOR[group],
                   edgecolors="white", linewidths=0.4, label=group, zorder=3)
    for _, r in per.iterrows():
        ax.annotate(r["patient"], (r["rna"], r["protein"]),
                    textcoords="offset points", xytext=(4, 3),
                    fontsize=4.5 * SCALE, color="#555555")

    rho = conc.loc["Summed", "Spearman rho"]
    p = conc.loc["Summed", "P"]
    ax.text(0.03, 0.95, f"Spearman $\\rho$ = {rho:.2f}, $P$ = {p:.3f}",
            transform=ax.transAxes, fontsize=5 * SCALE, ha="left", va="top")

    # Both axes are sums of two per-marker percentages, so both can exceed 100:
    # a cell or a region positive for both markers counts in each term. The two
    # modalities are therefore on the same composite scale.
    ax.set_xlabel("Summed CEACAM5 and CEACAM6 positive\nfractions (%), single cell",
                  fontsize=5 * SCALE)
    ax.set_ylabel("Summed staining (%),\nimmunohistochemistry", fontsize=5 * SCALE)
    ax.tick_params(labelsize=5 * SCALE, width=0.5, length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.5)
    ax.legend(loc="lower right", frameon=False, fontsize=5 * SCALE,
              handletextpad=0.3, borderpad=0.2)

    fig.subplots_adjust(left=0.20, right=0.98, top=0.96, bottom=0.30)
    stem = OUT_DIR / "S8_G_rna_protein_concordance"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
