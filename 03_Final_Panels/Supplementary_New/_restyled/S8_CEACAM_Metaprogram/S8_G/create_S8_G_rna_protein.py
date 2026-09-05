"""
Supplementary Figure S8G, RESTYLED (Version B) - transcript against protein in
the same eight tumors.

Copy of Supplementary_New/S8_CEACAM_Metaprogram/S8_G/create_S8_G_rna_protein.py
with the type taken from cnsplots via 00_Config/panel_style_cns.py, the canvas
set to the millimetre box the panel prints in, and the point-specified marks
rescaled so their size relative to the type is unchanged. Every value read,
computed and plotted is the same code.

Version A drew 26 x 12 cm at SCALE 4 and the assembler fitted it into a
65 x 26 mm box at 0.2167, so its smallest type - the patient labels at
4.5 * SCALE - printed at 3.90 pt.

Input : 04_Revision_Analyses/10_R1.3_CrossCohort_Convergence/outputs/
Output: S8_G_rna_protein_concordance.{svg,pdf,png}
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[5] / "00_Config"))
from paths import NEW_ANALYSES  # noqa: E402
import panel_style_cns as style  # noqa: E402

OUT_DIR = Path(__file__).parent
CONV = NEW_ANALYSES / "10_R1.3_CrossCohort_Convergence" / "outputs"

PANEL_W_MM = 65.0
PANEL_H_MM = 46.0

A_FIT = 0.2167
A_TYPE = 4.5 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2

GROUP_COLOR = {"NR": "#B2182B", "R": "#2166AC"}


def main():
    style.apply()
    per = pd.read_csv(CONV / "rna_protein_per_patient.csv")
    conc = pd.read_csv(CONV / "rna_protein_concordance.csv").set_index("Marker")
    per["rna"] = (per["CEACAM5_rna"] + per["CEACAM6_rna"]) * 100
    per["protein"] = per["CEACAM5_protein"] + per["CEACAM6_protein"]

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    for group, d in per.groupby("group_rna"):
        ax.scatter(d["rna"], d["protein"], s=34 * AREA, c=GROUP_COLOR[group],
                   edgecolors="white", linewidths=0.4 * MARK, label=group,
                   zorder=3)
    for _, r in per.iterrows():
        ax.annotate(r["patient"], (r["rna"], r["protein"]),
                    textcoords="offset points", xytext=(4 * MARK, 3 * MARK),
                    fontsize=style.tick_pt(), color="#555555")

    rho = conc.loc["Summed", "Spearman rho"]
    p = conc.loc["Summed", "P"]
    ax.text(0.03, 0.95, f"Spearman $\\rho$ = {rho:.2f}, $P$ = {p:.3f}",
            transform=ax.transAxes, fontsize=style.tick_pt(), ha="left",
            va="top")

    # Both axes are sums of two per-marker percentages, so both can exceed 100:
    # a cell or a region positive for both markers counts in each term. The two
    # modalities are therefore on the same composite scale.
    ax.set_xlabel("Summed CEACAM5 and CEACAM6 positive\nfractions (%), single cell")
    ax.set_ylabel("Summed staining (%),\nimmunohistochemistry")
    ax.tick_params(width=0.6, length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(loc="lower right", handletextpad=0.3, borderpad=0.2)

    style.margins_mm(fig, left=16, right=2, top=2, bottom=14)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUT_DIR / "S8_G_rna_protein_concordance")
    print(f"Saved S8_G_rna_protein_concordance at "
          f"{PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
