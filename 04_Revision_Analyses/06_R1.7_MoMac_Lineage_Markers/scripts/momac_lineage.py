"""
WP7 / Reviewer 1 point R1.7.

  "The cluster is described as 'IL-1b+ macrophages,' but its markers (IL1B, TNF,
   CXCL8, CCL3, CXCL2, and related genes) may also characterize inflammatory
   monocytes or recently recruited monocyte-derived cells. It would therefore be
   appropriate to include evidence supporting the annotation, such as macrophage
   versus monocyte lineage markers."

Approach: score every monocyte/macrophage state against a canonical monocyte
panel and a canonical macrophage panel, then compare the two scores within each
state. A cluster that is genuinely macrophage should sit clearly on the
macrophage side of the two reference states already present in the object
(C1_Mono_Classic_CD14 as the monocyte pole, C0_Mac_Classic_TREM2 as the
macrophage pole).

Note on the published figure: the marker dotplot in Figure 4E lists "CD16",
which is not a valid HGNC symbol. The gene is FCGR3A and is used here.

Input : Round_5/01_Raw_Inputs/01_H5AD/MoMac.h5ad
Output: momac_lineage_scores.csv, momac_lineage_tests.csv,
        momac_lineage_report.txt, panels S9_C and S9_D
"""

from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD, REVISED_PANELS  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S9 = REVISED_PANELS / "Supplementary_New" / "S9_Mechanism_Specificity"

SCALE, CM, DPI = 4, 1 / 2.54, 300

# Canonical lineage panels. The monocyte set is the classical/non-classical
# circulating programme; the macrophage set is the tissue-resident/differentiated
# programme including complement and lipid handling.
MONOCYTE = ["FCN1", "VCAN", "S100A8", "S100A9", "S100A12", "SELL", "CD14",
            "FCGR3A", "LYZ", "CSF3R"]
MACROPHAGE = ["CD68", "CD163", "MRC1", "C1QA", "C1QB", "C1QC", "APOE", "APOC1",
              "TREM2", "MERTK", "MSR1", "SEPP1", "SELENOP"]

MONO_POLE = "C1_Mono_Classic_CD14"
MAC_POLE = "C0_Mac_Classic_TREM2"
TARGET = "C3_Mac_Inflam_IL1B"
# The cluster is renamed in the revision from "macrophage" to "MoMac", because
# this analysis places it on the monocyte-macrophage continuum rather than at
# the macrophage pole. The stored label is left untouched so every downstream
# table still keys on it; only what the figure prints changes.
DISPLAY = {TARGET: "C3_MoMac_Inflam_IL1B"}


def shown(state):
    return DISPLAY.get(str(state), str(state))

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def main():
    ad = sc.read_h5ad(MOMAC_H5AD)

    present = set(ad.var_names)
    mono = [g for g in MONOCYTE if g in present]
    mac = [g for g in MACROPHAGE if g in present]
    missing = [g for g in MONOCYTE + MACROPHAGE if g not in present]

    sc.tl.score_genes(ad, mono, score_name="monocyte_score", random_state=0)
    sc.tl.score_genes(ad, mac, score_name="macrophage_score", random_state=0)
    ad.obs["lineage_index"] = ad.obs["macrophage_score"] - ad.obs["monocyte_score"]

    obs = ad.obs[["minor_cell_state", "monocyte_score", "macrophage_score",
                  "lineage_index"]].copy()
    obs["minor_cell_state"] = obs["minor_cell_state"].astype(str)

    summary = (obs.groupby("minor_cell_state")
               .agg(n_cells=("lineage_index", "size"),
                    monocyte_score=("monocyte_score", "mean"),
                    macrophage_score=("macrophage_score", "mean"),
                    lineage_index=("lineage_index", "mean"))
               .sort_values("lineage_index"))
    summary.to_csv(OUT / "momac_lineage_scores.csv")

    # Is the IL-1B cluster closer to the monocyte pole or the macrophage pole?
    tests = []
    tgt = obs.loc[obs["minor_cell_state"] == TARGET, "lineage_index"].values
    for pole, name in ((MONO_POLE, "monocyte pole"), (MAC_POLE, "macrophage pole")):
        ref = obs.loc[obs["minor_cell_state"] == pole, "lineage_index"].values
        u, p = stats.mannwhitneyu(tgt, ref, alternative="two-sided")
        tests.append(dict(
            comparison=f"{TARGET} vs {pole} ({name})",
            n_target=len(tgt), n_reference=len(ref),
            mean_target=float(tgt.mean()), mean_reference=float(ref.mean()),
            rank_biserial_r=float(2.0 * u / (len(tgt) * len(ref)) - 1.0),
            p_two_tailed=float(p)))
    tests = pd.DataFrame(tests)
    tests.to_csv(OUT / "momac_lineage_tests.csv", index=False)

    # Where does the IL-1B cluster sit on the pole-to-pole axis? 0 = monocyte
    # pole, 1 = macrophage pole, using the state means.
    m_lo = summary.loc[MONO_POLE, "lineage_index"]
    m_hi = summary.loc[MAC_POLE, "lineage_index"]
    pos = (summary.loc[TARGET, "lineage_index"] - m_lo) / (m_hi - m_lo)

    L = ["MONOCYTE VERSUS MACROPHAGE LINEAGE - Reviewer 1 point R1.7", "=" * 92, ""]
    L.append(f"Monocyte panel  ({len(mono)} genes): {', '.join(mono)}")
    L.append(f"Macrophage panel ({len(mac)} genes): {', '.join(mac)}")
    if missing:
        L.append(f"Not present in the object and therefore omitted: {', '.join(missing)}")
    L.append("")
    L.append("Mean scores per monocyte/macrophage state "
             "(lineage index = macrophage - monocyte):")
    L.append("-" * 92)
    L.append(f"  {'state':<34}{'n cells':>9}{'mono':>9}{'mac':>9}{'index':>9}")
    for state, r in summary.iterrows():
        mark = "   <-- IL-1B cluster" if state == TARGET else ""
        L.append(f"  {state:<34}{int(r['n_cells']):>9,}{r['monocyte_score']:>9.3f}"
                 f"{r['macrophage_score']:>9.3f}{r['lineage_index']:>9.3f}{mark}")
    L.append("")
    L.append("Position of the IL-1B cluster on the monocyte-to-macrophage axis")
    L.append("-" * 92)
    L.append(f"  0.00 = {MONO_POLE} (monocyte pole)")
    L.append(f"  1.00 = {MAC_POLE} (macrophage pole)")
    L.append(f"  {TARGET} sits at {pos:.2f}")
    L.append("")
    for _, t in tests.iterrows():
        L.append(f"  {t['comparison']}")
        L.append(f"     mean index {t['mean_target']:+.3f} vs {t['mean_reference']:+.3f}   "
                 f"r = {t['rank_biserial_r']:+.3f}   P = {t['p_two_tailed']:.3g}")
    L.append("")
    L.append("INTERPRETATION")
    L.append("-" * 92)
    if pos < 0.35:
        L.append("  The cluster sits nearer the monocyte pole than the macrophage pole.")
        L.append("  The annotation should be changed to 'IL-1b+ inflammatory monocyte-")
        L.append("  macrophage' and the manuscript should not claim a pure macrophage")
        L.append("  identity.")
    elif pos > 0.65:
        L.append("  The cluster sits nearer the macrophage pole, supporting the")
        L.append("  published annotation.")
    else:
        L.append("  The cluster sits between the two poles, consistent with a recently")
        L.append("  recruited monocyte-derived macrophage. The reviewer's alternative")
        L.append("  is not excluded, and the annotation is changed to 'IL-1b+")
        L.append("  inflammatory monocyte-macrophage' to reflect that.")

    report = "\n".join(L)
    (OUT / "momac_lineage_report.txt").write_text(report, encoding="utf-8")
    print(report)

    _panel_scatter(summary, obs)
    _panel_dotplot(ad, mono, mac)


def _panel_scatter(summary, obs):
    d = (S9 / "S9_C"); d.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.0 * SCALE * CM, 5.0 * SCALE * CM))
    for state, r in summary.iterrows():
        highlight = state == TARGET
        ax.scatter(r["monocyte_score"], r["macrophage_score"],
                   s=(70 if highlight else 45) * SCALE,
                   c="#B2182B" if highlight else "#4d4d4d",
                   edgecolors="white", linewidths=0.5 * SCALE, zorder=3)
        # C4 and C1 sit at almost the same height, and a label to the right of
        # C4 would end under C1's point and read as C1's. C4 is labelled on its
        # left instead, into empty space below the diagonal.
        label = shown(state).replace("_", " ")
        left = label.startswith("C4 ")
        ax.annotate(label, (r["monocyte_score"], r["macrophage_score"]),
                    textcoords="offset points",
                    xytext=((-6 if left else 6) * SCALE, 3 * SCALE),
                    ha="right" if left else "left",
                    fontsize=4.5 * SCALE,
                    color="#B2182B" if highlight else "#333333")
    lim = [min(summary["monocyte_score"].min(), summary["macrophage_score"].min()),
           max(summary["monocyte_score"].max(), summary["macrophage_score"].max())]
    pad = 0.12 * (lim[1] - lim[0])
    ax.plot([lim[0] - pad, lim[1] + pad], [lim[0] - pad, lim[1] + pad],
            color="#bbbbbb", linestyle="--", linewidth=0.8, zorder=1)
    # Labels sit to the right of their point and the longest is ~30 characters,
    # so the x axis needs room the data alone does not ask for.
    ax.set_xlim(lim[0] - pad, lim[1] + 5.5 * pad)
    ax.set_ylim(lim[0] - pad, lim[1] + pad)
    ax.set_xlabel("Monocyte signature score", fontsize=6 * SCALE)
    ax.set_ylabel("Macrophage signature score", fontsize=6 * SCALE)
    ax.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.subplots_adjust(left=0.17, right=0.97, top=0.96, bottom=0.13)
    stem = d / "S9_C_momac_lineage_scores"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


def _panel_dotplot(ad, mono, mac):
    d = (S9 / "S9_D"); d.mkdir(parents=True, exist_ok=True)
    genes = mono + mac
    ad = ad.copy()
    ad.obs["cell state"] = [shown(s) for s in ad.obs["minor_cell_state"].astype(str)]
    order = sorted(set(ad.obs["cell state"]))
    fig = sc.pl.dotplot(
        ad, genes, groupby="cell state", categories_order=order,
        var_group_positions=[(0, len(mono) - 1), (len(mono), len(genes) - 1)],
        var_group_labels=["Monocyte", "Macrophage"],
        standard_scale="var", show=False, return_fig=True,
        figsize=(len(genes) * 0.34, len(order) * 0.32))
    stem = d / "S9_D_momac_lineage_dotplot"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI)
    plt.close("all")
    print(f"Saved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
