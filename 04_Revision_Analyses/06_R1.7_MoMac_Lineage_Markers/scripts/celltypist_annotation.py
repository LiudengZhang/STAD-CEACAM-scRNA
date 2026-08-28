"""
WP7b / Reviewer 1 point R1.7 - independent annotation of the IL-1b+ cluster.

The marker-score analysis places the cluster between the monocyte and
macrophage poles. That is our own scoring of our own data, so it is worth asking
what an external, published reference calls these cells, with no input from us.

CellTypist Immune_All_Low is trained on annotated immune atlases and assigns a
label per cell independently of this study's clustering. It is already used in
the manuscript to annotate the two external validation cohorts, so it introduces
no new method.

The question is narrow: of the cells we call C3_Mac_Inflam_IL1B, what proportion
does an external reference label as macrophage, and what proportion as monocyte,
relative to the two reference states in the same object?

Input : Round_5/01_Raw_Inputs/01_H5AD/MoMac.h5ad
Output: celltypist_labels.csv, celltypist_report.txt, panel S11_D
"""

from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD, REVISED_PANELS  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S11 = REVISED_PANELS / "Supplementary_New" / "S11_Affirmative_Analyses"

SCALE, CM, DPI = 4, 1 / 2.54, 300
MODEL = "Immune_All_Low.pkl"
TARGET = "C3_Mac_Inflam_IL1B"
# See momac_lineage.py: the printed label is updated, the stored one is not.
DISPLAY = {TARGET: "C3_MoMac_Inflam_IL1B"}
MONO_POLE = "C1_Mono_Classic_CD14"
MAC_POLE = "C0_Mac_Classic_TREM2"


def main():
    import celltypist
    from celltypist import models

    try:
        models.download_models(model=[MODEL], force_update=False)
    except Exception as e:
        raise SystemExit(f"Could not obtain the CellTypist model {MODEL}: {e}")

    ad = sc.read_h5ad(MOMAC_H5AD)
    # CellTypist expects log1p of CPM-per-10k. Rebuild that from raw counts so
    # the input matches what the model was trained on.
    src = ad.raw.to_adata() if ad.raw is not None else ad.copy()
    src.var_names_make_unique()
    sc.pp.normalize_total(src, target_sum=1e4)
    sc.pp.log1p(src)
    src.obs["minor_cell_state"] = ad.obs["minor_cell_state"].astype(str).values

    pred = celltypist.annotate(src, model=MODEL, majority_voting=False)
    lab = pred.predicted_labels
    col = "predicted_labels" if "predicted_labels" in lab.columns else lab.columns[0]

    df = pd.DataFrame({"minor_cell_state": src.obs["minor_cell_state"].values,
                       "celltypist": lab[col].values})
    df.to_csv(OUT / "celltypist_labels.csv", index=False)

    tab = pd.crosstab(df["minor_cell_state"], df["celltypist"], normalize="index") * 100

    def frac(state, kind):
        cols = [c for c in tab.columns if kind in str(c).lower()]
        return float(tab.loc[state, cols].sum()) if state in tab.index and cols else np.nan

    L = ["INDEPENDENT ANNOTATION BY CELLTYPIST - Reviewer 1 point R1.7",
         "=" * 92, ""]
    L.append(f"Model: {MODEL} (the model already used for the external cohorts)")
    L.append(f"Cells annotated: {len(df):,}")
    L.append("")
    L.append("Percentage of each cluster assigned to macrophage- and "
             "monocyte-family labels:")
    L.append("")
    L.append(f"  {'cluster':<34}{'macrophage':>12}{'monocyte':>11}{'top label':>34}")
    for state in tab.index:
        top = tab.loc[state].idxmax()
        L.append(f"  {state:<34}{frac(state, 'macrophage'):>11.1f}%"
                 f"{frac(state, 'monocyte'):>10.1f}%"
                 f"{str(top)[:32]:>34}")
    L.append("")

    if TARGET in tab.index:
        mac_t, mono_t = frac(TARGET, "macrophage"), frac(TARGET, "monocyte")
        mac_ref = frac(MAC_POLE, "macrophage")
        mono_ref = frac(MONO_POLE, "monocyte")
        L.append("INTERPRETATION")
        L.append("-" * 92)
        L.append(f"  Reference states in the same object: {MAC_POLE} is "
                 f"{mac_ref:.1f}% macrophage-labelled;")
        L.append(f"  {MONO_POLE} is {mono_ref:.1f}% monocyte-labelled.")
        L.append(f"  {TARGET} is {mac_t:.1f}% macrophage and {mono_t:.1f}% monocyte.")
        if mac_t > mono_t:
            L.append("  An external reference therefore assigns the IL-1b+ cluster")
            L.append("  predominantly to the macrophage family, which supports the")
            L.append("  published annotation while remaining consistent with the")
            L.append("  intermediate marker-score position: a monocyte-derived")
            L.append("  macrophage.")
        else:
            L.append("  An external reference assigns the IL-1b+ cluster")
            L.append("  predominantly to the monocyte family, consistent with the")
            L.append("  reviewer's alternative and with the intermediate lineage")
            L.append("  index; the population is named accordingly.")

    report = "\n".join(L)
    (OUT / "celltypist_report.txt").write_text(report, encoding="utf-8")
    print(report)
    _panel(tab)


def _panel(tab):
    d = (S11 / "S11_D"); d.mkdir(parents=True, exist_ok=True)
    keep = tab.loc[:, tab.max(axis=0) >= 5]
    keep = keep.reindex(sorted(keep.index))
    fig, ax = plt.subplots(figsize=(9.0 * SCALE * CM, 4.6 * SCALE * CM))
    bottom = np.zeros(len(keep))
    cmap = plt.get_cmap("tab20")
    for i, c in enumerate(keep.columns):
        ax.barh(np.arange(len(keep)), keep[c], left=bottom, height=0.7,
                color=cmap(i % 20), edgecolor="white", linewidth=0.4,
                label=str(c)[:28])
        bottom += keep[c].values
    ax.set_yticks(np.arange(len(keep)))
    ax.set_yticklabels([DISPLAY.get(s, s).replace("_", " ") for s in keep.index],
                       fontsize=5 * SCALE)
    ax.set_xlabel("% of cells assigned by CellTypist", fontsize=6 * SCALE)
    ax.set_xlim(0, 100)
    ax.tick_params(axis="both", labelsize=5 * SCALE, width=0.8, length=3)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
              fontsize=4.5 * SCALE, labelspacing=0.35)
    fig.subplots_adjust(left=0.26, right=0.66, top=0.97, bottom=0.16)
    stem = d / "S11_D_celltypist_annotation"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
