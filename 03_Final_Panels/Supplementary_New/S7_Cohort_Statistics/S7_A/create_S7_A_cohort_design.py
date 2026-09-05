"""
Supplementary Figure S7A - cohort design and longitudinal pairing.

Reviewer 1 asked how many samples were longitudinally paired and whether the
pre- and post-treatment groups come from the same patients. This panel answers
both visually: one row per patient, one column per timepoint, with the single
paired patient joined by a line.

Input : 04_Revision_Analyses/01_R1.3_Cohort_Pairing/outputs/cohort_audit.csv
        Round_5/04_Manuscript/04_Tables/ST1_patient_sample_characteristics.csv
Output: S7_A_cohort_design.{svg,pdf,png}
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import MANUSCRIPT, NEW_ANALYSES  # noqa: E402

OUT_DIR = Path(__file__).parent
SCALE = 4
CM = 1 / 2.54
DPI = 300

PANEL_W_CM = 9.0 * SCALE
PANEL_H_CM = 6.5 * SCALE

# Standard four-group palette (blue = pre, pink = post; lighter = R).
# The two blues and the two pinks are deliberately close, so responder status is
# additionally encoded by marker shape: circle = R, square = NR.
COLORS = {
    ("Pre", "R"): "#bde0fe",
    ("Pre", "NR"): "#a2d2ff",
    ("Post", "R"): "#ffcfd2",
    ("Post", "NR"): "#f1c0e8",
}
MARKERS = {"R": "o", "NR": "s"}
UNLABELLED = "#e8e8e8"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def main():
    st1 = pd.read_csv(MANUSCRIPT / "04_Tables" / "ST1_patient_sample_characteristics.csv")
    st1.columns = [c.strip() for c in st1.columns]
    stomach = st1[st1["Anatomical site"] == "Stomach"].copy()

    # Patient order: labelled patients first (grouped by phase and response),
    # then the unlabelled ones, so the panel reads top-to-bottom as a cohort.
    stomach["label"] = stomach["R/NR Grouping"].fillna("Unlabelled")
    order_key = stomach.groupby("Patient ID").agg(
        has_pre=("Treatment phase", lambda s: "Pre" in set(s)),
        has_post=("Treatment phase", lambda s: "Post" in set(s)),
        labelled=("R/NR Grouping", lambda s: s.notna().any()),
        resp=("R/NR Grouping", lambda s: sorted(set(s.dropna()))[:1] or ["Z"]),
    )
    order_key["resp"] = order_key["resp"].apply(lambda v: v[0])
    order_key["paired"] = order_key["has_pre"] & order_key["has_post"]
    order_key["pnum"] = [int(p[1:]) for p in order_key.index]
    order_key = order_key.sort_values(
        ["paired", "labelled", "has_pre", "resp", "pnum"],
        ascending=[False, False, False, True, True],
    )
    patients = list(order_key.index)
    ypos = {p: len(patients) - i for i, p in enumerate(patients)}

    xpos = {"Pre": 0, "Post": 1}

    fig, ax = plt.subplots(figsize=(PANEL_W_CM * CM, PANEL_H_CM * CM))

    # Connector for any patient sampled at both timepoints.
    n_paired = 0
    for pid, row in order_key.iterrows():
        if row["paired"]:
            n_paired += 1
            ax.plot([0, 1], [ypos[pid]] * 2, color="#333333",
                    linewidth=1.2 * SCALE, zorder=1, solid_capstyle="round")

    for _, r in stomach.iterrows():
        phase = r["Treatment phase"]
        grp = r["R/NR Grouping"]
        color = COLORS.get((phase, grp), UNLABELLED)
        ax.scatter(xpos[phase], ypos[r["Patient ID"]], s=55 * SCALE,
                   c=color, marker=MARKERS.get(grp, "o"),
                   edgecolors="#444444", linewidths=0.5 * SCALE,
                   zorder=3, clip_on=False)

    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(0.2, len(patients) + 0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Pre-treatment", "Post-treatment"], fontsize=7 * SCALE)
    ax.set_yticks([ypos[p] for p in patients])
    ax.set_yticklabels(patients, fontsize=4.5 * SCALE)
    ax.set_ylabel("Patient", fontsize=7 * SCALE)
    ax.tick_params(axis="both", length=0, pad=3 * SCALE)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    # Legend
    handles = [
        plt.Line2D([], [], marker=MARKERS[k[1]], linestyle="none",
                   markersize=3.2 * SCALE, markerfacecolor=COLORS[k],
                   markeredgecolor="#444444", markeredgewidth=0.5 * SCALE,
                   label=f"{k[0]}-treatment {k[1]}")
        for k in [("Pre", "R"), ("Pre", "NR"), ("Post", "R"), ("Post", "NR")]
    ]
    handles.append(plt.Line2D([], [], marker="o", linestyle="none",
                              markersize=3.2 * SCALE, markerfacecolor=UNLABELLED,
                              markeredgecolor="#444444",
                              markeredgewidth=0.5 * SCALE,
                              label="No response designation"))
    handles.append(plt.Line2D([], [], color="#333333", linewidth=1.2 * SCALE,
                              label="Same patient, both timepoints"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.04, 1.0),
              frameon=False, fontsize=5.5 * SCALE, handletextpad=0.8,
              labelspacing=1.1, borderpad=0)

    n_pre = int((stomach["Treatment phase"] == "Pre").sum())
    n_post = int((stomach["Treatment phase"] == "Post").sum())
    fig.suptitle(
        f"{n_pre} pre-treatment and {n_post} post-treatment gastric specimens "
        f"from {stomach['Patient ID'].nunique()} patients\n"
        f"{n_paired} of {st1['Patient ID'].nunique()} patients were sampled at "
        f"both timepoints",
        fontsize=7 * SCALE, x=0.5, y=0.985, va="top")

    fig.subplots_adjust(left=0.10, right=0.58, top=0.86, bottom=0.07)
    stem = OUT_DIR / "S7_A_cohort_design"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}.[svg|pdf|png]  ({n_paired} paired patients)")


if __name__ == "__main__":
    main()
