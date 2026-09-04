"""
Supplementary Figure S7A, RESTYLED (Version B) - cohort design and longitudinal
pairing.

Copy of Supplementary_New/S7_Cohort_Statistics/S7_A/create_S7_A_cohort_design.py
with the type taken from cnsplots via 00_Config/panel_style_cns.py, the canvas
set to the millimetre box the panel prints in, and the point-specified marks
rescaled so their size relative to the type is unchanged. Every value read,
ordered and plotted is the same code.

S7A is the panel the assembler shrank least (fit 0.4538), so it is the one
figure in the set whose type was already close to the target: its patient labels
printed at 8.17 pt and its titles at 12.7 pt. Here the labels go to 7 pt and the
titles to 8 pt, so this panel's type gets *smaller*, not larger - which is the
point of setting one system rather than scaling each panel to taste.

Reviewer 1 asked how many samples were longitudinally paired and whether the
pre- and post-treatment groups come from the same patients. This panel answers
both visually: one row per patient, one column per timepoint, with the single
paired patient joined by a line.

Input : 02_New_Analyses/01_R1.3_Cohort_Pairing/outputs/cohort_audit.csv
        Round_5/04_Manuscript/04_Tables/ST1_patient_sample_characteristics.csv
Output: S7_A_cohort_design.{svg,pdf,png}
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[5] / "00_Config"))
from paths import MANUSCRIPT, NEW_ANALYSES  # noqa: E402
import panel_style_cns as style  # noqa: E402

OUT_DIR = Path(__file__).parent

PANEL_W_MM = 171.0
PANEL_H_MM = 134.0

A_FIT = 0.4538
A_TYPE = 4.5 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2

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


def main():
    style.apply()
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

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    # Connector for any patient sampled at both timepoints.
    n_paired = 0
    for pid, row in order_key.iterrows():
        if row["paired"]:
            n_paired += 1
            ax.plot([0, 1], [ypos[pid]] * 2, color="#333333",
                    linewidth=1.2 * 4 * MARK, zorder=1,
                    solid_capstyle="round")

    for _, r in stomach.iterrows():
        phase = r["Treatment phase"]
        grp = r["R/NR Grouping"]
        color = COLORS.get((phase, grp), UNLABELLED)
        ax.scatter(xpos[phase], ypos[r["Patient ID"]], s=55 * 4 * AREA,
                   c=color, marker=MARKERS.get(grp, "o"),
                   edgecolors="#444444", linewidths=0.5 * 4 * MARK,
                   zorder=3, clip_on=False)

    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(0.2, len(patients) + 0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Pre-treatment", "Post-treatment"],
                       fontsize=style.body_pt())
    ax.set_yticks([ypos[p] for p in patients])
    ax.set_yticklabels(patients)
    ax.set_ylabel("Patient")
    ax.tick_params(axis="both", length=0, pad=3 * 4 * MARK)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    # Legend
    # cnsplots halves legend keys (legend.markerscale 0.5); size each key so it
    # prints at the same diameter as the mark it stands for.
    key = 1.0 / plt.rcParams["legend.markerscale"]
    mk = (55 * 4 * AREA) ** 0.5 * key
    handles = [
        plt.Line2D([], [], marker=MARKERS[k[1]], linestyle="none",
                   markersize=mk, markerfacecolor=COLORS[k],
                   markeredgecolor="#444444", markeredgewidth=0.5 * 4 * MARK,
                   label=f"{k[0]}-treatment {k[1]}")
        for k in [("Pre", "R"), ("Pre", "NR"), ("Post", "R"), ("Post", "NR")]
    ]
    handles.append(plt.Line2D([], [], marker="o", linestyle="none",
                              markersize=mk, markerfacecolor=UNLABELLED,
                              markeredgecolor="#444444",
                              markeredgewidth=0.5 * 4 * MARK,
                              label="No response designation"))
    handles.append(plt.Line2D([], [], color="#333333",
                              linewidth=1.2 * 4 * MARK,
                              label="Same patient, both timepoints"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.04, 1.0),
              handletextpad=0.8, labelspacing=1.1, borderpad=0)

    n_pre = int((stomach["Treatment phase"] == "Pre").sum())
    n_post = int((stomach["Treatment phase"] == "Post").sum())
    fig.suptitle(
        f"{n_pre} pre-treatment and {n_post} post-treatment gastric specimens "
        f"from {stomach['Patient ID'].nunique()} patients\n"
        f"{n_paired} of {st1['Patient ID'].nunique()} patients were sampled at "
        f"both timepoints",
        fontsize=style.body_pt(), x=0.5, y=0.995, va="top")

    style.margins_mm(fig, left=15, right=48, top=12, bottom=11)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUT_DIR / "S7_A_cohort_design")
    print(f"Saved S7_A_cohort_design at {PANEL_W_MM} x {PANEL_H_MM} mm "
          f"({n_paired} paired patients)")


if __name__ == "__main__":
    main()
