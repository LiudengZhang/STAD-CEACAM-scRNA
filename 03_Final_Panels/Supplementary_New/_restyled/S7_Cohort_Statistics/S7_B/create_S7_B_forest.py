"""
Supplementary Figure S7B, RESTYLED (Version B) - two-sided sensitivity analysis.

Copy of Supplementary_New/S7_Cohort_Statistics/S7_B/create_S7_B_forest.py with
the type taken from cnsplots via 00_Config/panel_style_cns.py, the canvas set to
the millimetre box the panel prints in, and the point-specified marks rescaled
so their size relative to the type is unchanged. Every value read, flipped,
sorted and plotted is the same code.

Version A drew 44 x 32 cm at SCALE 4 and the assembler fitted it into a
171 x 74 mm box at 0.2313, so its row labels - set at 5 * SCALE - printed at
4.63 pt. Sixteen rows of labels up to about sixty characters need roughly 72 mm
of paper at 7 pt, against 47 mm at 4.63 pt, so the panel is taller and its label
column wider. The sixteen rows, their order, their effect sizes, their
intervals and both P columns are unchanged.

Reviewer 1 asked that every directional test be repeated two-sided with exact P
values, and that the conclusions be shown to hold. This forest plot puts all 16
comparisons on a single standardised scale (Hedges' g with a percentile
bootstrap 95% CI), so the reader can see the direction and magnitude of every
effect independently of the choice of tail, alongside the exact two-sided P.

Input : 02_New_Analyses/02_R1.3_TwoSided_Stats_Sweep/outputs/twosided_sweep.csv
Output: S7_B_twosided_forest.{svg,pdf,png}
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
SWEEP = (NEW_ANALYSES / "02_R1.3_TwoSided_Stats_Sweep" / "outputs"
         / "twosided_sweep.csv")

PANEL_W_MM = 171.0
PANEL_H_MM = 104.0

A_FIT = 0.2313
A_TYPE = 5.0 * 4 * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2

FAMILY_COLOR = {
    "CEACAM": "#2166AC",
    "Metaprogram": "#4393C3",
    "PD-L1": "#B2182B",
    "Spatial": "#7B3294",
    "Cytokine": "#D6604D",
    "Regulon": "#1B7837",
}

# The sweep keys panels by the figure code as it stands now, which was
# re-lettered after submission. Readers see the printed letters, so the plot
# does too. Kept in step with SUBMITTED_PANEL in 04_Manuscript_R1/04_Tables/
# build_tables.py, which does the same translation for Table S6.
SUBMITTED_PANEL = {
    "Fig 2 (N1)": "Fig 2K", "Fig 2 (N2)": "Fig 2K", "Fig 2 (O1)": "Fig 2L",
    "Fig 2 (O2)": "Fig 2L", "Fig 2 (Q)": "Fig 2N", "Fig 2 (K)": "Fig 2H",
    "Fig 5 (J)": "Fig 5J", "Fig 5 (K)": "Fig 5J", "Fig 5 (L)": "Fig 5J",
    "Fig 5 (M)": "Fig 5J", "Fig 3 (H)": "Fig 3H", "Fig 3 (I)": "Fig 3I",
    "Fig 3 (J)": "Fig 3J", "Fig 5 (regulon-BACH1)": "Fig 5D",
    "Fig 5 (regulon-NFKB1)": "Fig 5E", "Fig 5 (IL6-CD4)": "Fig 5L",
}


def short_label(row):
    """Panel tag plus a trimmed analysis name, so rows stay one line."""
    name = row["analysis"]
    name = name.replace(" expression, ", ", ").replace("CEACAM-high vs CEACAM-low spots", "high vs low")
    name = name.replace(", pre-treatment R vs all other groups", ", pre-R vs rest")
    if len(name) > 52:
        name = name[:49] + "..."
    panel = SUBMITTED_PANEL.get(row["panel"], row["panel"])
    return f"{panel}  {name}"


def main():
    style.apply()
    df = pd.read_csv(SWEEP)
    # Metaprogram and regulon effects are negative by construction, because the
    # responder group is the one being contrasted against the rest; flip them so
    # every row reads "higher in non-responders / CEACAM-high" and the zero line
    # means the same thing throughout.
    flip = df["family"].isin(("Metaprogram", "Regulon"))
    for col in ("hedges_g", "g_ci95_lo", "g_ci95_hi"):
        df.loc[flip, col] = -df.loc[flip, col]
    df.loc[flip, ["g_ci95_lo", "g_ci95_hi"]] = df.loc[
        flip, ["g_ci95_hi", "g_ci95_lo"]].values

    df = df.sort_values(["family", "panel"], ascending=[True, True]).reset_index(drop=True)
    df["row_label"] = df.apply(short_label, axis=1)
    y = np.arange(len(df))[::-1]

    fig, (ax, axp) = style.subplots_mm(
        PANEL_W_MM, PANEL_H_MM, 1, 2,
        gridspec_kw={"width_ratios": [3.0, 1.15], "wspace": 0.04})

    # A few n = 4 comparisons have bootstrap upper limits beyond g = 8, which
    # would flatten every other interval. The axis is clipped and those
    # intervals are drawn with an arrowhead instead.
    X_MAX = 5.0
    X_MIN = -1.0

    ax.axvline(0, color="#999999", linewidth=0.8 * MARK, linestyle="--",
                zorder=1)
    for i, (_, r) in enumerate(df.iterrows()):
        c = FAMILY_COLOR.get(r["family"], "#444444")
        lo = max(r["g_ci95_lo"], X_MIN)
        hi = min(r["g_ci95_hi"], X_MAX)
        ax.plot([lo, hi], [y[i]] * 2, color=c, linewidth=1.2 * MARK,
                solid_capstyle="round", zorder=2)
        if r["g_ci95_hi"] > X_MAX:
            ax.plot(X_MAX, y[i], marker=">", color=c, markersize=3 * MARK,
                    clip_on=False, zorder=3)
        if r["g_ci95_lo"] < X_MIN:
            ax.plot(X_MIN, y[i], marker="<", color=c, markersize=3 * MARK,
                    clip_on=False, zorder=3)
        ax.scatter(min(max(r["hedges_g"], X_MIN), X_MAX), y[i], s=26 * AREA,
                   c=c, edgecolors="white", linewidths=0.4 * MARK, zorder=4)
    ax.set_xlim(X_MIN - 0.15, X_MAX + 0.15)

    ax.set_yticks(y)
    ax.set_yticklabels(df["row_label"])
    ax.set_xlabel("Standardised effect size (Hedges' $g$), 95% bootstrap CI")
    ax.tick_params(axis="x", width=0.6, length=2)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.8, len(df) - 0.2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Right-hand column: exact P values, both tails.
    axp.set_ylim(ax.get_ylim())
    axp.axis("off")
    axp.text(0.0, len(df) - 0.2, "$P$ 1-tail", fontsize=style.tick_pt(),
             ha="left", va="bottom", style="italic")
    axp.text(0.52, len(df) - 0.2, "$P$ 2-tail", fontsize=style.tick_pt(),
             ha="left", va="bottom", style="italic")
    for i, (_, r) in enumerate(df.iterrows()):
        axp.text(0.0, y[i], f"{r['p_one_tailed']:.3f}",
                 fontsize=style.tick_pt(), ha="left", va="center",
                 color="#666666")
        bold = "bold" if r["p_two_tailed"] < 0.05 else "normal"
        axp.text(0.52, y[i], f"{r['p_two_tailed']:.3f}",
                 fontsize=style.tick_pt(), ha="left", va="center",
                 fontweight=bold,
                 color="#000000" if r["p_two_tailed"] < 0.05 else "#666666")
    axp.set_xlim(-0.05, 1.05)

    # cnsplots draws legend keys at legend.markerscale (0.5), which Version A
    # did not do, so the key is sized to print at the same diameter as the mark
    # it stands for once that scale is applied.
    key = 1.0 / plt.rcParams["legend.markerscale"]
    handles = [plt.Line2D([], [], color=c, linewidth=1.2 * MARK, marker="o",
                          markersize=(26 * AREA) ** 0.5 * key, label=f)
               for f, c in FAMILY_COLOR.items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.075),
              ncol=4, handletextpad=0.5, columnspacing=1.6)

    n_ci = int(((df["g_ci95_lo"] > 0) | (df["g_ci95_hi"] < 0)).sum())
    fig.suptitle(
        f"All {len(df)} directional comparisons repeated two-sided; "
        f"{n_ci}/{len(df)} retain a 95% CI excluding zero",
        fontsize=style.body_pt(), y=0.99)

    style.margins_mm(fig, left=70, right=16, top=6, bottom=15)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUT_DIR / "S7_B_twosided_forest")
    print(f"Saved S7_B_twosided_forest at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
