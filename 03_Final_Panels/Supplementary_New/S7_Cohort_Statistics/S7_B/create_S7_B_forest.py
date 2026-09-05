"""
Supplementary Figure S7B - two-sided sensitivity analysis.

Reviewer 1 asked that every directional test be repeated two-sided with exact P
values, and that the conclusions be shown to hold. This forest plot puts all 16
comparisons on a single standardised scale (Hedges' g with a percentile
bootstrap 95% CI), so the reader can see the direction and magnitude of every
effect independently of the choice of tail, alongside the exact two-sided P.

Input : 04_Revision_Analyses/02_R1.3_TwoSided_Stats_Sweep/outputs/twosided_sweep.csv
Output: S7_B_twosided_forest.{svg,pdf,png}
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import NEW_ANALYSES  # noqa: E402

OUT_DIR = Path(__file__).parent
SWEEP = (NEW_ANALYSES / "02_R1.3_TwoSided_Stats_Sweep" / "outputs"
         / "twosided_sweep.csv")

SCALE = 4
CM = 1 / 2.54
DPI = 300
PANEL_W_CM = 11.0 * SCALE
PANEL_H_CM = 8.0 * SCALE

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

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


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

    fig, (ax, axp) = plt.subplots(
        1, 2, figsize=(PANEL_W_CM * CM, PANEL_H_CM * CM),
        gridspec_kw={"width_ratios": [3.0, 1.0], "wspace": 0.04})

    # A few n = 4 comparisons have bootstrap upper limits beyond g = 8, which
    # would flatten every other interval. The axis is clipped and those
    # intervals are drawn with an arrowhead instead.
    X_MAX = 5.0
    X_MIN = -1.0

    ax.axvline(0, color="#999999", linewidth=0.8, linestyle="--", zorder=1)
    for i, (_, r) in enumerate(df.iterrows()):
        c = FAMILY_COLOR.get(r["family"], "#444444")
        lo = max(r["g_ci95_lo"], X_MIN)
        hi = min(r["g_ci95_hi"], X_MAX)
        ax.plot([lo, hi], [y[i]] * 2, color=c, linewidth=1.2,
                solid_capstyle="round", zorder=2)
        if r["g_ci95_hi"] > X_MAX:
            ax.plot(X_MAX, y[i], marker=">", color=c, markersize=3,
                    clip_on=False, zorder=3)
        if r["g_ci95_lo"] < X_MIN:
            ax.plot(X_MIN, y[i], marker="<", color=c, markersize=3,
                    clip_on=False, zorder=3)
        ax.scatter(min(max(r["hedges_g"], X_MIN), X_MAX), y[i], s=26, c=c,
                   edgecolors="white", linewidths=0.4, zorder=4)
    ax.set_xlim(X_MIN - 0.15, X_MAX + 0.15)

    ax.set_yticks(y)
    ax.set_yticklabels(df["row_label"], fontsize=5 * SCALE)
    ax.set_xlabel("Standardised effect size (Hedges' $g$), 95% bootstrap CI",
                  fontsize=6 * SCALE)
    ax.tick_params(axis="x", labelsize=5.5 * SCALE, width=0.5, length=2)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.8, len(df) - 0.2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)

    # Right-hand column: exact P values, both tails.
    axp.set_ylim(ax.get_ylim())
    axp.axis("off")
    axp.text(0.0, len(df) - 0.2, "$P$ 1-tail", fontsize=5.5 * SCALE,
             ha="left", va="bottom", style="italic")
    axp.text(0.52, len(df) - 0.2, "$P$ 2-tail", fontsize=5.5 * SCALE,
             ha="left", va="bottom", style="italic")
    for i, (_, r) in enumerate(df.iterrows()):
        axp.text(0.0, y[i], f"{r['p_one_tailed']:.3f}",
                 fontsize=5 * SCALE, ha="left", va="center", color="#666666")
        bold = "bold" if r["p_two_tailed"] < 0.05 else "normal"
        axp.text(0.52, y[i], f"{r['p_two_tailed']:.3f}", fontsize=5 * SCALE,
                 ha="left", va="center", fontweight=bold,
                 color="#000000" if r["p_two_tailed"] < 0.05 else "#666666")
    axp.set_xlim(-0.05, 1.05)

    handles = [plt.Line2D([], [], color=c, linewidth=1.2, marker="o",
                          markersize=3, label=f)
               for f, c in FAMILY_COLOR.items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.075),
              ncol=4, frameon=False, fontsize=5.5 * SCALE, handletextpad=0.5,
              columnspacing=1.6)

    n_ci = int(((df["g_ci95_lo"] > 0) | (df["g_ci95_hi"] < 0)).sum())
    fig.suptitle(
        f"All {len(df)} directional comparisons repeated two-sided; "
        f"{n_ci}/{len(df)} retain a 95% CI excluding zero",
        fontsize=6.5 * SCALE, y=0.985)

    fig.subplots_adjust(left=0.46, right=0.97, top=0.93, bottom=0.13)
    stem = OUT_DIR / "S7_B_twosided_forest"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
