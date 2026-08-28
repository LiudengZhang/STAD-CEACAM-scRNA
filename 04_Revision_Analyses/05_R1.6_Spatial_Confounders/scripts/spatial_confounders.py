"""
WP6 / Reviewer 1 point R1.6.

  "...this pattern may reflect several confounding factors, including tumor
   architecture, histological subtype, mucinous differentiation, or anatomical
   compartmentalization, rather than active immune exclusion. Controlling for
   total epithelial content is useful but insufficient to demonstrate that
   CEACAM expression itself drives immune exclusion."

What can be tested with the data in hand:
  1. Whether the CEACAM-distance association survives adjustment for BOTH total
     epithelial content and local epithelial density, in a mixed model with a
     random intercept per Visium sample.
  2. Whether it survives within strata of epithelial density, which removes the
     density confound non-parametrically rather than by assuming linearity.
  3. Whether it is reproducible across samples, i.e. how many of the 10 samples
     show the same sign.

What cannot be tested, and is stated as such: GSE251950 carries no annotation
for Lauren type, mucinous differentiation or tumour architecture, and the
spatial cohort is independent of the ICB-treated patients. The claim in the
manuscript is therefore reduced from immune exclusion being driven by CEACAM to
CEACAM-high regions co-localising with an immune-excluded architecture.

Input : Round_5/02_Preparation_for_Panels/Spatial/CEACAM_Deconvolution/spot_data.csv
Output: spatial_adjusted_models.csv, spatial_stratified.csv,
        spatial_confounders_report.txt, panels S9_A and S9_B
"""

from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA, REVISED_PANELS  # noqa: E402

warnings.filterwarnings("ignore")

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S9 = REVISED_PANELS / "Supplementary_New" / "S9_Mechanism_Specificity"

SCALE, CM, DPI = 4, 1 / 2.54, 300
OUTCOMES = {
    "distance_to_immune": "Distance to immune-rich regions",
    "distance_to_stroma": "Distance to stroma",
}
# Each model adds one adjustment, so the reader can see the CEACAM coefficient
# move (or not) as confounders are introduced.
MODELS = {
    "unadjusted": "{y} ~ CEACAM_ratio",
    "+ total epithelial content": "{y} ~ CEACAM_ratio + Total_Epi",
    "+ local epithelial density": (
        "{y} ~ CEACAM_ratio + Total_Epi + neighborhood_epi_density"),
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def main():
    spot = pd.read_csv(SPATIAL_SPOT_DATA)
    spot = spot.dropna(subset=list(OUTCOMES) + ["CEACAM_ratio", "Total_Epi",
                                                "neighborhood_epi_density", "sample"])
    # Distances are in pixels and vary in scale between samples; standardising
    # within sample makes the coefficients comparable and the random intercept
    # interpretable.
    for y in OUTCOMES:
        spot[y] = spot.groupby("sample")[y].transform(
            lambda s: (s - s.mean()) / s.std(ddof=1))

    rows = []
    for y, ylabel in OUTCOMES.items():
        for name, formula in MODELS.items():
            m = smf.mixedlm(formula.format(y=y), spot, groups=spot["sample"]).fit()
            ci = m.conf_int().loc["CEACAM_ratio"]
            rows.append(dict(
                outcome=ylabel, model=name,
                ceacam_coef=float(m.params["CEACAM_ratio"]),
                ceacam_se=float(m.bse["CEACAM_ratio"]),
                ceacam_p=float(m.pvalues["CEACAM_ratio"]),
                ci_low=float(ci[0]), ci_high=float(ci[1]),
                n_spots=int(m.nobs), n_samples=int(spot["sample"].nunique()),
            ))
    adjusted = pd.DataFrame(rows)
    adjusted.to_csv(OUT / "spatial_adjusted_models.csv", index=False)

    # ------------------------------------- stratified by epithelial density
    spot["density_tertile"] = pd.qcut(
        spot["neighborhood_epi_density"], 3,
        labels=["Low epithelial density", "Medium", "High epithelial density"])
    srows = []
    for y, ylabel in OUTCOMES.items():
        for tert, grp in spot.groupby("density_tertile", observed=True):
            m = smf.mixedlm(f"{y} ~ CEACAM_ratio", grp,
                            groups=grp["sample"]).fit()
            ci = m.conf_int().loc["CEACAM_ratio"]
            srows.append(dict(
                outcome=ylabel, stratum=str(tert),
                ceacam_coef=float(m.params["CEACAM_ratio"]),
                ceacam_p=float(m.pvalues["CEACAM_ratio"]),
                ci_low=float(ci[0]), ci_high=float(ci[1]), n_spots=int(m.nobs)))
    strat = pd.DataFrame(srows)
    strat.to_csv(OUT / "spatial_stratified.csv", index=False)

    # ------------------------------------------ per-sample sign consistency
    prows = []
    for y, ylabel in OUTCOMES.items():
        for s, grp in spot.groupby("sample", observed=True):
            if grp["CEACAM_ratio"].std(ddof=1) == 0:
                continue
            m = smf.ols(f"{y} ~ CEACAM_ratio + neighborhood_epi_density",
                        grp).fit()
            prows.append(dict(outcome=ylabel, sample=s,
                              ceacam_coef=float(m.params["CEACAM_ratio"]),
                              ceacam_p=float(m.pvalues["CEACAM_ratio"]),
                              n_spots=int(m.nobs)))
    per_sample = pd.DataFrame(prows)
    per_sample.to_csv(OUT / "spatial_per_sample.csv", index=False)

    # ------------------------------------------------------------- report
    L = ["SPATIAL CONFOUNDER ADJUSTMENT - Reviewer 1 point R1.6", "=" * 92, ""]
    L.append(f"Visium samples: {spot['sample'].nunique()}   "
             f"spots analysed: {len(spot):,}   (GSE251950, Korean Gut Atlas)")
    L.append("Outcomes are standardised within sample; models carry a random")
    L.append("intercept per sample.")
    L.append("")
    L.append("1. SEQUENTIAL ADJUSTMENT")
    L.append("-" * 92)
    for y in adjusted["outcome"].unique():
        L.append(f"  {y}")
        for _, r in adjusted[adjusted["outcome"] == y].iterrows():
            L.append(f"     {r['model']:<28} CEACAM beta = {r['ceacam_coef']:+.4f} "
                     f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]   "
                     f"P = {r['ceacam_p']:.3g}")
        L.append("")
    L.append("2. STRATIFIED BY LOCAL EPITHELIAL DENSITY")
    L.append("-" * 92)
    for y in strat["outcome"].unique():
        L.append(f"  {y}")
        for _, r in strat[strat["outcome"] == y].iterrows():
            L.append(f"     {r['stratum']:<28} CEACAM beta = {r['ceacam_coef']:+.4f} "
                     f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]   "
                     f"P = {r['ceacam_p']:.3g}   n = {r['n_spots']:,}")
        L.append("")
    L.append("3. CONSISTENCY ACROSS SAMPLES (density-adjusted, per sample)")
    L.append("-" * 92)
    for y in per_sample["outcome"].unique():
        sub = per_sample[per_sample["outcome"] == y]
        pos = int((sub["ceacam_coef"] > 0).sum())
        L.append(f"  {y}: positive coefficient in {pos}/{len(sub)} samples")
    L.append("")
    L.append("4. WHAT CANNOT BE CONTROLLED")
    L.append("-" * 92)
    L.append("  GSE251950 provides no annotation for Lauren classification, mucinous")
    L.append("  differentiation, or tumour architecture, and no spatial data exist for")
    L.append("  the chemo-immunotherapy cohort itself. Adjustment for epithelial")
    L.append("  content and density therefore addresses only part of the reviewer's")
    L.append("  concern. The manuscript claim is reduced accordingly, from CEACAM")
    L.append("  driving immune exclusion to CEACAM-high regions co-localising with an")
    L.append("  immune-excluded tissue architecture.")

    report = "\n".join(L)
    (OUT / "spatial_confounders_report.txt").write_text(report, encoding="utf-8")
    print(report)

    _panel_adjustment(adjusted)
    _panel_stratified(strat)


def _panel_adjustment(adjusted):
    d = (S9 / "S9_A"); d.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.0 * SCALE * CM, 3.6 * SCALE * CM))
    labels, y = [], []
    colors = {"Distance to immune-rich regions": "#7B3294",
              "Distance to stroma": "#1B7837"}
    k = 0
    for outcome in adjusted["outcome"].unique():
        for _, r in adjusted[adjusted["outcome"] == outcome].iterrows():
            ax.plot([r["ci_low"], r["ci_high"]], [k, k],
                    color=colors[outcome], linewidth=1.2)
            ax.scatter(r["ceacam_coef"], k, s=26, color=colors[outcome],
                       zorder=3, edgecolors="white", linewidths=0.4)
            # The three model labels repeat for each outcome, so the first row
            # of each block carries the outcome name to disambiguate them.
            short = ("immune" if "immune" in outcome else "stroma")
            labels.append(f"{short}:  {r['model']}")
            y.append(k)
            k += 1
        k += 0.8
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.5 * SCALE)
    ax.invert_yaxis()
    ax.set_xlabel("CEACAM ratio coefficient (SD of distance per unit)",
                  fontsize=6 * SCALE)
    ax.tick_params(axis="x", labelsize=5.5 * SCALE, width=0.8, length=3)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    handles = [plt.Line2D([], [], color=c, linewidth=1.2, marker="o",
                          markersize=3, label=o) for o, c in colors.items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.30),
              ncol=2, frameon=False, fontsize=5.5 * SCALE)
    fig.subplots_adjust(left=0.34, right=0.98, top=0.95, bottom=0.30)
    stem = d / "S9_A_spatial_adjusted_models"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


def _panel_stratified(strat):
    d = (S9 / "S9_B"); d.mkdir(parents=True, exist_ok=True)
    outcomes = list(strat["outcome"].unique())
    fig, axes = plt.subplots(1, len(outcomes),
                             figsize=(8.0 * SCALE * CM, 3.8 * SCALE * CM))
    for ax, outcome in zip(np.atleast_1d(axes), outcomes):
        sub = strat[strat["outcome"] == outcome]
        xs = np.arange(len(sub))
        ax.bar(xs, sub["ceacam_coef"], color="#7B3294", alpha=0.75,
               edgecolor="#444444", linewidth=0.5, width=0.6)
        ax.errorbar(xs, sub["ceacam_coef"],
                    yerr=[sub["ceacam_coef"] - sub["ci_low"],
                          sub["ci_high"] - sub["ceacam_coef"]],
                    fmt="none", ecolor="#333333", elinewidth=0.8, capsize=2)
        ax.axhline(0, color="#999999", linewidth=0.8)
        ax.set_xticks(xs)
        ax.set_xticklabels([s.replace(" epithelial density", "\nepith. density")
                            for s in sub["stratum"]], fontsize=5 * SCALE)
        ax.set_title(outcome, fontsize=6 * SCALE)
        ax.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    np.atleast_1d(axes)[0].set_ylabel("CEACAM ratio coefficient", fontsize=6 * SCALE)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.88, bottom=0.20, wspace=0.30)
    stem = d / "S9_B_spatial_density_stratified"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
