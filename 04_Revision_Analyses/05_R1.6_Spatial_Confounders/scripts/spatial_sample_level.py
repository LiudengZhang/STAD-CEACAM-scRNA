"""
Reviewer 1 point R1.6 - the stratified result with the section as the unit.

S9B fits one mixed model per density tertile over 7,777 spots with a random
intercept per section. The random intercept absorbs baseline differences between
sections but not differences in slope, so the P value is carried by the number of
spots rather than by the number of independent tissues, of which there are ten.
That is pseudo-replication and it is the first thing a statistically minded
reviewer would say.

This refits the same comparison one section at a time, inside each tertile, and
tests the ten coefficients as a single sample. Nothing about the model changes
except what counts as an observation.

The figure is for the response letter, not for the paper.

Inputs : Round_5/02_Preparation_for_Panels/Spatial/CEACAM_Deconvolution/spot_data.csv
Outputs: spatial_sample_level.csv, R_stratified_sample_level.[svg|pdf|png]
"""

from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA  # noqa: E402

warnings.filterwarnings("ignore")

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

SCALE, CM, DPI = 4, 1 / 2.54, 300
MIN_SPOTS = 30          # a section needs this many spots in a tertile to be fitted
TERTILES = ["Low", "Medium", "High"]
OUTCOMES = {"distance_to_immune": "Distance to immune-rich regions",
            "distance_to_stroma": "Distance to stroma"}
COLOR = "#7B3294"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def main():
    spot = pd.read_csv(SPATIAL_SPOT_DATA).dropna(
        subset=list(OUTCOMES) + ["CEACAM_ratio", "neighborhood_epi_density",
                                 "sample"])
    for y in OUTCOMES:
        spot[y] = spot.groupby("sample")[y].transform(
            lambda s: (s - s.mean()) / s.std(ddof=1))
    spot["tertile"] = pd.qcut(spot["neighborhood_epi_density"], 3, labels=TERTILES)

    per_section, summary = [], []
    for y, ylabel in OUTCOMES.items():
        for tert in TERTILES:
            sub = spot[spot["tertile"] == tert]
            coefs = {}
            for name, g in sub.groupby("sample"):
                if len(g) < MIN_SPOTS or g["CEACAM_ratio"].std(ddof=1) == 0:
                    continue
                coefs[name] = float(
                    smf.ols(f"{y} ~ CEACAM_ratio", g).fit().params["CEACAM_ratio"])
            c = pd.Series(coefs)
            for name, v in c.items():
                per_section.append(dict(outcome=ylabel, stratum=tert,
                                        section=name, ceacam_coef=v))
            summary.append(dict(
                outcome=ylabel, stratum=tert, n_sections=len(c),
                mean=float(c.mean()), median=float(c.median()),
                n_positive=int((c > 0).sum()),
                ci_lo=float(c.mean() - 1.96 * c.sem()),
                ci_hi=float(c.mean() + 1.96 * c.sem()),
                p_wilcoxon=float(stats.wilcoxon(c)[1]) if len(c) > 5 else np.nan,
                p_ttest=float(stats.ttest_1samp(c, 0).pvalue)))

    sec = pd.DataFrame(per_section)
    summ = pd.DataFrame(summary)
    sec.to_csv(OUT / "spatial_sample_level_sections.csv", index=False)
    summ.to_csv(OUT / "spatial_sample_level.csv", index=False)
    print(summ.round(4).to_string(index=False))

    _panel(sec, summ)


def _panel(sec, summ):
    fig, axes = plt.subplots(1, 2, figsize=(12.0 * SCALE * CM, 5.2 * SCALE * CM),
                             sharey=False)
    rng = np.random.default_rng(0)

    for ax, (y, ylabel) in zip(axes, OUTCOMES.items()):
        s = sec[sec["outcome"] == ylabel]
        m = summ[summ["outcome"] == ylabel].set_index("stratum")
        for i, tert in enumerate(TERTILES):
            v = s.loc[s["stratum"] == tert, "ceacam_coef"].values
            jitter = rng.uniform(-0.11, 0.11, len(v))
            ax.scatter(np.full(len(v), i) + jitter, v, s=14 * SCALE, c=COLOR,
                       alpha=0.75, edgecolors="white", linewidths=0.4 * SCALE,
                       zorder=3)
            r = m.loc[tert]
            ax.plot([i - 0.25, i + 0.25], [r["mean"]] * 2, color="#333333",
                    linewidth=1.2, zorder=4)
            ax.plot([i, i], [r["ci_lo"], r["ci_hi"]], color="#333333",
                    linewidth=0.9, zorder=4)
        ax.axhline(0, color="#999999", linewidth=0.8, zorder=1)
        # Headroom first, so the annotations sit above the points rather than on
        # them; they are placed in axes fraction, the points in data units.
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.20 * (hi - lo))
        for i, tert in enumerate(TERTILES):
            r = m.loc[tert]
            ax.annotate(
                f"P = {r['p_wilcoxon']:.3g}\n"
                f"{int(r['n_positive'])}/{int(r['n_sections'])} sections",
                (i, 0.99), xycoords=("data", "axes fraction"),
                ha="center", va="top", fontsize=5 * SCALE, color="#333333")
        ax.set_xticks(range(3))
        ax.set_xticklabels(["Low", "Medium", "High"], fontsize=6 * SCALE)
        ax.set_xlabel("Local epithelial density", fontsize=6 * SCALE)
        ax.set_ylabel("CEACAM ratio coefficient, per section",
                      fontsize=6 * SCALE)
        ax.set_title(ylabel, fontsize=6.5 * SCALE)
        ax.set_xlim(-0.55, 2.55)
        ax.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    fig.subplots_adjust(left=0.09, right=0.985, top=0.86, bottom=0.17, wspace=0.30)
    stem = OUT / "R_stratified_sample_level"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
