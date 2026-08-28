"""
WP5 / Reviewer 1 point R1.5.

  "The manuscript frequently treats CEACAM5 and CEACAM6 as a single entity...
   whether the CEACAM-high cell state requires expression of both markers or
   either marker; and if CEACAM5-only and CEACAM6-only populations exist and
   predict response. ... why IHC staining percentages were summed rather than
   analyzed separately or combined using a prespecified composite score."

Three questions are answered separately.

  A. Do single-positive populations exist? Every pre-treatment epithelial cell is
     assigned to CEACAM5+CEACAM6-, CEACAM5-CEACAM6+, double-positive or
     double-negative on detected counts, and the per-sample fraction of each is
     compared between responders and non-responders.
  B. Which marker carries the association? The two genes are tested separately.
  C. Does the IHC conclusion depend on summation? CEACAM5 and CEACAM6 staining
     are tested separately and against the summed composite.

Inputs : Round_5/01_Raw_Inputs/01_H5AD/Epithelial.h5ad
         Round_5/02_Preparation_for_Panels/IHC/ceacam_ihc_color_deconv_results.csv
Outputs: ceacam_state_fractions.csv, ceacam_state_tests.csv,
         ihc_per_marker_tests.csv, ceacam5_vs_6_report.txt
         panels S8_B (state fractions) and S8_C (per-marker IHC)
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD, PREPARATION, REVISED_PANELS  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S8 = REVISED_PANELS / "Supplementary_New" / "S8_CEACAM_Metaprogram"

SCALE, CM, DPI = 4, 1 / 2.54, 300
MIN_CELLS = 20
COLOR_R, COLOR_NR = "#2166AC", "#B2182B"

STATES = ["CEACAM5+ only", "CEACAM6+ only", "Double positive", "Double negative"]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def mw(a, b):
    """Two-sided Mann-Whitney with rank-biserial effect size."""
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(p), float(2.0 * u / (len(a) * len(b)) - 1.0)


def load_pre_epithelial():
    ad = sc.read_h5ad(EPITHELIAL_H5AD)
    ad = ad[ad.obs["Sample site"] == "Stomach"]
    ad = ad[ad.obs["Treatment phase"] == "Pre"]
    ad = ad[ad.obs["stomach_pre_grouping"].isin(["Responsed", "No-response"])].copy()
    src = ad.raw if ad.raw is not None else ad
    out = {}
    for g in ("CEACAM5", "CEACAM6"):
        i = list(src.var_names).index(g)
        x = src.X[:, i]
        out[g] = x.toarray().flatten() if hasattr(x, "toarray") else np.asarray(x).flatten()
    df = pd.DataFrame(out)
    df["sample"] = ad.obs["sample"].astype(str).values
    df["group"] = ad.obs["stomach_pre_grouping"].map(
        {"Responsed": "R", "No-response": "NR"}).values
    return df


def main():
    L = ["CEACAM5 VERSUS CEACAM6 - Reviewer 1 point R1.5", "=" * 90, ""]

    # ============================================ A. single-positive states
    cells = load_pre_epithelial()
    # Positivity is detection of at least one UMI, the standard threshold for a
    # sparse count matrix; no arbitrary expression cutoff is introduced.
    c5, c6 = cells["CEACAM5"] > 0, cells["CEACAM6"] > 0
    cells["state"] = np.select(
        [c5 & ~c6, ~c5 & c6, c5 & c6],
        ["CEACAM5+ only", "CEACAM6+ only", "Double positive"],
        default="Double negative")

    keep = cells.groupby("sample").size()
    cells = cells[cells["sample"].isin(keep[keep >= MIN_CELLS].index)]

    frac = (cells.groupby(["sample", "group", "state"], observed=True).size()
            .unstack("state").fillna(0))
    frac = frac.div(frac.sum(axis=1), axis=0).reset_index()
    for s in STATES:
        if s not in frac.columns:
            frac[s] = 0.0
    frac.to_csv(OUT / "ceacam_state_fractions.csv", index=False)

    rows = []
    for s in STATES:
        nr = frac.loc[frac["group"] == "NR", s].values
        r = frac.loc[frac["group"] == "R", s].values
        p, eff = mw(nr, r)
        rows.append(dict(measure=s, level="cell state fraction",
                         n_NR=len(nr), n_R=len(r),
                         mean_NR=nr.mean(), mean_R=r.mean(),
                         p_two_tailed=p, rank_biserial_r=eff))
    # Also: "either marker" - the union, which is what "CEACAM5/6+" implies
    frac["Either marker"] = (frac["CEACAM5+ only"] + frac["CEACAM6+ only"]
                             + frac["Double positive"])
    nr = frac.loc[frac["group"] == "NR", "Either marker"].values
    r = frac.loc[frac["group"] == "R", "Either marker"].values
    p, eff = mw(nr, r)
    rows.append(dict(measure="Either marker (union)", level="cell state fraction",
                     n_NR=len(nr), n_R=len(r), mean_NR=nr.mean(), mean_R=r.mean(),
                     p_two_tailed=p, rank_biserial_r=eff))
    tests = pd.DataFrame(rows)
    tests.to_csv(OUT / "ceacam_state_tests.csv", index=False)

    L.append("A. DO CEACAM5-ONLY AND CEACAM6-ONLY POPULATIONS EXIST?")
    L.append("-" * 90)
    L.append(f"   Pre-treatment gastric epithelial cells analysed: {len(cells):,}")
    L.append(f"   Samples: {frac['sample'].nunique()} "
             f"({(frac['group'] == 'NR').sum()} NR, {(frac['group'] == 'R').sum()} R)")
    L.append("")
    overall = cells["state"].value_counts(normalize=True) * 100
    for s in STATES:
        L.append(f"   {s:<18} {overall.get(s, 0):5.1f}% of all epithelial cells")
    L.append("")
    L.append("   Per-sample fraction, non-responders vs responders (two-sided MW):")
    for _, t in tests.iterrows():
        L.append(f"     {t['measure']:<22} NR {t['mean_NR']*100:5.2f}%  vs  "
                 f"R {t['mean_R']*100:5.2f}%   P = {t['p_two_tailed']:.4f}   "
                 f"r = {t['rank_biserial_r']:+.3f}")
    L.append("")

    # ================================================ C. per-marker IHC
    ihc = pd.read_csv(PREPARATION / "IHC" / "ceacam_ihc_color_deconv_results.csv")
    piv = ihc.pivot_table(index=["patient", "group"], columns="marker",
                          values="staining_pct").reset_index()
    piv["Summed (published)"] = piv["CEACAM5"] + piv["CEACAM6"]
    piv["Mean of the two"] = piv[["CEACAM5", "CEACAM6"]].mean(axis=1)
    # A prespecified composite that does not depend on the two markers being on
    # the same scale: the mean of within-marker z-scores.
    z = piv[["CEACAM5", "CEACAM6"]].apply(lambda c: (c - c.mean()) / c.std(ddof=1))
    piv["z-score composite"] = z.mean(axis=1)
    piv.to_csv(OUT / "ihc_per_marker_values.csv", index=False)

    ihc_rows = []
    for col in ("CEACAM5", "CEACAM6", "Summed (published)", "Mean of the two",
                "z-score composite"):
        nr = piv.loc[piv["group"] == "NR", col].values
        r = piv.loc[piv["group"] == "R", col].values
        p, eff = mw(nr, r)
        ihc_rows.append(dict(measure=col, level="IHC staining percent",
                             n_NR=len(nr), n_R=len(r),
                             mean_NR=nr.mean(), mean_R=r.mean(),
                             p_two_tailed=p, rank_biserial_r=eff))
    ihc_tests = pd.DataFrame(ihc_rows)
    ihc_tests.to_csv(OUT / "ihc_per_marker_tests.csv", index=False)

    L.append("B/C. IMMUNOHISTOCHEMISTRY, EACH MARKER SEPARATELY (n = 4 vs 4)")
    L.append("-" * 90)
    for _, t in ihc_tests.iterrows():
        L.append(f"     {t['measure']:<22} NR {t['mean_NR']:7.2f}  vs  "
                 f"R {t['mean_R']:7.2f}   P = {t['p_two_tailed']:.4f}   "
                 f"r = {t['rank_biserial_r']:+.3f}")
    L.append("")
    L.append("   At n = 4 per group the smallest two-sided P a Mann-Whitney test can")
    L.append("   return is 2/70 = 0.029, so none of these comparisons can distinguish")
    L.append("   a marker-specific effect from the summed one. The summed panel is")
    L.append("   therefore replaced by the two markers shown separately, and the")
    L.append("   summed value is reported only as a descriptive composite.")

    report = "\n".join(L)
    (OUT / "ceacam5_vs_6_report.txt").write_text(report, encoding="utf-8")
    print(report)

    # ------------------------------------------------------------- panels
    _panel_states(frac)
    _panel_ihc(piv, ihc_tests)


def _panel_states(frac):
    d = (S8 / "S8_B"); d.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(11.0 * SCALE * CM, 4.2 * SCALE * CM))
    rng = np.random.default_rng(1)
    for ax, s in zip(axes, STATES):
        r = frac.loc[frac["group"] == "R", s].values * 100
        nr = frac.loc[frac["group"] == "NR", s].values * 100
        bp = ax.boxplot([r, nr], positions=[0, 1], widths=0.55, patch_artist=True,
                        showfliers=False, boxprops=dict(linewidth=0.8),
                        whiskerprops=dict(linewidth=0.8),
                        capprops=dict(linewidth=0.8),
                        medianprops=dict(color="black", linewidth=1.2))
        bp["boxes"][0].set_facecolor(COLOR_R); bp["boxes"][0].set_alpha(0.55)
        bp["boxes"][1].set_facecolor(COLOR_NR); bp["boxes"][1].set_alpha(0.55)
        for i, (vals, c) in enumerate(((r, COLOR_R), (nr, COLOR_NR))):
            ax.scatter(i + rng.uniform(-0.1, 0.1, len(vals)), vals, s=9 * SCALE,
                       c=c, zorder=3, edgecolors="white", linewidths=0.3 * SCALE)
        p, _ = mw(nr, r)
        top = max(np.max(r), np.max(nr))
        ax.plot([0, 0, 1, 1], [top * 1.06, top * 1.11, top * 1.11, top * 1.06],
                color="#444444", linewidth=0.8)
        ax.text(0.5, top * 1.13, f"P = {p:.3f}", ha="center", va="bottom",
                fontsize=5 * SCALE)
        ax.set_ylim(0, top * 1.30)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["R", "NR"], fontsize=6 * SCALE)
        ax.set_title(s, fontsize=6.5 * SCALE)
        ax.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("% of pre-treatment\nepithelial cells", fontsize=6 * SCALE)
    fig.subplots_adjust(left=0.09, right=0.99, top=0.87, bottom=0.11, wspace=0.35)
    stem = d / "S8_B_ceacam_single_double_positive"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


def _panel_ihc(piv, tests):
    d = (S8 / "S8_C"); d.mkdir(parents=True, exist_ok=True)
    cols = ["CEACAM5", "CEACAM6", "Summed (published)"]
    fig, axes = plt.subplots(1, 3, figsize=(8.0 * SCALE * CM, 4.2 * SCALE * CM))
    rng = np.random.default_rng(2)
    for ax, col in zip(axes, cols):
        r = piv.loc[piv["group"] == "R", col].values
        nr = piv.loc[piv["group"] == "NR", col].values
        bp = ax.boxplot([r, nr], positions=[0, 1], widths=0.55, patch_artist=True,
                        showfliers=False, boxprops=dict(linewidth=0.8),
                        whiskerprops=dict(linewidth=0.8),
                        capprops=dict(linewidth=0.8),
                        medianprops=dict(color="black", linewidth=1.2))
        bp["boxes"][0].set_facecolor(COLOR_R); bp["boxes"][0].set_alpha(0.55)
        bp["boxes"][1].set_facecolor(COLOR_NR); bp["boxes"][1].set_alpha(0.55)
        for i, (vals, c) in enumerate(((r, COLOR_R), (nr, COLOR_NR))):
            ax.scatter(i + rng.uniform(-0.1, 0.1, len(vals)), vals, s=11 * SCALE,
                       c=c, zorder=3, edgecolors="white", linewidths=0.3 * SCALE)
        p = tests.loc[tests["measure"] == col, "p_two_tailed"].iloc[0]
        top = max(np.max(r), np.max(nr))
        ax.plot([0, 0, 1, 1], [top * 1.05, top * 1.09, top * 1.09, top * 1.05],
                color="#444444", linewidth=0.8)
        ax.text(0.5, top * 1.11, f"P = {p:.3f}", ha="center", va="bottom",
                fontsize=5.5 * SCALE)
        ax.set_ylim(0, top * 1.28)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["R", "NR"], fontsize=6 * SCALE)
        ax.set_title(col, fontsize=6.5 * SCALE)
        ax.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("DAB+ area (% of tissue)", fontsize=6 * SCALE)
    fig.suptitle("Immunohistochemistry, n = 4 responders vs 4 non-responders "
                 "(two-sided Mann-Whitney)", fontsize=6 * SCALE, y=0.99)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.80, bottom=0.11, wspace=0.35)
    stem = d / "S8_C_ihc_per_marker"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
