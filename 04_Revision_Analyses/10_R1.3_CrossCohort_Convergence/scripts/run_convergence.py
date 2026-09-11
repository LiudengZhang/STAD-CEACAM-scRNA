"""
Reviewer 1 point R1.3 - what does the pre-treatment CEACAM5/6 result look like
when the cohorts are read together rather than one at a time?

The submitted paper reports three CEACAM measurements separately, each landing
near P = 0.057, and each therefore reads as unconvincing on its own. Two of them
are statistically independent and can be combined; the third cannot, and this
script is explicit about why.

  scRNA discovery cohort     4 NR vs 4 R, pre-treatment gastric epithelium
  PRJEB25780 (TIGER)         33 NR vs 12 R, deconvolved bulk epithelium
  immunohistochemistry       the SAME eight biopsies as the scRNA cohort: each
                             block is the fixed portion of the specimen that was
                             sequenced, so it is an orthogonal measurement of one
                             specimen, not an independent cohort, and is excluded
                             from every combination below

Four analyses, all from numbers already on disk:

  1. combination   Stouffer and Fisher across the two independent cohorts
  2. forest        Hedges' g with bootstrap CI for every measurement
  3. stability     leave-one-patient-out refits of the four-versus-four test
  4. concordance   per-patient transcript fraction against protein staining

Outputs (04_Revision_Analyses/10_R1.3_CrossCohort_Convergence/outputs/)
  pretx_sample_means.csv        per-sample CEACAM5/6 means, the Fig. 2K quantity
  crosscohort_combination.csv   one row per gene per combination method
  forest_effect_sizes.csv       one row per measurement
  loo_stability.csv             one row per dropped patient per gene
  rna_protein_concordance.csv   one row per marker
  convergence_report.txt
"""

from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD, NEW_ANALYSES  # noqa: E402
from shared.sample_ids import sample_id_map, to_study_ids  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
SWEEP = NEW_ANALYSES / "02_R1.3_TwoSided_Stats_Sweep" / "outputs" / "twosided_sweep.csv"
PAIRING = NEW_ANALYSES / "01_R1.3_Cohort_Pairing" / "outputs"
CEACAM = NEW_ANALYSES / "04_R1.5_CEACAM5_vs_CEACAM6" / "outputs"

RNG = np.random.default_rng(42)
N_BOOT = 10000
MIN_CELLS = 20

# The two rows of the sweep that come from patient sets with no overlap. The IHC
# rows are deliberately absent; see the module docstring.
INDEPENDENT = {
    "CEACAM6": ["CEACAM6 expression, pre-treatment epithelium",
                "CEACAM6 expression, PRJEB25780 deconvolved epithelium"],
    "CEACAM5": ["CEACAM5 expression, pre-treatment epithelium",
                "CEACAM5 expression, PRJEB25780 deconvolved epithelium"],
}
IHC_EXCLUDED = ("immunohistochemistry: same eight patients as the scRNA cohort, "
                "so not independent evidence")


# ----------------------------------------------------------------- statistics
def hedges_g(a, b):
    """Standardised mean difference with the small-sample correction."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n1, n2 = len(a), len(b)
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1))
                 / (n1 + n2 - 2))
    g = (a.mean() - b.mean()) / sp if sp > 0 else 0.0
    n = n1 + n2
    return float(g * (1 - 3 / (4 * n - 9))) if n > 3 else float(g)


def mann_whitney(nr, r):
    """Two-sided P and rank-biserial r, the same convention as the sweep."""
    u, p = stats.mannwhitneyu(nr, r, alternative="two-sided")
    return float(p), float(2 * u / (len(nr) * len(r)) - 1)


def stouffer(p_one, weights=None):
    """Combined one-tailed P; weights of None gives the unweighted form."""
    z = stats.norm.isf(np.asarray(p_one, float))
    w = np.ones_like(z) if weights is None else np.asarray(weights, float)
    return float(stats.norm.sf((w * z).sum() / np.sqrt((w ** 2).sum())))


def fisher(p_one):
    """Fisher's combined one-tailed P, reported as a check on Stouffer."""
    chi2 = -2 * np.log(np.asarray(p_one, float)).sum()
    return float(stats.chi2.sf(chi2, 2 * len(p_one)))


def spearman_ci(x, y, n=N_BOOT):
    """Percentile bootstrap CI for Spearman's rho over paired observations."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    idx = RNG.integers(0, len(x), size=(n, len(x)))
    out = [stats.spearmanr(x[i], y[i]).statistic for i in idx]
    out = np.asarray(out)
    out = out[np.isfinite(out)]
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


# ------------------------------------------------------------------- the data
def sample_means_from_h5ad(genes, phase="Pre", group_col="stomach_pre_grouping",
                           site="Stomach", min_cells=MIN_CELLS):
    """
    Sample-level mean expression, copied from run_twosided_sweep.py so that the
    leave-one-out refits below act on exactly the quantity Figure 2K plots. Both
    genes are pulled from one read; the object is five gigabytes.
    """
    ad = sc.read_h5ad(EPITHELIAL_H5AD)
    # Resolved before subsetting, so the crosswalk is checked against every
    # specimen in the object rather than only the pre-treatment stomach ones.
    ids = sample_id_map(ad.obs)
    if site is not None and "Sample site" in ad.obs:
        ad = ad[ad.obs["Sample site"] == site]
    ad = ad[ad.obs["Treatment phase"] == phase]
    ad = ad[ad.obs[group_col].isin(["Responsed", "No-response"])].copy()

    frames = []
    for gene in genes:
        src = ad.raw if (ad.raw is not None and gene in ad.raw.var_names) else ad
        idx = list(src.var_names).index(gene)
        x = src.X[:, idx]
        x = x.toarray().flatten() if hasattr(x, "toarray") else np.asarray(x).flatten()

        df = pd.DataFrame({"sample": ad.obs["sample"].astype(str).values,
                           "group": ad.obs[group_col].astype(str).values,
                           "value": x})
        keep = df.groupby("sample").size()
        df = df[df["sample"].isin(keep[keep >= min_cells].index)]
        out = df.groupby(["sample", "group"], observed=True)["value"].mean().reset_index()
        out["gene"] = gene
        out["group"] = out["group"].map({"No-response": "NR", "Responsed": "R"})
        # Specimens are named the way Supplementary Table 1 names them. The
        # relabel is in place, after the grouping, so no row moves and no value
        # changes.
        out["sample"] = to_study_ids(out["sample"], ids)
        frames.append(out)
    return pd.concat(frames, ignore_index=True)


def main():
    sweep = pd.read_csv(SWEEP).set_index("analysis")
    audit = pd.read_csv(PAIRING / "cohort_audit.csv")
    # Keyed by the study sample ID: every scRNA table read below is written with
    # the Supplementary Table 1 labels, and cohort_audit.csv carries that column
    # beside the study patient ID.
    crosswalk = dict(zip(audit["Sample"].astype(str), audit["Patient ID"]))

    # ------------------------------------------------------- 1. combination
    print("[1/4] combining the two independent cohorts ...")
    comb_rows = []
    for gene, labels in INDEPENDENT.items():
        sub = sweep.loc[labels]
        if not (sub.effect_size > 0).all():
            raise SystemExit(f"{gene}: the cohorts do not agree in direction, so "
                             "a one-tailed combination is not defensible")
        p_one = (sub.p_two_tailed / 2).values
        n = (sub.n_hi + sub.n_lo).values
        for method, p_combined in (
                ("Stouffer, unweighted", stouffer(p_one)),
                ("Stouffer, weighted by sqrt(n)", stouffer(p_one, np.sqrt(n))),
                ("Fisher", fisher(p_one))):
            comb_rows.append(dict(
                Gene=gene, Method=method,
                Cohorts="; ".join(f"{l.split(', ')[1]} (n = {int(a)} vs {int(b)}, "
                                  f"P = {p:.3f})"
                                  for l, a, b, p in zip(labels, sub.n_hi, sub.n_lo,
                                                        sub.p_two_tailed)),
                **{"P, combined two-sided": round(2 * p_combined, 5)},
                Excluded=IHC_EXCLUDED))
    combination = pd.DataFrame(comb_rows)
    combination.to_csv(OUT / "crosscohort_combination.csv", index=False)

    # ------------------------------------------------------------ 2. forest
    print("[2/4] assembling the effect-size forest ...")
    # The forest contrasts the CEACAM measurements that come from different
    # patient sets: transcript expression in the discovery epithelium, the same
    # in the PRJEB25780 cohort, and the immunohistochemistry on the discovery
    # patients themselves. The tumour-content-adjusted proportion added to the
    # sweep for Fig. S2D is not one of those: it is the same four-versus-four
    # discovery patients the expression rows already contribute, so it would
    # enter the forest as a second helping of the same evidence - and, under
    # the string rule below, would be filed as independent of the very cohort
    # it comes from. It belongs in Table S6, which is where it is.
    PROPORTION_ROWS = ("CEACAM5/6+ proportion, tumor-content adjusted",)
    ceacam = sweep[sweep.family == "CEACAM"]
    ceacam = ceacam[~ceacam.index.to_series().astype(str).str.startswith(PROPORTION_ROWS)
                    ] if ceacam.index.name == "analysis" else ceacam[
        ~ceacam["analysis"].astype(str).str.startswith(PROPORTION_ROWS)]
    forest = ceacam.reset_index()[
        ["analysis", "n_hi", "n_lo", "p_two_tailed", "effect_size",
         "hedges_g", "g_ci95_lo", "g_ci95_hi"]]
    forest.columns = ["Measurement", "n (NR)", "n (R)", "P, two-sided",
                      "Rank-biserial r", "Hedges g", "g 95% CI low",
                      "g 95% CI high"]
    forest["Independent of the scRNA cohort"] = [
        "no - same patients" if m.startswith("IHC") else
        ("discovery cohort" if "pre-treatment epithelium" in m else "yes")
        for m in forest["Measurement"]]
    forest.to_csv(OUT / "forest_effect_sizes.csv", index=False)

    # --------------------------------------------------------- 3. stability
    print("[3/4] leave-one-patient-out refits (reading the epithelial object) ...")
    means = sample_means_from_h5ad(("CEACAM6", "CEACAM5"))
    means["patient"] = means["sample"].map(crosswalk)
    if means["patient"].isna().any():
        raise SystemExit("unmapped scRNA samples: "
                         f"{sorted(means.loc[means.patient.isna(), 'sample'])}")
    means.to_csv(OUT / "pretx_sample_means.csv", index=False)

    loo_rows = []
    for gene, d in means.groupby("gene"):
        full_p, full_r = mann_whitney(d.loc[d.group == "NR", "value"],
                                      d.loc[d.group == "R", "value"])
        loo_rows.append(dict(Gene=gene, Dropped="none (as published)",
                             **{"n (NR)": int((d.group == "NR").sum()),
                                "n (R)": int((d.group == "R").sum()),
                                "P, two-sided": round(full_p, 4),
                                "Rank-biserial r": round(full_r, 3),
                                "Hedges g": round(hedges_g(
                                    d.loc[d.group == "NR", "value"],
                                    d.loc[d.group == "R", "value"]), 3)}))
        for patient in sorted(d.patient.unique()):
            k = d[d.patient != patient]
            nr, r = k.loc[k.group == "NR", "value"], k.loc[k.group == "R", "value"]
            p, rb = mann_whitney(nr, r)
            loo_rows.append(dict(Gene=gene, Dropped=patient,
                                 **{"n (NR)": len(nr), "n (R)": len(r),
                                    "P, two-sided": round(p, 4),
                                    "Rank-biserial r": round(rb, 3),
                                    "Hedges g": round(hedges_g(nr, r), 3)}))
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(OUT / "loo_stability.csv", index=False)

    # ------------------------------------------------------- 4. concordance
    print("[4/4] transcript against protein, per patient ...")
    frac = pd.read_csv(CEACAM / "ceacam_state_fractions.csv")
    frac["patient"] = frac["sample"].map(crosswalk)
    frac["CEACAM5"] = frac["CEACAM5+ only"] + frac["Double positive"]
    frac["CEACAM6"] = frac["CEACAM6+ only"] + frac["Double positive"]
    frac["Summed"] = frac["CEACAM5"] + frac["CEACAM6"]

    ihc = pd.read_csv(CEACAM / "ihc_per_marker_values.csv")
    ihc["patient"] = ihc["patient"].str.replace(r"^P0", "P", regex=True)
    ihc = ihc.rename(columns={"Summed (published)": "Summed"})

    merged = frac.merge(ihc, on="patient", suffixes=("_rna", "_protein"))
    if len(merged) != 8:
        raise SystemExit(f"expected eight patients in both modalities, got {len(merged)}")

    conc_rows = []
    for marker in ("CEACAM5", "CEACAM6", "Summed"):
        x = merged[f"{marker}_rna"].values
        y = merged[f"{marker}_protein"].values
        rho = stats.spearmanr(x, y)
        lo, hi = spearman_ci(x, y)
        conc_rows.append(dict(
            Marker=marker, n=len(x),
            **{"Spearman rho": round(float(rho.statistic), 3),
               "P": round(float(rho.pvalue), 4),
               "rho 95% CI low": round(lo, 3), "rho 95% CI high": round(hi, 3)}))
    concordance = pd.DataFrame(conc_rows)
    concordance.to_csv(OUT / "rna_protein_concordance.csv", index=False)
    merged[["patient", "group_rna", "CEACAM5_rna", "CEACAM6_rna",
            "CEACAM5_protein", "CEACAM6_protein"]].to_csv(
        OUT / "rna_protein_per_patient.csv", index=False)

    # ------------------------------------------------------------- report
    lines = [
        "CROSS-COHORT CONVERGENCE OF THE PRE-TREATMENT CEACAM5/6 RESULT",
        "Reviewer 1, point R1.3", "=" * 94, "",
        "A. WHICH COHORTS MAY BE COMBINED", "-" * 94,
        "   scRNA discovery      4 NR vs 4 R, pre-treatment gastric epithelium",
        "   PRJEB25780 (TIGER)   33 NR vs 12 R, deconvolved bulk epithelium",
        "   These share no patients, so their P values combine.", "",
        "   The immunohistochemistry is NOT combined with them. It was performed",
        "   on the same eight patients as the scRNA cohort, so it is an",
        "   orthogonal measurement of the same tumours rather than independent",
        "   evidence. Section D tests it as such.", "",
        "B. COMBINED ACROSS THE TWO INDEPENDENT COHORTS", "-" * 94]
    for gene in INDEPENDENT:
        lines.append(f"   {gene}")
        for _, row in combination[combination.Gene == gene].iterrows():
            lines.append(f"      {row.Method:<32} two-sided P = "
                         f"{row['P, combined two-sided']:.4f}")
        lines.append("")
    lines += ["C. LEAVE-ONE-PATIENT-OUT, FOUR VERSUS FOUR", "-" * 94]
    for gene, d in loo.groupby("Gene"):
        drops = d[d.Dropped != "none (as published)"]
        pub = d[d.Dropped == "none (as published)"].iloc[0]
        lines += [f"   {gene}   as published P = {pub['P, two-sided']:.4f}, "
                  f"g = {pub['Hedges g']:+.2f}",
                  f"      dropping one patient: P from {drops['P, two-sided'].min():.4f} "
                  f"to {drops['P, two-sided'].max():.4f}, "
                  f"g from {drops['Hedges g'].min():+.2f} to "
                  f"{drops['Hedges g'].max():+.2f}",
                  f"      direction preserved in {(drops['Hedges g'] > 0).sum()} "
                  f"of {len(drops)} refits", ""]
    lines += ["D. TRANSCRIPT AGAINST PROTEIN, SAME EIGHT PATIENTS", "-" * 94]
    for _, row in concordance.iterrows():
        lines.append(f"   {row.Marker:<10} rho = {row['Spearman rho']:+.3f} "
                     f"[{row['rho 95% CI low']:+.3f}, {row['rho 95% CI high']:+.3f}]"
                     f"   P = {row.P:.4f}")
    lines += ["",
              "   Provenance. The stained blocks are identified by hospital pathology",
              "   accession and the sequenced biopsies by internal sequencing ID. The",
              "   two are linked in the treating team's specimen worksheet, which",
              "   carries patient identifiers and is therefore not deposited: for all",
              "   eight patients the accession of the stained block is the accession",
              "   of the pre-treatment gastroscopic biopsy that was sequenced. The two",
              "   modalities therefore measure one specimen per patient, which is why",
              "   the immunohistochemistry is excluded from the combination in",
              "   section B and tested for concordance here instead.", ""]
    (OUT / "convergence_report.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"wrote six files to {OUT}")


if __name__ == "__main__":
    main()
