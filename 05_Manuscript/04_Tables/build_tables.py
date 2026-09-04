"""
Supplementary tables for the revision.

  ST1  expanded with the per-patient timepoint, RECIST trajectory and regimen
       that Reviewer 1 asked for (R1.3a, R1.3b)
  ST6  new: every directional comparison reported with both tails, an effect
       size, a bootstrap CI and the BH-adjusted P (R1.3c)
  ST7  new: CEACAM5-only / CEACAM6-only / double-positive fractions and the
       per-marker IHC values (R1.5)
  ST8  new: every pre-treatment comparison repeated with the two specimens the
       cohort audit flagged reclassified or excluded (R1.3b)
  ST9  new: effect sizes, the combination of the two cohorts that share no
       patients, leave-one-patient-out refits and transcript-protein
       concordance (R1.3)

ST2-ST5 are carried over unchanged from the submission.
"""

from pathlib import Path
import shutil
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from paths import NEW_ANALYSES, REVIEWER_MATERIALS  # noqa: E402

OUT = Path(__file__).parent
SRC = REVIEWER_MATERIALS / "tables_submitted"


# The sweep labels panels by the figure code as it stands today, which was
# re-lettered after submission: it calls Figure 2 A-Q with split N1/N2 and O1/O2
# panels. The paper prints Figure 2 as A-N. Passing the sweep's labels through
# would send a reader checking a P value to the wrong panel, so they are
# translated here, in one place.
SUBMITTED_PANEL = {
    "Fig 2 (N1)": "Fig 2K left (CEACAM6)",
    "Fig 2 (N2)": "Fig 2K right (CEACAM5)",
    "Fig 2 (O1)": "Fig 2L left (CEACAM6)",
    "Fig 2 (O2)": "Fig 2L right (CEACAM5)",
    "Fig 2 (Q)": "Fig 2N",
    "Fig 2 (K)": "Fig 2H",
    "Fig 5 (J)": "Fig 5J monocytes/macrophages",
    "Fig 5 (K)": "Fig 5J epithelial cells",
    "Fig 5 (L)": "Fig 5J fibroblasts",
    "Fig 5 (M)": "Fig 5J dendritic cells",
    "Fig 3 (H)": "Fig 3H",
    "Fig 3 (I)": "Fig 3I",
    "Fig 3 (J)": "Fig 3J",
    "Fig 5 (regulon-BACH1)": "Fig 5D",
    "Fig 5 (regulon-NFKB1)": "Fig 5E",
    "Fig 5 (IL6-CD4)": "Fig 5L",
}


def main():
    # ---------------------------------------------- carry over ST2-ST5 as-is
    for f in sorted(SRC.glob("ST*.csv")):
        if f.name.startswith("ST1_"):
            continue
        # The frozen source is read-only and copy2 preserves the mode, so a
        # second run would fail trying to overwrite its own output.
        dst = OUT / f.name
        if dst.exists():
            dst.chmod(0o644)
        shutil.copy2(f, dst)
        dst.chmod(0o644)

    # ------------------------------------------- ST2 gains the sixth cohort
    # Figure S11A is titled "GSE246011 replication" and that accession appeared
    # nowhere in the manuscript, the Methods or this table, which listed the
    # other five external cohorts and GSE251950 in its place. A dataset shown in
    # a figure and defined nowhere is the kind of thing a reviewer finds first.
    #
    # The usage wording is written from the numbers, not from the analysis
    # script's own summary text: gse246011_replication.csv gives a positive
    # CEACAM coefficient both unadjusted (+0.0505) and adjusted for epithelial
    # content (+0.0333), while the script prints an unconditional caveat block
    # written for a negative coefficient. The three definitional differences it
    # lists are real and are carried into the note; the sign claim in it is not.
    st2 = pd.read_csv(OUT / "ST2_external_cohort_details.csv")
    if "GSE246011" not in set(st2["Dataset ID"]):
        st2.loc[len(st2)] = {
            "Dataset ID": "GSE246011",
            "Year": 2024,
            "N patients": 4,
            "N samples": 4,
            "Data type": "Spatial transcriptomics (Visium)",
            "Platform": "10x Genomics Visium",
            "Usage in paper": (
                "Supporting spatial analysis of the CEACAM-immune distance "
                "association in an independent gastric cohort (Fig. S11A)"),
            "Note": (
                "Not a like-for-like replication of the primary spatial "
                "analysis: immune-rich spots are defined by the top PTPRC "
                "quartile rather than a deconvolved immune fraction, distance "
                "is to the single nearest immune spot rather than the mean of "
                "the five nearest, and CEACAM exposure is the raw per-spot sum "
                "with epithelial content as a covariate."),
        }
        st2.to_csv(OUT / "ST2_external_cohort_details.csv", index=False)

    # ------------------------------------------------------- ST1, expanded
    st1 = pd.read_csv(SRC / "ST1_patient_sample_characteristics.csv")
    audit = pd.read_csv(NEW_ANALYSES / "01_R1.3_Cohort_Pairing" / "outputs"
                        / "cohort_audit.csv")
    pairing = pd.read_csv(NEW_ANALYSES / "01_R1.3_Cohort_Pairing" / "outputs"
                          / "pairing_summary.csv")

    st1 = st1.merge(
        audit[["Sample", "regimen", "recist_raw", "recist_best", "recist_change"]],
        on="Sample", how="left")
    st1 = st1.merge(
        pairing[["Patient ID", "paired_any_site"]], on="Patient ID", how="left")
    st1 = st1.rename(columns={
        "regimen": "Treatment regimen and cycles",
        "recist_raw": "RECIST 1.1 trajectory",
        "recist_best": "Best response",
        "recist_change": "Change after best response",
        "paired_any_site": "Sampled at both timepoints",
    })
    st1["Sampled at both timepoints"] = st1["Sampled at both timepoints"].map(
        {True: "Yes", False: "No"})
    st1.to_csv(OUT / "ST1_patient_sample_characteristics.csv", index=False)

    # --------------------------------------------------------- ST6, new
    sweep = pd.read_csv(NEW_ANALYSES / "02_R1.3_TwoSided_Stats_Sweep" / "outputs"
                        / "twosided_sweep.csv")
    st6 = sweep[[
        "panel", "analysis", "family", "test", "group_hi", "n_hi", "group_lo",
        "n_lo", "mean_hi", "mean_lo", "p_one_tailed", "p_two_tailed",
        "p_two_tailed_BH", "p_two_tailed_floor", "effect_size_name",
        "effect_size", "hedges_g", "g_ci95_lo", "g_ci95_hi", "diff_of_means",
        "ci95_lo", "ci95_hi", "verdict",
    ]].rename(columns={
        "panel": "Figure panel",
        "analysis": "Comparison",
        "family": "Hypothesis family",
        "test": "Test",
        "group_hi": "Group 1", "n_hi": "n (group 1)",
        "group_lo": "Group 2", "n_lo": "n (group 2)",
        "mean_hi": "Mean (group 1)", "mean_lo": "Mean (group 2)",
        "p_one_tailed": "P, one-tailed (as submitted)",
        "p_two_tailed": "P, two-tailed (revised)",
        "p_two_tailed_BH": "P, two-tailed, BH within family",
        "p_two_tailed_floor": "Smallest two-tailed P attainable at these n",
        "effect_size_name": "Effect size measure", "effect_size": "Effect size",
        "hedges_g": "Hedges g", "g_ci95_lo": "Hedges g 95% CI lower",
        "g_ci95_hi": "Hedges g 95% CI upper",
        "diff_of_means": "Difference of means",
        "ci95_lo": "Difference 95% CI lower", "ci95_hi": "Difference 95% CI upper",
        "verdict": "Robustness under two-sided testing",
    })
    unknown = set(st6["Figure panel"]) - set(SUBMITTED_PANEL)
    if unknown:
        raise SystemExit(f"sweep panel labels with no submitted equivalent: {unknown}")
    st6["Figure panel"] = st6["Figure panel"].map(SUBMITTED_PANEL)
    st6.to_csv(OUT / "ST6_two_sided_sensitivity.csv", index=False)

    # --------------------------------------------------------- ST7, new
    states = pd.read_csv(NEW_ANALYSES / "04_R1.5_CEACAM5_vs_CEACAM6" / "outputs"
                         / "ceacam_state_tests.csv")
    ihc = pd.read_csv(NEW_ANALYSES / "04_R1.5_CEACAM5_vs_CEACAM6" / "outputs"
                      / "ihc_per_marker_tests.csv")
    st7 = pd.concat([states, ihc], ignore_index=True).rename(columns={
        "measure": "Measure", "level": "Data level",
        "n_NR": "n non-responders", "n_R": "n responders",
        "mean_NR": "Mean, non-responders", "mean_R": "Mean, responders",
        "p_two_tailed": "P, two-tailed",
        "rank_biserial_r": "Rank-biserial correlation",
    })
    st7.to_csv(OUT / "ST7_ceacam5_vs_ceacam6.csv", index=False)

    # --------------------------------------------------------- ST8, new
    # How far the pre-treatment comparisons depend on the response label of the
    # two specimens the cohort audit could not reconcile with the SD/PD rule.
    st8 = pd.read_csv(NEW_ANALYSES / "01_R1.3_Cohort_Pairing" / "outputs"
                      / "response_label_sensitivity.csv")
    st8.to_csv(OUT / "ST8_response_label_sensitivity.csv", index=False)

    # --------------------------------------------------------- ST9, new
    # The four convergence analyses have different natural shapes, so they are
    # stacked into one schema: what was analysed, on what, and the statistic.
    conv = NEW_ANALYSES / "10_R1.3_CrossCohort_Convergence" / "outputs"
    forest = pd.read_csv(conv / "forest_effect_sizes.csv")
    comb = pd.read_csv(conv / "crosscohort_combination.csv")
    loo = pd.read_csv(conv / "loo_stability.csv")
    conc = pd.read_csv(conv / "rna_protein_concordance.csv")

    def block(analysis, measurement, comparison, n_nr, n_r, stat, value,
              lo, hi, p):
        return pd.DataFrame({
            "Analysis": analysis, "Measurement": measurement,
            "Comparison": comparison, "n non-responders": n_nr,
            "n responders": n_r, "Statistic": stat, "Value": value,
            "95% CI low": lo, "95% CI high": hi, "P, two-sided": p})

    st9 = pd.concat([
        block("Effect size", forest["Measurement"],
              forest["Independent of the scRNA cohort"].map(
                  {"discovery cohort": "discovery cohort",
                   "yes": "independent of the discovery cohort",
                   "no - same patients": "same patients as the discovery cohort"}),
              forest["n (NR)"], forest["n (R)"], "Hedges g",
              forest["Hedges g"].round(3), forest["g 95% CI low"].round(3),
              forest["g 95% CI high"].round(3), forest["P, two-sided"].round(4)),
        block("Combined across cohorts sharing no patients", comb["Gene"],
              comb["Cohorts"], None, None, comb["Method"], None, None, None,
              comb["P, combined two-sided"]),
        block("Leave-one-patient-out", loo["Gene"],
              loo["Dropped"].replace({"none (as published)": "none (as published)"}),
              loo["n (NR)"], loo["n (R)"], "Hedges g", loo["Hedges g"],
              None, None, loo["P, two-sided"]),
        block("Transcript fraction against protein staining", conc["Marker"],
              "same eight patients, both modalities", None, None,
              "Spearman rho", conc["Spearman rho"], conc["rho 95% CI low"],
              conc["rho 95% CI high"], conc["P"]),
    ], ignore_index=True)
    st9.to_csv(OUT / "ST9_crosscohort_convergence.csv", index=False)

    print("Supplementary tables written to", OUT)
    for f in sorted(OUT.glob("ST*.csv")):
        d = pd.read_csv(f)
        print(f"   {f.name:<48} {d.shape[0]:>4} rows x {d.shape[1]:>2} cols")


if __name__ == "__main__":
    main()
