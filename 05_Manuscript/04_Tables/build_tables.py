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
  ST10 new: the sample-level (pseudobulk) sensitivity analysis of the NF-kB
       enrichment beside the adopted per-cell values (R1.8). The per-cell
       analysis is primary; the pseudobulk one is a cited sensitivity check

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
    # Supplementary Figure S2 panel D, whose printed P was one-sided.
    "Fig S2 (D)": "Fig S2D",
}

# Two of the sweep's panel labels are worn by more than one row, and the rows
# they cover do not all print in the same panel. Mapping on the panel label
# alone therefore gave three rows and two rows one printed pointer each, and
# three of those five pointers were wrong:
#
#   Fig 2 (Q)  covers the summed IHC comparison AND the two per-marker ones.
#              Figure 2N is 02_M/create_ihc_combined_boxplot.py, the COMBINED
#              CEACAM5+CEACAM6 boxplot, so it is right for the summed row only.
#              The per-marker comparison prints in Supplementary Figure S7
#              panel C - measured off the shipped S7 PDF, whose panel C is the
#              "DAB+ area (% of tissue)" row and prints P = 0.200 for CEACAM5
#              and P = 0.200 for CEACAM6, the two-tailed values these rows
#              carry; PROVENANCE.csv rows S7,B and S7,C put both under
#              _drivers/draw_ceacam5_vs_ceacam6.py, and the response letter's
#              own caption for S7B is the epithelial-percentage panel above it.
#
#   Fig 2 (K)  covers S-MP4 and S-MP5. Printed 2H and 2I are both drawn by
#              02_G/create_mp45_horizontal_boxplot.py, one script and two
#              panels: PROVENANCE.csv gives 2H mp4_horizontal_boxplot.svg and
#              2I mp5_horizontal_boxplot.svg, and reading the shipped Figure 2
#              at the printed_rect_mm of each confirms it - 2H prints "S-MP4"
#              with P = 0.083, 2I prints "S-MP5" with P = 0.188. So S-MP4 was
#              already right and S-MP5 was pointing at its neighbour.
#
# Rows listed here take their printed panel from the comparison rather than
# from the sweep's panel label. check_panel_labels() below refuses to build if
# a shared label is left to resolve by itself again.
PANEL_BY_ANALYSIS = {
    "IHC staining, CEACAM5 + CEACAM6 summed": "Fig 2N",
    "IHC staining, CEACAM5 only": "Fig S7C (CEACAM5)",
    "IHC staining, CEACAM6 only": "Fig S7C (CEACAM6)",
    "S-MP4 score, pre-treatment R vs all other groups": "Fig 2H",
    "S-MP5 score, pre-treatment R vs all other groups": "Fig 2I",
}


def check_panel_labels(sweep):
    """Refuse to build while a shared sweep panel label is resolved by itself.

    The defect this guards against is not a wrong entry in a table; it is a
    label that two different comparisons share, quietly collapsing to whichever
    printed panel the map happens to name. Nothing said so - the map was
    complete, every key resolved, and the table built. So the guard is on the
    shape of the input rather than on the values: if a sweep panel label is
    worn by more than one row, every one of those rows must say for itself
    which printed panel it belongs to.
    """
    unknown = set(sweep["panel"]) - set(SUBMITTED_PANEL)
    if unknown:
        raise SystemExit(f"sweep panel labels with no submitted equivalent: {unknown}")
    shared = sweep.groupby("panel")["analysis"].apply(list)
    loose = sorted(
        f"{panel!r} is shared by {len(rows)} comparisons; "
        f"{a!r} is not in PANEL_BY_ANALYSIS"
        for panel, rows in shared.items() if len(rows) > 1
        for a in rows if a not in PANEL_BY_ANALYSIS)
    if loose:
        raise SystemExit("\n".join(
            ["a sweep panel label covers comparisons that may print in "
             "different panels, and at least one of them has not been told "
             "which is its own:", *[f"  - {m}" for m in loose],
             "", "Add it to PANEL_BY_ANALYSIS, after establishing the printed "
             "panel from PROVENANCE.csv and the shipped figure."]))
    orphan = sorted(set(PANEL_BY_ANALYSIS) - set(sweep["analysis"]))
    if orphan:
        raise SystemExit(f"PANEL_BY_ANALYSIS names comparisons the sweep does "
                         f"not have: {orphan}")


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

    # ----------------------------------------- ST2 loses two retired strings
    # Both are inherited verbatim from the frozen source and describe things
    # that no longer print anywhere:
    #   "C3_Mac" is the pre-rename cluster name, in the GSE239676 row, which
    #   describes Figure 4H - and Figure 4H prints "MoMac_IL1B Proportion".
    #   "NF-kB" is the ASCII spelling, in the TCGA-STAD row. The manuscript
    #   uses the Greek kappa throughout and ASCII zero times.
    st2_path = OUT / "ST2_external_cohort_details.csv"
    st2_text = st2_path.read_text(encoding="utf-8")
    for old, new in (("C3_Mac", "MoMac_IL1B"), ("NF-kB", "NF-\u03baB")):
        if old not in st2_text:
            raise SystemExit(
                f"build_tables(): {old!r} is no longer in ST2. It was there "
                f"when this fix was written; if the frozen source changed, "
                f"re-check this fix.")
        st2_text = st2_text.replace(old, new)
    st2_path.write_text(st2_text, encoding="utf-8")

    # ------------------------------------ ST2 "Usage in paper" matches the paper
    # Two inherited entries described analyses the manuscript does not contain:
    # PRJEB25780 was never used for survival (the only survival analysis is
    # TCGA-STAD, Fig. S2G), and GSE239676 was never used for a response
    # association (it is tumour versus adjacent normal, Fig. 4H). Each entry now
    # names the panels the cohort actually appears in.
    st2 = pd.read_csv(OUT / "ST2_external_cohort_details.csv")
    USAGE = {
        "PRJEB25780": (
            "BayesPrism deconvolution of epithelial expression; CEACAM5 and "
            "CEACAM6 in responders versus non-responders and the cross-cohort "
            "combination (Fig. 2L, Table S9); epithelial metaprogram scores by "
            "response; purity-adjusted association of epithelial CEACAM5/6 "
            "with immune infiltration"),
        "GSE239676": (
            "Independent validation of IL-1\u03b2+ MoMac (MoMac_IL1B) enrichment in "
            "primary tumor versus adjacent normal tissue (Fig. 4H)"),
    }
    for ds, text in USAGE.items():
        hit = st2["Dataset ID"] == ds
        if hit.sum() != 1:
            raise SystemExit(f"build_tables(): expected one ST2 row for {ds}, "
                             f"found {int(hit.sum())}")
        st2.loc[hit, "Usage in paper"] = text

    # ------------------------------------- ST2's GSE251950 patient count is 9
    # The frozen source says N patients = 10, N samples = 10. The project's own
    # figure says otherwise, and it is the figure the paper prints: the six
    # S4_Spatial_Validation panel scripts each carry the same SAMPLE_INFO map,
    #
    #   sample_01..sample_10 -> GC1 GC2 GC3 GC4 GC5 GC6 GC6-PM GC7 GC8 GC9
    #
    # which is ten sections from NINE patients, GC6 contributing both a primary
    # and its paired peritoneal metastasis. The shipped page agrees: the text
    # layer of 03_Supplementary_Figures/S4_Spatial_Validation.pdf carries those
    # ten labels, six times each, across nine distinct patient names.
    #
    # So N samples stays 10 - ten sections is right, and it is what the Results'
    # "a pattern seen in all 10 samples" and the Methods' "10 Visium samples"
    # are counting. Only the patient count moves. Author's ruling 2026-09-10:
    # this one field, and no accompanying sentence in the Methods or the S4
    # legend.
    #
    # Corrected here rather than in the CSV because 01_Reviewer_Materials/
    # tables_submitted/ is the frozen submission and stays as submitted;
    # this function is the one owner of what ST2 ships as.
    hit = st2["Dataset ID"] == "GSE251950"
    if hit.sum() != 1:
        raise SystemExit(f"build_tables(): expected one ST2 row for GSE251950, "
                         f"found {int(hit.sum())}")
    was = st2.loc[hit, "N patients"].iloc[0]
    if int(was) not in (9, 10):
        raise SystemExit(
            f"build_tables(): ST2's GSE251950 N patients is {was!r}. It was 10 "
            f"when this fix was written and the correction takes it to 9; if "
            f"the frozen source changed, re-check this against the "
            f"S4_Spatial_Validation SAMPLE_INFO map before overwriting it.")
    if int(st2.loc[hit, "N samples"].iloc[0]) != 10:
        raise SystemExit(
            "build_tables(): ST2's GSE251950 N samples is not 10. Ten sections "
            "is the count the manuscript and Fig. S4 both make; a change there "
            "is not this fix's to absorb.")
    st2.loc[hit, "N patients"] = 9
    st2.to_csv(OUT / "ST2_external_cohort_details.csv", index=False)

    # What the correction deliberately does NOT reach, recorded 2026-09-10 so
    # that a later reader finds it stated rather than discovers it.
    #
    # ST6's three GSE251950 rows - Fig 3H, Fig 3I, Fig 3J - keep n = 10 in both
    # group-size columns, and the paired Wilcoxon signed-rank statistics behind
    # them are unchanged. Ten paired sections is what the test was run on and
    # what it is right to report, and two of those ten sections (GC6 and
    # GC6-PM) come from one patient. The author ruled that this is a stated
    # fact and not outstanding work: the patient count is corrected in
    # Supplementary Table 2 and nowhere else, no disclosure sentence is added
    # to the Methods, the S4 legend or the Results, and no statistic is re-run.
    # Do not "finish" this by re-running the spatial tests on nine patients.

    # ------------------------------------------------ ST3 names and Tex lists
    # "(C3_Mac)" is the pre-rename cluster code; the manuscript's name for the
    # state is IL-1beta+ MoMac. Two exhaustion signatures are in use and the
    # table listed one: the 19-gene state score (Fig. 5K,
    # 05_O/create_tex_nfkb_scatter_cd8.py STATE_GENES, scanpy score_genes) and
    # the 18-gene ssGSEA signature of the TCGA panel (Fig. 3E,
    # 03_D/create_tcga_ceacam_scatter.py TEX_GENES), whose legend says
    # "18-gene signature". Both are now listed with the panel each serves.
    st3 = pd.read_csv(OUT / "ST3_gene_signatures.csv")
    if not st3["Signature"].str.contains("(C3_Mac)", regex=False).any():
        raise SystemExit("build_tables(): '(C3_Mac)' is no longer in ST3; re-check")
    st3["Signature"] = st3["Signature"].str.replace(
        "(C3_Mac)", "(IL-1\u03b2+ MoMac)", regex=False)
    tex = st3["Signature"] == "T cell exhaustion (Tex)"
    if tex.sum() != 1:
        raise SystemExit("build_tables(): expected one Tex row in ST3")
    if int(st3.loc[tex, "N_genes"].iloc[0]) != 19:
        raise SystemExit("build_tables(): the ST3 Tex row no longer has 19 genes")
    st3.loc[tex, "Signature"] = "T cell exhaustion (Tex), 19-gene state score (Fig. 5K)"
    tex18 = pd.DataFrame([{
        "Signature": "T cell exhaustion (Tex), 18-gene ssGSEA signature (Fig. 3E)",
        "N_genes": 18,
        "Genes": ("PDCD1; HAVCR2; LAG3; TIGIT; CTLA4; TOX; ENTPD1; CXCL13; LAYN; "
                  "CD38; BATF; IRF4; PRDM1; TOX2; ITGAE; NR4A1; NR4A2; NR4A3"),
        "Reference": "Wherry & Kurachi, Nat Rev Immunol 2015; this study",
    }])
    # The composite Tex score of Fig. S3D is a third list: the ten-gene
    # TEX_GENES of create_S3_D_cd8_tex_score_umap.py (scanpy score_genes).
    tex10 = pd.DataFrame([{
        "Signature": "T cell exhaustion (Tex), 10-gene composite score (Fig. S3D)",
        "N_genes": 10,
        "Genes": "PDCD1; HAVCR2; LAG3; TIGIT; CTLA4; TOX; ENTPD1; LAYN; CXCL13; BATF",
        "Reference": "Wherry & Kurachi, Nat Rev Immunol 2015; this study",
    }])
    i = int(st3.index[tex][0]) + 1
    st3 = pd.concat([st3.iloc[:i], tex18, tex10, st3.iloc[i:]], ignore_index=True)
    # The IL-1beta+ state signature: the top 50 Wilcoxon markers of
    # the state against the other monocyte/macrophage states, written by
    # 08_R2.1_PreTx_Inflammatory/scripts (il1b_signature_genes.csv). Distinct
    # from the 15-gene external-validation signature described in the Methods.
    il1b = pd.read_csv(NEW_ANALYSES / "08_R2.1_PreTx_Inflammatory" / "outputs"
                       / "il1b_signature_genes.csv")["gene"].tolist()
    if len(il1b) != 50:
        raise SystemExit(f"build_tables(): expected 50 IL-1beta state markers, "
                         f"found {len(il1b)}")
    st3 = pd.concat([st3, pd.DataFrame([{
        # The Methods sentence that defined this signature described a
        # withdrawn panel and went with it, so the definition lives here, in
        # the table that carries the genes.
        "Signature": "IL-1\u03b2+ MoMac state signature: the 50 strongest "
                     "Wilcoxon markers of the state against the other "
                     "monocyte/macrophage states",
        "N_genes": 50, "Genes": "; ".join(il1b), "Reference": "This study",
    }])], ignore_index=True)
    st3.to_csv(OUT / "ST3_gene_signatures.csv", index=False)

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
        # paired_any_site is one row PER PATIENT in pairing_summary.csv, and the
        # merge is on Patient ID, so the flag is broadcast onto every specimen
        # row of that patient. Under the old name a reader counted three "Yes"
        # rows and read three pairs; they are three specimens of ONE patient
        # (P01), and one of them is a liver sample. The column name now says
        # so: 1 of 35 patients was sampled at both timepoints.
        "paired_any_site": "Patient sampled at both timepoints (any site)",
    })
    col = "Patient sampled at both timepoints (any site)"
    st1[col] = st1[col].map({True: "Yes", False: "No"})
    st1.to_csv(OUT / "ST1_patient_sample_characteristics.csv", index=False)

    # --------------------------------------------------------- ST6, new
    sweep = pd.read_csv(NEW_ANALYSES / "02_R1.3_TwoSided_Stats_Sweep" / "outputs"
                        / "twosided_sweep.csv")
    st6 = sweep[[
        "panel", "analysis", "family", "test", "group_hi", "n_hi", "group_lo",
        "n_lo", "mean_hi", "mean_lo", "p_one_tailed", "p_two_tailed",
        "p_two_tailed_BH", "effect_size_name",
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
        # 2026-09-11, author's ruling: the "Smallest two-tailed P attainable at
        # these n" column is dropped. Its legend gloss was removed on
        # 2026-09-10 for the same reason the Discussion's floor clause was - a
        # limitation is acknowledged, its floor price is not quoted - which left
        # the column shipping unexplained. The ruling is that the column goes
        # rather than the gloss coming back. The quantity itself is not lost: it
        # is p_two_tailed_floor in
        # 02_R1.3_TwoSided_Stats_Sweep/outputs/twosided_sweep.csv, which is what
        # verify_numbers.py's "n=4 vs 4 two-sided floor" check reads, and the
        # response letter states the value in prose without citing this table.
        # Every other column and every value is unchanged.
        "effect_size_name": "Effect size measure", "effect_size": "Effect size",
        "hedges_g": "Hedges g", "g_ci95_lo": "Hedges g 95% CI lower",
        "g_ci95_hi": "Hedges g 95% CI upper",
        "diff_of_means": "Difference of means",
        "ci95_lo": "Difference 95% CI lower", "ci95_hi": "Difference 95% CI upper",
        "verdict": "Robustness under two-sided testing",
    })
    check_panel_labels(sweep)
    st6["Figure panel"] = [
        PANEL_BY_ANALYSIS.get(a, SUBMITTED_PANEL[p])
        for p, a in zip(st6["Figure panel"], st6["Comparison"])]
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

    # --------------------------------------------------------- ST10, new
    # The per-cell Welch analysis of the NF-kB enrichment is primary; the
    # sample-level pseudobulk analysis (15_Pseudobulk_Sample_Level: summed UMIs
    # per sample, limma-voom and DESeq2, the same prerank GSEA) ships as a
    # table and is cited as a sensitivity check rather than replacing it. The
    # ttest_percell_* columns of the source are sound13, the adopted table, so
    # the per-cell values printed here are those of Fig. S9C and the Results.
    pb = pd.read_csv(NEW_ANALYSES / "15_Pseudobulk_Sample_Level" / "outputs"
                     / "nfkb_comparison.csv")
    if not bool(pb["testable"].all()) or len(pb) != 26:
        raise SystemExit("build_tables(): expected 26 testable contrasts in "
                         "nfkb_comparison.csv")
    labels = {"B_cells": "B cells", "DC_cells": "DC",
              "Endothelial_cells": "Endothelial", "Epithelial": "Epithelial",
              "Fibroblast": "Fibroblast", "Mast_cells": "Mast", "MoMac": "MoMac",
              "Neutrophils": "Neutrophils", "NK_cells": "NK",
              "Pericyte": "Pericyte", "Plasma_cells": "Plasma",
              "TCD4_cells": "CD4+ T", "TCD8_cells": "CD8+ T"}
    st10 = pd.DataFrame({
        "Cell type": pb["cell_type"].map(labels),
        "Timepoint": pb["phase"].map({"pre": "Pre-treatment",
                                      "post": "Post-treatment"}),
        "n samples": pb["n_samples"], "n responders": pb["n_R"],
        "n non-responders": pb["n_NR"], "n cells": pb["n_cells"],
        "Median cells per sample": pb["median_cells_per_sample"],
        "Samples dropped (fewer than 10 cells)": pb["samples_dropped"],
        "limma-voom NES": pb["limma_nes"].round(3),
        "limma-voom nominal P": pb["limma_p"].round(4),
        "limma-voom FDR q": pb["limma_q"].round(4),
        "limma-voom rank": pb["limma_rank"],
        "limma-voom Hallmark sets tested": pb["limma_nsets"],
        "DESeq2 NES": pb["deseq2_nes"].round(3),
        "DESeq2 nominal P": pb["deseq2_p"].round(4),
        "DESeq2 FDR q": pb["deseq2_q"].round(4),
        "DESeq2 rank": pb["deseq2_rank"],
        "DESeq2 Hallmark sets tested": pb["deseq2_nsets"],
        "Per-cell NES (primary analysis; Fig. S9C)": pb["ttest_percell_nes"].round(3),
        "Per-cell nominal P": pb["ttest_percell_p"].round(4),
        "Per-cell FDR q": pb["ttest_percell_q"].round(4),
        "Per-cell rank": pb["ttest_percell_rank"],
        "Per-cell Hallmark sets tested": pb["ttest_percell_nsets"],
    })
    if st10["Cell type"].isna().any():
        raise SystemExit("build_tables(): unlabelled cell type in ST10")
    st10.to_csv(OUT / "ST10_nfkb_pseudobulk_sensitivity.csv", index=False)

    print("Supplementary tables written to", OUT)
    for f in sorted(OUT.glob("ST*.csv")):
        d = pd.read_csv(f)
        print(f"   {f.name:<48} {d.shape[0]:>4} rows x {d.shape[1]:>2} cols")


if __name__ == "__main__":
    main()
