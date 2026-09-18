"""
Supplementary tables for the revision.

  ST1  expanded with the per-patient timepoint, adjudicated RECIST 1.1 response
       that Reviewer 1 asked for (R1.3a, R1.3b)
  ST6  new: every directional comparison reported two-sided, with an effect
       size, a bootstrap CI and the BH-adjusted P (R1.3c)
ST2-ST5 are carried over unchanged from the submission. The former ST7 and
ST8 audit outputs remain reproducible from 04_Revision_Analyses but are not shipped
as supplementary tables.
"""

from pathlib import Path
import re
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
    # Supplementary Figure S2 panel D as submitted, whose printed P was
    # one-sided; S3 D since the renumbering of 2026-09-16 (see RENUMBERED_S).
    "Fig S2 (D)": "Fig S3D",
}

# SUPPLEMENTARY FIGURES RENUMBERED ON 2026-09-16 (the author's fifth reading):
# the submitted S1 was split into S1 (QC) and S2 (cell-type annotation), and
# the former S2-S9 are S3-S10. The analysis outputs the tables are built from
# name panels in their comparison labels ("... (Fig. S2D)") and are frozen
# (freeze_baseline.py), so the translation is made here, in the one place
# that writes the shipped tables, and nowhere else. renumber_supp() refuses a
# label that names a supplementary panel this map does not know, so a new
# analysis label cannot ship under a stale number.
RENUMBERED_S = {"S2D": "S3D", "S3D": "S4D", "S7A": "S8A", "S7B": "S8B",
                "S7C": "S8C", "S9C": "S10C"}
_SUPP_PANEL = re.compile(r"(Fig\.? ?)S(\d+)([A-Z])\b")


def renumber_supp(text):
    """A comparison label with its supplementary panel under today's number."""
    if not isinstance(text, str):
        return text

    def one(m):
        key = f"S{m.group(2)}{m.group(3)}"
        if key not in RENUMBERED_S:
            raise KeyError(f"{key!r} in {text!r}: not in RENUMBERED_S - "
                           "decide its post-2026-09-16 number before shipping")
        return m.group(1) + RENUMBERED_S[key]
    return _SUPP_PANEL.sub(one, text)

# Two of the sweep's panel labels are worn by more than one row, and the rows
# they cover do not all print in the same panel. Mapping on the panel label
# alone therefore gave three rows and two rows one printed pointer each, and
# three of those five pointers were wrong:
#
#   Fig 2 (Q)  covers the summed IHC comparison AND the two per-marker ones.
#              Figure 2N is 02_M/create_ihc_combined_boxplot.py, the COMBINED
#              CEACAM5+CEACAM6 boxplot, so it is right for the summed row only.
#              The per-marker comparison prints in Supplementary Figure S7
#              panel C (S8 C since 2026-09-16) - measured off the shipped S7 PDF, whose panel C is the
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
    "IHC staining, CEACAM5 only": "Fig S8C (CEACAM5)",
    "IHC staining, CEACAM6 only": "Fig S8C (CEACAM6)",
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


# --------------------------------------------- ST1's adjudicated response
# Supplementary Table 1 reports the fixed RECIST 1.1 assessment in one column.
# Calls are keyed by patient so every specimen from that patient receives the
# same assessment. Patients without a definitive call remain blank.
RECIST_RESPONSE = {
    "P1": "SD", "P2": "PD", "P3": "PR", "P4": "PR", "P5": "PR",
    "P6": "PD", "P7": "", "P10": "SD", "P11": "PR", "P12": "PR",
    "P13": "PR", "P14": "PD", "P16": "PD", "P21": "PR", "P22": "PR",
    "P23": "SD", "P24": "PD", "P25": "PD", "P26": "PD", "P33": "",
    "P34": "PD",
}


def _check_recist_response(st1):
    """Gate the adjudicated call against ST1's own R/NR Grouping.

    RECIST 1.1 makes a responder CR or PR. Every response-labelled row must
    therefore carry PR against R and PD against NR; SD is the one call the
    grouping cannot be predicted from, because it was decided on the narrative
    tendency, and blank rows have no call to check.

    This is a gate rather than an audit because the two quantities live in
    different columns of the same published table and are read side by side: if
    they ever disagree again, the table must not be written. The mutation that
    shows the gate can fail is in the module's __main__ block.
    """
    bad = [v for v in RECIST_RESPONSE.values() if v not in ("PR", "SD", "PD", "")]
    if bad:
        raise SystemExit(f"build_tables(): RECIST_RESPONSE holds calls that are "
                         f"not a RECIST 1.1 response: {sorted(set(bad))}")

    missing = sorted(set(st1.loc[st1["R/NR Grouping"].notna(), "Patient ID"])
                     - set(RECIST_RESPONSE))
    if missing:
        raise SystemExit(f"build_tables(): response-labelled patients with no "
                         f"adjudicated RECIST call: {missing}")

    expected = {"PR": "R", "PD": "NR"}
    labelled = st1[st1["R/NR Grouping"].notna()]
    clash = [f"{r['Sample']}: RECIST {r['RECIST 1.1 response']!r} against "
             f"R/NR Grouping {r['R/NR Grouping']!r}"
             for _, r in labelled.iterrows()
             if r["RECIST 1.1 response"] in expected
             and expected[r["RECIST 1.1 response"]] != r["R/NR Grouping"]]
    if clash:
        raise SystemExit("\n".join(
            ["build_tables(): the adjudicated RECIST response contradicts ST1's "
             "own R/NR Grouping:", *[f"  - {c}" for c in clash]]))

    # The composition the Methods and the response letter both quote. It is
    # asserted here because this is the table those sentences cite.
    split = (labelled[labelled["Anatomical site"] == "Stomach"]
             .groupby(["Treatment phase", "R/NR Grouping"]).size().to_dict())
    if split != {("Post", "NR"): 6, ("Post", "R"): 5,
                 ("Pre", "NR"): 4, ("Pre", "R"): 4}:
        raise SystemExit(f"build_tables(): the response-labelled gastric split is "
                         f"no longer 4 R / 4 NR pre and 5 R / 6 NR post: {split}")


def main():
    # Guard against retired tables surviving from an older build. Neither table
    # belongs to the fixed-label ST1-ST8 layout.
    for obsolete in (
        "ST7_ceacam5_vs_ceacam6.csv",
        "ST8_nfkb_pseudobulk_sensitivity.csv",
        "ST9_crosscohort_convergence.csv",
        "ST10_nfkb_pseudobulk_sensitivity.csv",
    ):
        (OUT / obsolete).unlink(missing_ok=True)

    # ---------------------------------------------- carry over ST2-ST5 as-is
    for f in sorted(SRC.glob("ST*.csv")):
        table_number = int(re.match(r"ST(\d+)_", f.name).group(1))
        if table_number > 6 or f.name.startswith("ST1_"):
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
            "CEACAM6 in responders versus non-responders (Fig. 2L); epithelial "
            "metaprogram scores by "
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
    # The composite Tex score of Fig. S4D (S3D until 2026-09-16) is a third list: the ten-gene
    # TEX_GENES of create_S3_D_cd8_tex_score_umap.py (scanpy score_genes).
    tex10 = pd.DataFrame([{
        "Signature": "T cell exhaustion (Tex), 10-gene composite score (Fig. S4D)",
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

    st1 = st1.merge(audit[["Sample", "regimen"]], on="Sample", how="left")
    st1 = st1.merge(
        pairing[["Patient ID", "paired_any_site"]], on="Patient ID", how="left")
    st1 = st1.rename(columns={
        "regimen": "Treatment regimen and cycles",
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
    # Placed where the three columns it replaces used to sit - after the
    # regimen and before the pairing flag - so a reader of the previous version
    # finds the response in the same place.
    st1.insert(st1.columns.get_loc(col),
               "RECIST 1.1 response",
               st1["Patient ID"].map(RECIST_RESPONSE).fillna(""))
    # Patient-level response grouping and treatment apply to every specimen.
    # These 13 rows were blank in the clinician-supplied sheet even though a
    # matched specimen from the same patient carried the value. Use explicit
    # source specimens so P01's pre-treatment regimen is not confused with its
    # later four-cycle post-treatment record.
    fill_from = {
        "P01-M1": "P01-P1", "P02-M1": "P02-P1", "P03-M1": "P03-P1",
        "P04-M1": "P04-P1", "P10-M2": "P10-P2", "P14-M2": "P14-P2",
        "P21-M1": "P21-P1", "P22-M1": "P22-P1", "P25-M1": "P25-P1",
        "P25-B1": "P25-P1", "P26-M1": "P26-P1", "P26-B1": "P26-P1",
        "P34-B2": "P34-P2",
    }
    original_columns = st1.columns.tolist()
    indexed = st1.set_index("Sample")
    for target, source in fill_from.items():
        for field in ("R/NR Grouping", "Treatment regimen and cycles"):
            if pd.isna(indexed.at[target, field]) or indexed.at[target, field] == "":
                indexed.at[target, field] = indexed.at[source, field]
    st1 = indexed.reset_index()[original_columns]
    _check_recist_response(st1)
    st1.to_csv(OUT / "ST1_patient_sample_characteristics.csv", index=False)

    # --------------------------------------------------------- ST6, new
    sweep = pd.read_csv(NEW_ANALYSES / "02_R1.3_TwoSided_Stats_Sweep" / "outputs"
                        / "twosided_sweep.csv")
    st6 = sweep[[
        "panel", "analysis", "family", "test", "group_hi", "n_hi", "group_lo",
        "n_lo", "mean_hi", "mean_lo", "p_two_tailed",
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
    st6["Comparison"] = st6["Comparison"].map(renumber_supp)
    st6.to_csv(OUT / "ST6_two_sided_sensitivity.csv", index=False)

    print("Supplementary tables written to", OUT)
    for f in sorted(OUT.glob("ST*.csv")):
        d = pd.read_csv(f)
        print(f"   {f.name:<48} {d.shape[0]:>4} rows x {d.shape[1]:>2} cols")


def _mutation_test():
    """Show that _check_recist_response() can fail, in the run that uses it.

    Three mutations, one per branch of the gate, each MUST raise: an impossible
    call, a response-labelled patient with no adjudicated call, and a PR
    standing against an NR grouping.
    """
    st1 = pd.read_csv(OUT / "ST1_patient_sample_characteristics.csv")
    st1["RECIST 1.1 response"] = st1["Patient ID"].map(RECIST_RESPONSE).fillna("")
    real = dict(RECIST_RESPONSE)
    mutants = [
        ("an impossible RECIST call", lambda d, f: d.__setitem__("P3", "CRPR")),
        ("a labelled patient with no call", lambda d, f: d.pop("P6")),
        ("a PR call against an NR grouping",
         lambda d, f: (d.__setitem__("P6", "PR"),
                       f.__setitem__("RECIST 1.1 response",
                                     f["Patient ID"].map(d).fillna("")))),
    ]
    failed = []
    for name, mutate in mutants:
        RECIST_RESPONSE.clear()
        RECIST_RESPONSE.update(real)
        frame = st1.copy()
        mutate(RECIST_RESPONSE, frame)
        try:
            _check_recist_response(frame)
        except SystemExit:
            print(f"   mutation raises, as it must: {name}")
        else:
            failed.append(name)
    RECIST_RESPONSE.clear()
    RECIST_RESPONSE.update(real)
    if failed:
        raise SystemExit("\n".join(
            ["build_tables(): _check_recist_response() PASSED a mutant, so it is "
             "not checking what it claims to check:",
             *[f"  - {m}" for m in failed]]))


if __name__ == "__main__":
    main()
    _mutation_test()
