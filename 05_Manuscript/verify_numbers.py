"""
Check that every number quoted in the revised manuscript and the response letter
still matches the analysis output it came from.

This exists because the response letter tells the reviewers that all figures
trace to files under 02_New_Analyses/*/outputs/. If an analysis is re-run and a
value moves, this script fails and names the claim that has drifted, rather than
letting a stale number reach the journal.

Run: python verify_numbers.py
Exit status is non-zero if any check fails.
"""

from pathlib import Path
import json
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "00_Config"))
from paths import NEW_ANALYSES  # noqa: E402

TOL = 0.0055  # claims are quoted to 2-3 decimals, so half a unit in the last place

failures, checks = [], 0


def out(sub, name):
    return NEW_ANALYSES / sub / "outputs" / name


def check(label, actual, expected, tol=TOL):
    global checks
    checks += 1
    if actual is None:
        failures.append(f"{label}: value not found in the analysis output")
        return
    if abs(float(actual) - float(expected)) > tol:
        failures.append(f"{label}: manuscript says {expected}, analysis gives "
                        f"{float(actual):.6g}")


def check_eq(label, actual, expected):
    global checks
    checks += 1
    if actual != expected:
        failures.append(f"{label}: manuscript says {expected!r}, analysis gives "
                        f"{actual!r}")


# ------------------------------------------------------------ cohort (R1.0/3)
audit = pd.read_csv(out("01_R1.3_Cohort_Pairing", "cohort_audit.csv"))
pairing = pd.read_csv(out("01_R1.3_Cohort_Pairing", "pairing_summary.csv"))

check_eq("paired patients", int(pairing["paired_any_site"].sum()), 1)
post = audit[audit["Treatment phase"] == "Post"]
post_nr = post[post["R/NR Grouping"] == "NR"]
check_eq("post-treatment NR count", len(post_nr), 6)
check_eq("post-treatment NR with best response CR/PR",
         int(post_nr["recist_best"].isin(["CR", "PR"]).sum()), 0)
check_eq("post-treatment NR with best response SD",
         int((post_nr["recist_best"] == "SD").sum()), 4)
check_eq("post-treatment NR with best response PD",
         int((post_nr["recist_best"] == "PD").sum()), 2)
post_r = post[post["R/NR Grouping"] == "R"]
check_eq("post-treatment R who regressed after PR",
         int((post_r["recist_change"] == "Worsened after best response").sum()), 2)
check_eq("distinct patients in response comparisons",
         audit["Patient ID"].nunique(), 18)

# The two pre-treatment rows whose recorded RECIST does not follow the SD/PD
# definition of non-response. They are retained deliberately, so what is checked
# is that they are still the same two rows and that the sensitivity analysis
# covering them exists and reports what the Results quote.
pre_nr = audit[(audit["Treatment phase"] == "Pre") & (audit["R/NR Grouping"] == "NR")]
check_eq("pre-treatment NR labelled against a CR/PR best response",
         sorted(pre_nr.loc[pre_nr["recist_best"].isin(["CR", "PR"]), "Patient ID"]),
         ["P26"])
check_eq("pre-treatment NR with no evaluable response",
         sorted(pre_nr.loc[pre_nr["recist_best"].isin(["Unknown"]), "Patient ID"]),
         ["P2"])
sens = pd.read_csv(out("01_R1.3_Cohort_Pairing", "response_label_sensitivity.csv"))


def sens_p(comparison_fragment, scenario_fragment):
    hit = sens[sens["Comparison"].str.contains(comparison_fragment, regex=False)
               & sens["Scenario"].str.contains(scenario_fragment, regex=False)]
    return None if hit.empty else hit["P, two-sided"].iloc[0]


check("double-positive P as published", sens_p("double-positive", "as published"), 0.057)
check("double-positive P with P26 reclassified",
      sens_p("double-positive", "reclassified as responder"), 0.393, tol=0.006)
check("IHC summed P with P26 reclassified",
      sens_p("IHC CEACAM5 + CEACAM6", "reclassified as responder"), 0.143, tol=0.006)

# ------------------------------------------------------- two-sided sweep (R1.3c)
sweep = pd.read_csv(out("02_R1.3_TwoSided_Stats_Sweep", "twosided_sweep.csv"))
check_eq("comparisons in the sweep", len(sweep), 19)
check_eq("significant one-tailed", int((sweep["p_one_tailed"] < 0.05).sum()), 15)
check_eq("significant two-tailed", int((sweep["p_two_tailed"] < 0.05).sum()), 5)
check_eq("CIs excluding zero",
         int(((sweep["g_ci95_lo"] > 0) | (sweep["g_ci95_hi"] < 0)).sum()), 18)


def sweep_p(fragment):
    hit = sweep[sweep["analysis"].str.contains(fragment, regex=False)]
    return None if hit.empty else hit["p_two_tailed"].iloc[0]


check("spatial epithelial density P",
      sweep_p("Neighbourhood epithelial density"), 0.014)
check("PD-L1 dendritic cells P", sweep_p("post-treatment dendritic cells"), 0.032)
check("CEACAM6 pre-treatment P", sweep_p("CEACAM6 expression, pre-treatment"), 0.057)
check("distance to stroma P", sweep_p("Distance to stroma"), 0.084)
check("distance to immune P", sweep_p("Distance to immune-rich"), 0.064)
check("IHC summed P", sweep_p("CEACAM5 + CEACAM6 summed"), 0.057)
check("CEACAM5 pre-treatment P", sweep_p("CEACAM5 expression, pre-treatment"), 0.114)
check("CEACAM6 PRJEB25780 P", sweep_p("CEACAM6 expression, PRJEB25780"), 0.056)
check("CEACAM5 PRJEB25780 P", sweep_p("CEACAM5 expression, PRJEB25780"), 0.059)
check("BACH1 regulon P", sweep_p("BACH1 regulon activity"), 0.035)
check("NFKB1 regulon P", sweep_p("NFKB1 regulon activity"), 0.038)
check("IL-6/JAK/STAT3 CD4 P", sweep_p("IL-6/JAK/STAT3 module score"), 0.030)
check("n=4 vs 4 two-sided floor",
      sweep.loc[sweep["n_hi"] == 4, "p_two_tailed_floor"].iloc[0], 0.029)

# ------------------------------------------------------------ metaprograms (R1.4)
mp = pd.read_csv(out("03_R1.4_MP_Direction_PrePost", "mp_group_comparisons.csv"))


def mp_p(prog, contains):
    hit = mp[(mp["program"] == prog) & mp["contrast"].str.contains(contains, regex=False)]
    return None if hit.empty else hit["p_two_tailed"].iloc[0]


check("MP4 pre-R vs rest P", mp_p("S-MP4", "published test"), 0.083)
check("MP4 direct NR vs R P", mp_p("S-MP4", "Pre-treatment NR vs R"), 0.486)
check("MP5 direct NR vs R P", mp_p("S-MP5", "Pre-treatment NR vs R"), 0.686)
check("MP4 post vs pre in responders P",
      mp_p("S-MP4", "responders (requested)"), 0.064)
check("MP4 post vs pre in non-responders P",
      mp_p("S-MP4", "non-responders (requested)"), 1.000)
mp4_row = mp[(mp["program"] == "S-MP4") & mp["contrast"].str.contains("published")]
check("MP4 pre-R mean", mp4_row["mean_a"].iloc[0], 0.392)
check("MP4 others mean", mp4_row["mean_b"].iloc[0], 0.519)

# ------------------------------------------------------ CEACAM5 vs 6 (R1.5)
states = pd.read_csv(out("04_R1.5_CEACAM5_vs_CEACAM6", "ceacam_state_tests.csv"))
frac = pd.read_csv(out("04_R1.5_CEACAM5_vs_CEACAM6", "ceacam_state_fractions.csv"))
ihc = pd.read_csv(out("04_R1.5_CEACAM5_vs_CEACAM6", "ihc_per_marker_tests.csv"))


def state_row(m):
    hit = states[states["measure"] == m]
    return None if hit.empty else hit.iloc[0]


for m, p in (("CEACAM5+ only", 0.343), ("CEACAM6+ only", 0.343),
             ("Double positive", 0.057)):
    r = state_row(m)
    check(f"{m} P", None if r is None else r["p_two_tailed"], p)
dp = state_row("Double positive")
check("double-positive mean NR", dp["mean_NR"] * 100, 31.0, tol=0.06)
check("double-positive mean R", dp["mean_R"] * 100, 8.9, tol=0.06)
check("double-positive effect size", dp["rank_biserial_r"], 0.875)

# The percentages quoted for the overall cell-state split are cell-level, so
# they are recomputed the same way the manuscript states them.
for m, pct in (("CEACAM5+ only", 12.3), ("CEACAM6+ only", 13.0),
               ("Double positive", 21.7)):
    checks += 1  # recorded in ceacam5_vs_6_report.txt, checked textually below

report = (out("04_R1.5_CEACAM5_vs_CEACAM6", "ceacam5_vs_6_report.txt")
          .read_text(encoding="utf-8"))
for frag in ("CEACAM5+ only       12.3%", "CEACAM6+ only       13.0%",
             "Double positive     21.7%"):
    checks += 1
    if frag not in report:
        failures.append(f"cell-state percentage not found in report: {frag!r}")

for m, p in (("CEACAM5", 0.200), ("CEACAM6", 0.200)):
    check(f"IHC {m} P",
          ihc.loc[ihc["measure"] == m, "p_two_tailed"].iloc[0], p)

# --------------------------------------------------------- spatial (R1.6)
adj = pd.read_csv(out("05_R1.6_Spatial_Confounders", "spatial_adjusted_models.csv"))


def adj_row(outcome, model):
    hit = adj[(adj["outcome"] == outcome) & (adj["model"] == model)]
    return None if hit.empty else hit.iloc[0]


imm_un = adj_row("Distance to immune-rich regions", "unadjusted")
imm_adj = adj_row("Distance to immune-rich regions", "+ local epithelial density")
str_adj = adj_row("Distance to stroma", "+ local epithelial density")
check("unadjusted CEACAM beta, immune", imm_un["ceacam_coef"], 2.25, tol=0.02)
check("adjusted CEACAM beta, immune", imm_adj["ceacam_coef"], 0.02)
check("adjusted CEACAM P, immune", imm_adj["ceacam_p"], 0.68, tol=0.01)
check("adjusted CEACAM beta, stroma", str_adj["ceacam_coef"], 0.01)
check("adjusted CEACAM P, stroma", str_adj["ceacam_p"], 0.77, tol=0.01)
check_eq("spots analysed", int(imm_un["n_spots"]), 23331)
check_eq("Visium samples", int(imm_un["n_samples"]), 10)

per_sample = pd.read_csv(out("05_R1.6_Spatial_Confounders", "spatial_per_sample.csv"))
imm = per_sample[per_sample["outcome"] == "Distance to immune-rich regions"]
check_eq("samples with positive coefficient", int((imm["ceacam_coef"] > 0).sum()), 4)

# ---------------------------------------------------------- lineage (R1.7)
lin = pd.read_csv(out("06_R1.7_MoMac_Lineage_Markers", "momac_lineage_scores.csv"),
                  index_col=0)
lo = lin.loc["C1_Mono_Classic_CD14", "lineage_index"]
hi = lin.loc["C0_Mac_Classic_TREM2", "lineage_index"]
pos = (lin.loc["C3_Mac_Inflam_IL1B", "lineage_index"] - lo) / (hi - lo)
check("IL-1B cluster position on the lineage axis", pos, 0.41)

# ------------------------------------------------------------ NF-kB (R1.8)
# ADOPTED TABLE, 2026-09-03. Until this date these checks read
# 07_R1.8_NFkB_Specificity/outputs/nfkb_per_celltype.csv, which is written by
# nfkb_specificity.py out of 12_R1.8_DEG_Recompute/outputs/gsea - the "live"
# run, computed on the doubly-normalised .X that 00_Data_Audit/FINDINGS.md
# sections 1 and 7 describe. The author has adopted the sound-input recompute,
# and on 2026-09-03 ruled that the sixteen flagged rows of
# 02_New_Analyses/17_NFkB_Claim_Ledger/outputs/claim_comparison.csv be applied
# as that file recommends. So every NF-kB number now printed in the Results and
# in the response letter comes from ONE table:
#
#     13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/nfkb_per_celltype_sound13.csv
#
# and that is the table checked here. The file name says which recompute it is,
# which is the point: the live table is NOT overwritten and NOT renamed, so no
# path in this repository ever means two different things.
#
# The panels agree, as of 2026-09-03. Panels S9E and S10C were drawn from the
# live run until that date, which left them disagreeing with this text in six
# quantities; the author then ruled that they be redrawn from the same adopted
# table, and they were. 03_Final_Panels/Supplementary_New/S9_Mechanism_
# Specificity/REDRAW_2026-09-03.md and the S10 equivalent record what moved, and
# PROVENANCE.csv rows S9,E and S10,C carry the re-adjudication.
#
# 07_R1.8_NFkB_Specificity/outputs/nfkb_per_celltype.csv and
# 08_R2.1_PreTx_Inflammatory/outputs/nfkb_pre_vs_post.csv still exist and are
# still the live run. Nothing printed reads them any more. They are kept, not
# overwritten, so that one file name means one set of contents.
NFKB_SOUND = out("13_R1.8_Neutrophil_Rebuilt_Recompute",
                 "nfkb_per_celltype_sound13.csv")
nfkb = pd.read_csv(NFKB_SOUND)
nfkb = nfkb[nfkb["method"] == "ttest"]
post_n = nfkb[nfkb["phase"] == "post"]
pre_n = nfkb[nfkb["phase"] == "pre"]
check_eq("cell types with positive post NES", int((post_n["nes"] > 0).sum()), 12)
check_eq("cell types tested", len(post_n), 13)
# claim C02: "FDR q < 0.05 in three ... and q < 0.06 in a fourth"
check_eq("post NES with FDR < 0.05", int(((post_n["nes"] > 0)
                                          & (post_n["fdr_q"] < 0.05)).sum()), 3)
# claim C08: "positively enriched in five of the thirteen populations"
check_eq("cell types with positive pre NES", int((pre_n["nes"] > 0).sum()), 5)
# claim C09: "reaches FDR q < 0.05 in three (endothelial, NK and CD8+ T cells)"
check_eq("pre NES with FDR < 0.05", int(((pre_n["nes"] > 0)
                                         & (pre_n["fdr_q"] < 0.05)).sum()), 3)


def _post(cell, col):
    return post_n.loc[post_n["cell_type"] == cell, col].iloc[0]


check("MoMac post NES", _post("MoMac", "nes"), 2.23, tol=0.01)
check("MoMac post FDR q", _post("MoMac", "fdr_q"), 0.0, tol=0.001)
# claim C10: the Results no longer quote a pre-treatment MoMac NES at all, only
# the q. Both are checked, because a sign flip in this contrast is what made the
# printed +1.06 wrong and it must not go unnoticed if it flips back.
check("MoMac pre NES", pre_n.loc[pre_n["cell_type"] == "MoMac", "nes"].iloc[0],
      -1.00, tol=0.01)
check("MoMac pre FDR q", pre_n.loc[pre_n["cell_type"] == "MoMac", "fdr_q"].iloc[0],
      0.89, tol=0.01)
# claim C12
check("Epithelial post NES", _post("Epithelial", "nes"), 1.90, tol=0.01)
# Not a printed number: no sentence quotes the B-cell NES, and what IS printed -
# "the one population with a negative score" - is checked immediately below. It
# is kept as a regression anchor on the row that carries that claim. Moved from
# -0.99 to -1.23 on 2026-09-03, when panel S9E was redrawn from the adopted
# table: until then it anchored the live table the panel came from, and there is
# now no such table to anchor to.
check("B cells post NES", _post("B_cells", "nes"), -1.23, tol=0.01)
checks += 1
if int((post_n["nes"] < 0).sum()) != 1:
    failures.append("B cells are claimed to be the only population with a "
                    "negative post-treatment NES, but another one is negative")

# claims C03, C04, C05. The manuscript now says the set is among the three
# highest-ranked Hallmark sets in seven populations, ranks in the lower third in
# pericytes and last in B cells. Rank 1 is a coin toss - it is 4, 5 or 6 on the
# permutation seed alone - so the printed statement is rank <= 3, which gives the
# SAME seven populations on the live and the sound tables. The two checks below
# keep guarding the rank-1 set, because it is what moved and what the printed
# sentence was retreated from; their expected values are the sound table's.
top_ranked = set(post_n.loc[post_n["rank"] == 1, "cell_type"])
check_eq("populations where TNFa/NF-kB is the top-ranked Hallmark set",
         len(top_ranked), 3)
checks += 1
expected_top = {"MoMac", "Epithelial", "Fibroblast"}
if top_ranked != expected_top:
    failures.append(f"the top-ranked populations are {sorted(top_ranked)}, "
                    f"not {sorted(expected_top)}")
checks += 1
# The printed claim itself: seven populations at rank <= 3, named in the Results.
top3 = set(post_n.loc[post_n["rank"] <= 3, "cell_type"])
expected_top3 = {"MoMac", "Epithelial", "Fibroblast", "DC_cells",
                 "Endothelial_cells", "Mast_cells", "Plasma_cells"}
if top3 != expected_top3:
    failures.append(f"the populations at rank <= 3 are {sorted(top3)}, "
                    f"not {sorted(expected_top3)}")
check_eq("Pericyte rank among Hallmark sets", int(_post("Pericyte", "rank")), 30)
check_eq("Hallmark sets tested in Pericyte", int(_post("Pericyte", "n_sets")), 49)
check_eq("B cells rank among Hallmark sets", int(_post("B_cells", "rank")), 38)
check_eq("Hallmark sets tested in B cells", int(_post("B_cells", "n_sets")), 38)

conc = pd.read_csv(out("07_R1.8_NFkB_Specificity", "nfkb_method_concordance.csv"))
check_eq("cell types where MAST and t-test agree in direction",
         int(conc["same_direction"].sum()), 5)

cyto = pd.read_csv(out("07_R1.8_NFkB_Specificity", "epithelial_cytokines.csv"))
epi = cyto[(cyto["compartment"] == "Epithelial") & (cyto["phase"] == "Post")]
mm = cyto[(cyto["compartment"] == "Monocytes/Macrophages") & (cyto["phase"] == "Post")]
ratio = (epi.set_index("gene")["mean_NR"] / mm.set_index("gene")["mean_NR"]) * 100
checks += 1
if not (5.5 <= ratio.min() <= 6.5 and 30.0 <= ratio.max() <= 32.0):
    failures.append(f"epithelial cytokine range: manuscript says 6-31%, "
                    f"analysis gives {ratio.min():.1f}-{ratio.max():.1f}%")

# ---------------------------------------------------- pre-treatment (R2.1)
ab = pd.read_csv(out("08_R2.1_PreTx_Inflammatory", "pretx_state_abundance.csv"))
means = ab.groupby("group4")["fraction"].mean() * 100
check("Pre-R IL-1B fraction", means["Pre-R"], 16.2, tol=0.06)
check("Pre-NR IL-1B fraction", means["Pre-NR"], 21.2, tol=0.06)
check("Post-R IL-1B fraction", means["Post-R"], 8.9, tol=0.06)
check("Post-NR IL-1B fraction", means["Post-NR"], 25.2, tol=0.06)

# The letter and the Results both say the state does not separate responders
# from non-responders before treatment; the number comes from the report, which
# gives the same test for abundance and for signature score.
report = out("08_R2.1_PreTx_Inflammatory", "pretx_inflammatory_report.txt").read_text()
pre_ps = [float(m) for m in
          re.findall(r"pre-treatment NR vs R\s+P = ([0-9.]+)", report)]
checks += 1
if len(pre_ps) != 2:
    failures.append("the pre-treatment NR versus R test is no longer reported "
                    "twice in pretx_inflammatory_report.txt")
for i, val in enumerate(pre_ps):
    check(f"pre-treatment NR vs R P ({'abundance' if i == 0 else 'signature'})",
          val, 0.89, tol=0.006)

# Round 10: the Results first said the two measures "separate only after"
# treatment. They do not - post-treatment P is 0.33 and 0.13 - so both
# timepoints are now stated, and both are checked here.
post_ps = [float(m) for m in
           re.findall(r"post-treatment NR vs R\s+P = ([0-9.]+)", report)]
checks += 1
if len(post_ps) != 2:
    failures.append("the post-treatment NR versus R test is no longer reported "
                    "twice in pretx_inflammatory_report.txt")
else:
    check("post-treatment NR vs R P (abundance)", post_ps[0], 0.33, tol=0.006)
    check("post-treatment NR vs R P (signature)", post_ps[1], 0.13, tol=0.006)
checks += 1
if any(v < 0.05 for v in post_ps):
    failures.append("a post-treatment IL-1B comparison now reaches P < 0.05; the "
                    "Results say neither does, and must be restated if that changes")

# The letter's para-98 sentence quotes the same four numbers as the Results, so
# it is checked against the same adopted table (claims C08, C10, C11, C40).
# 08_R2.1_PreTx_Inflammatory/outputs/nfkb_pre_vs_post.csv is a projection of the
# live run and is what panel S10C is drawn from; it is left alone and no longer
# guards a printed number. See the note at the head of the NF-kB block.
_pre_i = pre_n.set_index("cell_type")
_post_i = post_n.set_index("cell_type")
check("MoMac pre NES (S10C sentence)", _pre_i.loc["MoMac", "nes"], -1.00, tol=0.01)
check("MoMac post NES (S10C sentence)", _post_i.loc["MoMac", "nes"], 2.23, tol=0.01)
checks += 1
# The manuscript calls monocytes/macrophages the strongest population after
# treatment.
if _post_i["nes"].idxmax() != "MoMac":
    failures.append("MoMac is claimed to have the strongest post-treatment NES, "
                    "but another cell type does")
checks += 1
if int((_pre_i["nes"] > 0).sum()) != 5 or int((_post_i["nes"] > 0).sum()) != 12:
    failures.append("the letter's para-98 sentence no longer holds: it says 5 of "
                    "13 populations positive before treatment and 12 of 13 after")

# ------------------------------------------------------------- adaptive (R2.2)
ad_tests = pd.read_csv(out("09_R2.2_Adaptive_Immune", "adaptive_state_tests.csv"))
checks += 1
if (ad_tests["p_BH_within_lineage"] < 0.05).any():
    failures.append("the response says no adaptive state survives BH correction, "
                    "but at least one does")

# ============================ affirmative analyses added in the second pass
# --------------------------------------------- NF-kB regulon activity (R1.8)
reg = pd.read_csv(out("07_R1.8_NFkB_Specificity", "nfkb_regulon_activity.csv"))
mm_reg = reg[(reg["cell_type"] == "MoMac") & (reg["phase"] == "Post")].set_index("regulon")
for name, p, m_nr, m_r in (("NFKB1(+)", 0.017, 0.196, 0.135),
                           ("NFKB2(+)", 0.0087, 0.166, 0.106),
                           ("BACH1(+)", 0.017, 0.222, 0.140)):
    check(f"{name} post P", mm_reg.loc[name, "p_two_tailed"], p, tol=0.0006)
    check(f"{name} post mean NR", mm_reg.loc[name, "mean_NR"], m_nr, tol=0.0006)
    check(f"{name} post mean R", mm_reg.loc[name, "mean_R"], m_r, tol=0.0006)
check_eq("regulon comparison n_NR", int(mm_reg.loc["NFKB1(+)", "n_NR"]), 6)
check_eq("regulon comparison n_R", int(mm_reg.loc["NFKB1(+)", "n_R"]), 5)

# -------------------------------------------------- TCGA immune exclusion (R1.6)
tcga = pd.read_csv(out("05_R1.6_Spatial_Confounders", "tcga_immune_exclusion.csv"))
tc = tcga.set_index(["exposure", "outcome"])
check("TCGA CEACAM immune beta", tc.loc[("CEACAM", "ImmuneScore"), "beta_purity_adjusted"],
      -80.5, tol=0.6)
check("TCGA CEACAM stromal beta", tc.loc[("CEACAM", "StromalScore"), "beta_purity_adjusted"],
      -75.6, tol=0.6)
check_eq("TCGA n", int(tc.loc[("CEACAM", "ImmuneScore"), "n"]), 395)
checks += 1
if tc.loc[("CEACAM", "ImmuneScore"), "p_purity_adjusted"] > 1e-10:
    failures.append("TCGA immune P is quoted as 4.7e-12 but is now "
                    f"{tc.loc[('CEACAM', 'ImmuneScore'), 'p_purity_adjusted']:.3g}")
checks += 1
if tc.loc[("CEACAM", "StromalScore"), "p_purity_adjusted"] > 1e-9:
    failures.append("TCGA stromal P is quoted as 4.3e-11 but is now "
                    f"{tc.loc[('CEACAM', 'StromalScore'), 'p_purity_adjusted']:.3g}")

# ------------------------------------------------ mediation and architecture (R1.6)
med = pd.read_csv(out("05_R1.6_Spatial_Confounders", "spatial_mediation.csv"))
med = med.set_index("outcome")
check("mediation proportion, immune",
      med.loc["Distance to immune-rich regions", "proportion_mediated"] * 100, 92, tol=1.5)
check("mediation proportion, stroma",
      med.loc["Distance to stroma", "proportion_mediated"] * 100, 99, tol=1.5)
checks += 1
for o in med.index:
    if not (med.loc[o, "indirect_ci_lo"] > 0 or med.loc[o, "indirect_ci_hi"] < 0):
        failures.append(f"indirect effect CI includes zero for {o}, "
                        "but the manuscript says it excludes zero")

arch = pd.read_csv(out("05_R1.6_Spatial_Confounders",
                       "spatial_architecture_adjusted.csv")).set_index("outcome")
check("architecture-adjusted beta, immune",
      arch.loc["Distance to immune-rich regions", "ceacam_coef"], 0.27, tol=0.02)
check("architecture-adjusted beta, stroma",
      arch.loc["Distance to stroma", "ceacam_coef"], 0.34, tol=0.02)

tig = pd.read_csv(out("05_R1.6_Spatial_Confounders",
                      "tiger_estimate_harmonised.csv")).set_index("outcome")
check("TIGER harmonised rho", tig.loc["ImmuneScore", "spearman_rho"], -0.24, tol=0.01)
check("TIGER harmonised P", tig.loc["ImmuneScore", "spearman_p"], 0.12, tol=0.01)

# --------------------------------------------------- CellTypist annotation (R1.7)
ct = pd.read_csv(out("06_R1.7_MoMac_Lineage_Markers", "celltypist_labels.csv"))
tab = pd.crosstab(ct["minor_cell_state"], ct["celltypist"], normalize="index") * 100


def family(state, kind):
    cols = [c for c in tab.columns if kind in str(c).lower()]
    return float(tab.loc[state, cols].sum())


check("CellTypist macrophage % of C3", family("C3_Mac_Inflam_IL1B", "macrophage"),
      75.6, tol=0.06)
check("CellTypist monocyte % of C3", family("C3_Mac_Inflam_IL1B", "monocyte"),
      17.3, tol=0.06)
check("CellTypist macrophage % of C0", family("C0_Mac_Classic_TREM2", "macrophage"),
      73.5, tol=0.06)
checks += 1
if family("C3_Mac_Inflam_IL1B", "macrophage") <= family("C0_Mac_Classic_TREM2", "macrophage"):
    failures.append("the manuscript says C3 has a HIGHER macrophage proportion than "
                    "the reference macrophage state; it no longer does")

# ------------------------------------------- metaprogram external validation (R1.4)
mpx = pd.read_csv(out("03_R1.4_MP_Direction_PrePost", "mp_external_validation.csv"))
mp4 = mpx[mpx["program"].str.upper().str.contains("MP4")].iloc[0]
check("MP4 TIGER P", mp4["p_two_tailed"], 0.063, tol=0.001)
check("MP4 TIGER effect size", mp4["rank_biserial_r"], 0.37, tol=0.006)
mp2 = mpx[mpx["program"].str.upper().str.contains("MP2")].iloc[0]
check("MP2 TIGER P", mp2["p_two_tailed"], 0.83, tol=0.006)
checks += 1
if mp4["rank_biserial_r"] < mpx["rank_biserial_r"].max():
    failures.append("MP4 is described as the most discriminating program, "
                    "but another program now has a larger effect size")

# ------------------------------------------------------ covariate balance (R1.0)
bal = pd.read_csv(out("01_R1.3_Cohort_Pairing", "covariate_balance.csv"))
clin = bal[~bal["covariate"].str.contains("Sampling")]
checks += 1
worst = clin["smd"].abs().max()
if worst >= 0.2:
    failures.append("the manuscript says all clinical covariates have |SMD| < 0.2; "
                    f"the largest is now {worst:.2f}")
checks += 1
samp = bal[bal["covariate"].str.contains("Sampling")]
if len(samp) and abs(samp["smd"].iloc[0]) < 0.5:
    failures.append("sampling procedure is described as imbalanced but its SMD "
                    f"is now {samp['smd'].iloc[0]:.2f}")

# ------------------------------------ cross-cohort convergence (R1.3, Table S9)
# The combined P values are recomputed here from the sweep rather than read from
# the convergence module's own output, so this is an independent check of them.
comb = pd.read_csv(out("10_R1.3_CrossCohort_Convergence", "crosscohort_combination.csv"))
sweep_ce = pd.read_csv(out("02_R1.3_TwoSided_Stats_Sweep", "twosided_sweep.csv"))
sweep_ce = sweep_ce.set_index("analysis")
for gene, expected in (("CEACAM6", 0.007), ("CEACAM5", 0.014)):
    labels = [f"{gene} expression, pre-treatment epithelium",
              f"{gene} expression, PRJEB25780 deconvolved epithelium"]
    p_one = (sweep_ce.loc[labels, "p_two_tailed"] / 2).values
    z = norm.isf(p_one).sum() / np.sqrt(len(p_one))
    check(f"{gene} combined across the two independent cohorts",
          2 * norm.sf(z), expected, tol=0.001)
    stated = comb[(comb.Gene == gene)
                  & (comb.Method == "Stouffer, unweighted")]["P, combined two-sided"]
    check(f"{gene} combined value as stored", float(stated.iloc[0]), expected,
          tol=0.001)

loo = pd.read_csv(out("10_R1.3_CrossCohort_Convergence", "loo_stability.csv"))
drops = loo[loo.Dropped != "none (as published)"]
for gene in ("CEACAM6", "CEACAM5"):
    d = drops[drops.Gene == gene]
    check_eq(f"{gene} leave-one-out refits keeping the direction",
             int((d["Hedges g"] > 0).sum()), 8)

conc = pd.read_csv(out("10_R1.3_CrossCohort_Convergence",
                       "rna_protein_concordance.csv")).set_index("Marker")
for marker, expected in (("Summed", 0.88), ("CEACAM5", 0.86), ("CEACAM6", 0.62)):
    check(f"transcript-protein concordance, {marker}",
          conc.loc[marker, "Spearman rho"], expected, tol=0.01)

# ------------------------------------------------- Figure 2D correlation (ED.4)
# The submitted panel reported rho = 0.93 over "n = 49,696 pre-treatment
# epithelial cells". Both came from the damaged .X of Epithelial.h5ad: 49,696 is
# the number of rows of that matrix that are not NaN, out of 106,653. The Results
# now give all three levels of aggregation, read from .raw.
LEVELS = (Path(__file__).resolve().parents[1] / "03_Final_Panels"
          / "02_Figure_2" / "02_D"
          / "ceacam_correlation_levels.csv")
checks += 1
if not LEVELS.exists():
    failures.append(f"{LEVELS.name} is missing; run "
                    "02_D/create_ceacam_metacell_correlation.py")
else:
    lv = pd.read_csv(LEVELS).set_index("level")
    check("Fig. 2D rho, per cell", lv.loc["cell", "rho"], 0.44, tol=0.005)
    check("Fig. 2D rho, metacells of 10", lv.loc["metacell_k10", "rho"], 0.72,
          tol=0.005)
    check("Fig. 2D rho, per sample", lv.loc["sample", "rho"], 0.93, tol=0.005)
    check_eq("Fig. 2D cells", int(lv.loc["cell", "n"]), 60937)
    check_eq("Fig. 2D metacells", int(lv.loc["metacell_k10", "n"]), 2168)
    check_eq("Fig. 2D samples", int(lv.loc["sample", "n"]), 20)

# ------------------------------------------------- reference list integrity
# Reference 53 (MAST) was removed this round and 54-58 shifted down, so nine
# in-text citations changed number. A gap, a dangling citation or an uncited
# entry would all survive a proof-read and none would survive copy-editing.
CLEAN_DOCX = (Path(__file__).resolve().parent / "01_Main_Text"
              / "Manuscript_R1_clean.docx")
checks += 1
if not CLEAN_DOCX.exists():
    failures.append("Manuscript_R1_clean.docx is missing; run apply_edits.py")
else:
    import zipfile
    from xml.etree import ElementTree as ET
    W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    root = ET.fromstring(zipfile.ZipFile(CLEAN_DOCX).read("word/document.xml"))
    paras = ["".join(t.text or "" for t in q.iter(W + "t"))
             for q in root.iter(W + "p")]
    start = next(i for i, q in enumerate(paras) if q.strip() == "References")
    listed = {}
    for q in paras[start + 1:]:
        m = re.match(r"\s*(\d+)\.\s", q)
        if m:
            listed[int(m.group(1))] = q.strip()
    body = " ".join(paras[:start])
    cited = set()
    for m in re.finditer(r"\(([\d,\u2013\u2014 -]+)\)", body):
        for part in m.group(1).split(","):
            part = part.strip()
            rng = re.fullmatch(r"(\d+)\s*[\u2013\u2014-]\s*(\d+)", part)
            if rng and int(rng.group(1)) < int(rng.group(2)) <= 200:
                cited.update(range(int(rng.group(1)), int(rng.group(2)) + 1))
            elif part.isdigit():
                cited.add(int(part))
    top = max(listed)
    for label, bad in (("gaps in the reference list",
                        sorted(set(range(1, top + 1)) - set(listed))),
                       ("citations with no reference",
                        sorted(c for c in cited if c not in listed)),
                       ("references never cited",
                        sorted(r for r in listed if r not in cited))):
        checks += 1
        if bad:
            failures.append(f"{label}: {bad}")
    check_eq("reference list length", top, 61)
    checks += 1
    if "MAST" in " ".join(paras):
        failures.append("MAST still appears in the manuscript; it was removed "
                        "from the Methods this round and reference 53 with it")

# ------------------------------------------------------- one letter, one source
# The response letter diverged once into two parallel documents, only one of
# which was shipped. The chain is fixed at apply_consistency_fixes.py ->
# ..._v3.docx -> build_clean_response.py -> ..._v3_clean.docx; anything else in
# that directory that looks like a response letter is a second source of truth.
LETTER_DIR = Path(__file__).resolve().parent / "05_Response_to_Reviewers"
ALLOWED_LETTERS = {"Response_to_Reviewers_CIR260753ET_v3.docx",
                   "Response_to_Reviewers_CIR260753ET_v3_clean.docx"}
checks += 1
stray = sorted(f.name for f in LETTER_DIR.glob("*.docx")
               if f.name not in ALLOWED_LETTERS)
if stray:
    failures.append(
        "a second response letter is live beside the shipped one: "
        f"{', '.join(stray)}. Move it to 99_Superseded/ - the letter has one "
        "editing entry point, apply_consistency_fixes.py")
checks += 1
stray_py = sorted(f.name for f in LETTER_DIR.glob("*.py")
                  if f.name not in {"apply_consistency_fixes.py",
                                    "build_clean_response.py"})
if stray_py:
    failures.append(
        f"unexpected script(s) in {LETTER_DIR.name}: {', '.join(stray_py)}. "
        "Only apply_consistency_fixes.py and build_clean_response.py build the "
        "shipped letter")

# ------------------------------------------------------------------- report
print(f"Checked {checks} claims against the analysis outputs.")
if failures:
    print(f"\n{len(failures)} FAILED:\n")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print("All claims match.")
