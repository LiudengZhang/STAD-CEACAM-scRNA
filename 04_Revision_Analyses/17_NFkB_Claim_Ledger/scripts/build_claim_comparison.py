#!/usr/bin/env python
"""Build CLAIM_COMPARISON.md and outputs/claim_comparison.csv from one row list.

Read-only with respect to everything outside
02_New_Analyses/17_NFkB_Claim_Ledger/. Every number in ROWS was read from a table
already on disk; nothing here re-runs GSEA, DEG or a panel.

Sources, all read not written:
  02_New_Analyses/13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/
        nfkb_per_celltype_sound13.csv        the sound-input recompute (axis 1)
        nfkb_per_celltype_live.csv           the shipped table
        counted_claims.csv                   the counted claims on both
  02_New_Analyses/17_NFkB_Claim_Ledger/outputs/
        signal.csv, signal_runs.csv          the gate and the seed sweeps (axis 2)
        stability.csv, counted_claims_by_seed.csv, ordinal_denominator.csv
        para70_summary.csv                   Fig. 5N / Results para 71
  02_New_Analyses/14_MAST_Specification/outputs/nfkb_by_specification.csv  (reserve)
  03_Revised_Panels/SUPPLEMENTARY_AUDIT.md   fault 1, the S10E measurement
  04_Manuscript_R1/verify_numbers.py         the guard column
  04_Manuscript_R1/01_Main_Text/Manuscript_R1_clean.docx
  04_Manuscript_R1/05_Response_to_Reviewers/Response_to_Reviewers_CIR260753ET_v3_clean.docx
"""
import csv
import pathlib

HERE = pathlib.Path(__file__).resolve().parents[1]
CSV_OUT = HERE / "outputs" / "claim_comparison.csv"
MD_OUT = HERE / "CLAIM_COMPARISON.md"

S13 = ("13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/nfkb_per_celltype_sound13.csv")
CC = ("13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/counted_claims.csv")
SIG = "17_NFkB_Claim_Ledger/outputs/signal.csv"
P70 = "17_NFkB_Claim_Ledger/outputs/para70_summary.csv"
AUDIT = "03_Revised_Panels/SUPPLEMENTARY_AUDIT.md fault 1"

F = ["claim_id", "ledger_id", "source", "location", "verbatim", "printed_value",
     "sound_input_value", "sound_input_table", "axis1_input", "axis2_stability",
     "gate", "verify_numbers_guard", "guard_expected_value", "recommendation",
     "proposed_wording", "supporting_measurement", "note"]

R = []


def row(**kw):
    d = {k: "" for k in F}
    d.update(kw)
    missing = set(kw) - set(F)
    assert not missing, missing
    R.append(d)


# ============================================================ Results, NF-kB para
# clean.docx paragraph 69 (1-based over all w:p). ledger.csv calls it para 68;
# ledger.csv's paragraph numbers are 0-based. Both are given throughout.
row(claim_id="C01", ledger_id="M01 / R01a", source="main text",
    location="clean.docx para 69 s2 (ledger para 68); edits.py:361-364",
    verbatim="Systematic GSEA across the 13 major cell types showed positive "
             "enrichment of the Hallmark TNFα signaling via NF-κB gene set in "
             "non-responders in 12 of 13 populations after treatment,",
    printed_value="12 of 13", sound_input_value="12 of 13",
    sound_input_table=S13 + " (ttest, phase=post); " + CC + " row sound13",
    axis1_input="unchanged", axis2_stability="safe",
    gate="8 of 13 post contrasts signal-bearing on sound13; of the 12 positives, "
         "8 yes / 2 weak / 2 no",
    verify_numbers_guard="yes - verify_numbers.py:228 check_eq('cell types with "
                         "positive post NES'); :229 ('cell types tested')",
    guard_expected_value="12; 13", recommendation="keep verbatim",
    proposed_wording="",
    supporting_measurement="12 of 13 on live, on sound13, at all 7 seeds, under all "
                           "four ranking metrics, and 13 of 13 at the sample level",
    note="Survives adoption with no edit and no change to its guard.")

row(claim_id="C02", ledger_id="M02 / R01b", source="main text",
    location="clean.docx para 69 s2 (cont.); edits.py:364-366",
    verbatim="reaching FDR q < 0.05 in four (monocytes/macrophages, dendritic "
             "cells, epithelial cells and fibroblasts; Fig. S9E).",
    printed_value="four (MoMac, DC, Epithelial, Fibroblast)",
    sound_input_value="three (MoMac, Epithelial, Fibroblast); DC_cells q = 0.0517",
    sound_input_table=S13 + " (ttest, post, fdr_q); " + CC + " post_n_positive_q05",
    axis1_input="moved (4 -> 3; DC_cells drops out at q = 0.0517)",
    axis2_stability="floor - 4 / 6 / 8 / 9 by ranking metric, 8 at sample level "
                    "(both methods), 3 on sound inputs, 3 under MAST spec A",
    gate="all four named contrasts are signal-bearing on both tables",
    verify_numbers_guard="yes - verify_numbers.py:231 check_eq('post NES with FDR "
                         "< 0.05')",
    guard_expected_value="4  -> must become 3 on adoption",
    recommendation="restate",
    proposed_wording="reaching FDR q < 0.05 in three (monocytes/macrophages, "
                     "epithelial cells and fibroblasts) and q < 0.06 in a fourth "
                     "(dendritic cells) under the ranking metric used here; "
                     "alternative metrics and a sample-level analysis place the "
                     "count between six and nine (Fig. S9E).",
    supporting_measurement=S13 + " gives fdr_q 0.0000 (MoMac), 0.0000 (Epithelial), "
                           "0.0000 (Fibroblast), 0.0517 (DC_cells); the 6/8/9 "
                           "alternatives are 16_GSEA_Metric_Sensitivity and "
                           "15_Pseudobulk_Sample_Level",
    note="Adoption lowers the count. The metric clause is what stops the number "
         "moving again: every alternative metric raises it.")

row(claim_id="C03", ledger_id="M03 / R02", source="main text",
    location="clean.docx para 69 s5; edits.py:366-368",
    verbatim="The enrichment was graded rather than uniform: it was the top-ranked "
             "Hallmark gene set in those four populations and in mast cells,",
    printed_value="five populations at rank 1",
    sound_input_value="three (Epithelial, Fibroblast, MoMac)",
    sound_input_table=S13 + " (ttest, post, rank); " + CC + " post_n_rank1",
    axis1_input="moved (5 -> 3)",
    axis2_stability="unstable - 6,5,4,6,5,5,5 across seeds 1/7/13/42/101/777/2026 "
                    "and 6 at the published seed 42 in a fresh run; 6 under all four "
                    "metrics; 3-4 across seeds on sound13",
    gate="plasma cells tie into first place on a 0.0079 NES gap while the seed alone "
         "moves that contrast by 0.044",
    verify_numbers_guard="yes - verify_numbers.py:257 check_eq('populations where "
                         "TNFa/NF-kB is the top-ranked Hallmark set') and :261 the "
                         "named set",
    guard_expected_value="5; {MoMac, DC_cells, Epithelial, Fibroblast, Mast_cells}  "
                         "-> must become 3; {MoMac, Epithelial, Fibroblast}",
    recommendation="restate",
    proposed_wording="The enrichment was graded rather than uniform: it was among "
                     "the three highest-ranked Hallmark gene sets in seven of the "
                     "thirteen populations - monocytes/macrophages, epithelial "
                     "cells, fibroblasts, dendritic cells, endothelial cells, mast "
                     "cells and plasma cells - but ranked in the lower third of the "
                     "sets tested in pericytes and last in B cells, the one "
                     "population with a negative score, which argues against a "
                     "global inflammatory artifact affecting all populations "
                     "equally.",
    supporting_measurement="rank <= 3 after treatment gives the SAME seven cell types "
                           "on the shipped live table and on sound13 "
                           "(live ranks DC 1, Endo 2, Epi 1, Fib 1, Mast 1, MoMac 1, "
                           "Plasma 2; sound13 DC 2, Endo 2, Epi 1, Fib 1, Mast 3, "
                           "MoMac 1, Plasma 3) - measured here from "
                           "nfkb_per_celltype_live.csv and _sound13.csv",
    note="Rank 1 is a coin toss; rank <= 3 is not. This one replacement sentence "
         "also covers C04, C05 and C06. FINDINGS.md section 9 proposed rank <= 3 and "
         "asserted it is stable at every seed and metric; that seed claim is carried "
         "forward from the ledger and was NOT re-measured here (re-measuring means "
         "re-running GSEA). What is measured here is that it is identical on the "
         "shipped and the sound tables.")

row(claim_id="C04", ledger_id="M04 / R03a", source="main text",
    location="clean.docx para 69 s5 (cont.); edits.py:368-369",
    verbatim="but ranked 31st of 49 sets in pericytes",
    printed_value="31st of 49", sound_input_value="30th of 49",
    sound_input_table=S13 + " (ttest, post, Pericyte rank/n_sets); " + CC,
    axis1_input="moved (31 -> 30)",
    axis2_stability="ordinal - rank 31-32 across seeds; 32 / 22 / 12 / 12 across the "
                    "four ranking metrics; 23 (limma-voom) and 3 (DESeq2) at the "
                    "sample level. The DENOMINATOR is a min_size artefact: 33 of 50 "
                    "at min_size 5, 32 of 49 at 10, 32 of 49 at 15, 28 of 43 at 25",
    gate="pericyte post is signal-bearing on both tables",
    verify_numbers_guard="yes - verify_numbers.py:264 check_eq('Pericyte rank among "
                         "Hallmark sets'); :265 ('Hallmark sets tested in Pericyte')",
    guard_expected_value="31; 49  -> 30; 49 on adoption",
    recommendation="restate",
    proposed_wording="see C03 - 'but ranked in the lower third of the sets tested in "
                     "pericytes'. If the ordinal is kept it becomes 30th of 49.",
    supporting_measurement=S13 + " Pericyte post rank 30, n_sets 49; "
                           "ordinal_denominator.csv for the min_size sweep",
    note="The NES itself does not move with min_size (0.959275 at 5, 10, 15 and 25 "
         "to six decimals) - only the denominator does. Dropping the ordinal removes "
         "a number that is fragile on three separate axes.")

row(claim_id="C05", ledger_id="M05 / R03b", source="main text",
    location="clean.docx para 69 s5 (cont.); edits.py:369",
    verbatim="and 37th of 38 in B cells,",
    printed_value="37th of 38", sound_input_value="38th of 38 (last)",
    sound_input_table=S13 + " (ttest, post, B_cells rank/n_sets); " + CC,
    axis1_input="moved (37 -> 38; the set becomes last rather than second-to-last)",
    axis2_stability="ordinal - 37 of 38 at all 7 seeds on live and 38 of 38 at all 7 "
                    "on sound13; 37 or 38 across metrics; DISAGREES at the sample "
                    "level, where B cells post is POSITIVE at rank 11 of 44 "
                    "(limma-voom) and 14 of 44 (DESeq2)",
    gate="CHANGED BY ADOPTION - B cells post is has_signal=yes on live but "
         "has_signal=NO on sound13 (max |NES| 1.418 < 1.5, zero sets at FDR < 0.25). "
         "signal.csv",
    verify_numbers_guard="yes - verify_numbers.py:266 check_eq('B cells rank among "
                         "Hallmark sets'); :267 ('Hallmark sets tested in B cells')",
    guard_expected_value="37; 38  -> 38; 38 on adoption",
    recommendation="restate",
    proposed_wording="see C03 - 'and last in B cells'. That is literally what the "
                     "sound table gives (38th of 38) and it removes a "
                     "second-to-last/last distinction the data does not support.",
    supporting_measurement=S13 + " B_cells post rank 38 of 38, NES -1.225; "
                           "signal.csv B_cells post sound13 has_signal=no",
    note="Correction to stability.csv, which records contrast_has_signal=yes for "
         "this claim. That was scored on the live table. On the table adoption "
         "moves to, the B-cells post contrast carries no gene-set signal at all.")

row(claim_id="C06", ledger_id="M06", source="main text",
    location="clean.docx para 69 s5 (cont.); edits.py:369-371",
    verbatim="the one population with a negative score, which argues against a "
             "global inflammatory artifact affecting all populations equally.",
    printed_value="1 negative (B cells)", sound_input_value="1 negative (B cells)",
    sound_input_table=S13 + " (ttest, post, nes < 0)",
    axis1_input="unchanged",
    axis2_stability="safe within the per-cell Welch t-test - 1 at all 7 seeds, under "
                    "all four metrics and on sound inputs. NOT a property of the "
                    "data: at the sample level 0 of 13 are negative and B cells post "
                    "is positive (+1.32 limma-voom, +1.14 DESeq2)",
    gate="B cells post carries no gene-set signal on sound13 (see C05)",
    verify_numbers_guard="yes - verify_numbers.py:248 check('B cells post NES', "
                         "-0.99, tol 0.01) and :250-253 the 'only one negative' block",
    guard_expected_value="-0.99  -> -1.225 on adoption (the tol=0.01 check FAILS "
                         "unless the expected value moves)",
    recommendation="keep verbatim",
    supporting_measurement=S13 + " B_cells post NES -1.225 (still the only negative)",
    proposed_wording="",
    note="The sentence survives; the number inside its guard does not. Do not "
         "strengthen this into a biological claim about B cells - it is a true "
         "statement about the per-cell t-test and nothing more.")

row(claim_id="C07", ledger_id="M12", source="main text",
    location="clean.docx para 69 s3-s4 (ledger para 68); edits.py:372-375",
    verbatim="Consistent with this, minor-state composition of the CD4+ T, CD8+ T, "
             "NK, B and plasma compartments showed no difference surviving "
             "Benjamini-Hochberg correction, while pathway-level analysis separated "
             "the groups, with interferon-γ programs enriched in responders and "
             "proliferation programs in non-responders (Fig. S10D, S10E).",
    printed_value="direction only - interferon-gamma in R, proliferation in NR",
    sound_input_value="does not survive: 2 of 40 top-5 Hallmark slots shared between "
                      "the shipped tables and the sound recompute; G2-M Checkpoint, "
                      "Mitotic Spindle and E2F Targets - the whole proliferation "
                      "reading - are in the top 5 for neither CD4 nor CD8 in either "
                      "phase",
    sound_input_table="03_Revised_Panels/SUPPLEMENTARY_AUDIT.md fault 1 (measured "
                      "against 12_R1.8_DEG_Recompute/outputs/gsea)",
    axis1_input="moved - almost completely",
    axis2_stability="unstable independently of the input: MAST-recompute vs "
                    "t-test-recompute on the SAME sound matrix also share 0-1 of 5, "
                    "so the lineage top-5 is not stable to the test either",
    gate="not scored - a different prerank (09_R2.2_Adaptive_Immune)",
    verify_numbers_guard="no - verify_numbers.py's R2.2 block reads "
                         "adaptive_state_tests.csv only (:337-341); nothing checks "
                         "adaptive_hallmark_top.csv",
    guard_expected_value="",
    recommendation="drop",
    proposed_wording="Consistent with this, minor-state composition of the CD4+ T, "
                     "CD8+ T, NK, B and plasma compartments showed no difference "
                     "surviving Benjamini-Hochberg correction, and lineage-level "
                     "Hallmark rankings were not stable across differential-"
                     "expression methods, so no pathway-level separation of "
                     "responders from non-responders is claimed here (Fig. S10D, "
                     "S10E).",
    supporting_measurement="SUPPLEMENTARY_AUDIT.md fault 1: 0/5, 0/5, 0/5, 0/5, 1/5, "
                           "0/5, 1/5, 0/5 top-5 overlap across the eight "
                           "lineage x phase cells; and 0-1/5 between MAST and t-test "
                           "on the same sound matrix",
    note="CORRECTION TO THE LEDGER. FINDINGS.md section 8 records Fig. S10D/E as "
         "'the one remaining unaudited GSEA claim' and stability.csv scores M12 "
         "'untestable' on all five axes. It has since been measured, in "
         "SUPPLEMENTARY_AUDIT.md fault 1, and it is the claim that moves most. The "
         "author's decision to change S10E makes this the sentence that follows it.")

row(claim_id="C08", ledger_id="M07 / R04a", source="main text",
    location="clean.docx para 69 s10; edits.py:378-379",
    verbatim="before treatment the gene set is positively enriched in only 6 of 13 "
             "populations",
    printed_value="6 of 13", sound_input_value="5 of 13",
    sound_input_table=S13 + " (ttest, phase=pre, nes > 0); " + CC + " pre_n_positive",
    axis1_input="moved (6 -> 5; MoMac leaves the positive set)",
    axis2_stability="unstable across designs - 6 at all 7 seeds but 6/7/7/7 across "
                    "metrics, 8 (limma-voom) and 9 (DESeq2) at the sample level, 5 on "
                    "sound inputs, 6 under MAST spec A",
    gate="5 of the 13 pre contrasts carry no gene-set signal on either table. On "
         "live, 3 of the 6 positives are signal-bearing; on sound13, 4 of the 5 are",
    verify_numbers_guard="yes - verify_numbers.py:232 check_eq('cell types with "
                         "positive pre NES'); also :332 the S10C 6-and-12 check",
    guard_expected_value="6  -> 5 on adoption (and :332 must move with it)",
    recommendation="restate",
    proposed_wording="before treatment the gene set is positively enriched in five "
                     "of the thirteen populations, and the pre-treatment direction "
                     "is not consistent across analyses (five to nine positive "
                     "depending on the unit of replication and the ranking metric),",
    supporting_measurement=CC + " sound13 pre_n_positive = 5; the 5-to-9 range is "
                           "15_Pseudobulk_Sample_Level and 16_GSEA_Metric_Sensitivity",
    note="Adoption improves the gate here: 4 of the 5 sound-input positives sit in "
         "signal-bearing contrasts against 3 of the 6 on the shipped table.")

row(claim_id="C09", ledger_id="M08", source="main text",
    location="clean.docx para 69 s10 (cont.); edits.py:379",
    verbatim="and reaches q < 0.05 in one,",
    printed_value="one (endothelial cells)",
    sound_input_value="three (endothelial cells, NK cells, CD8+ T cells)",
    sound_input_table=S13 + " (ttest, pre, nes > 0 & fdr_q < 0.05); " + CC +
                      " pre_n_positive_q05, pre_types_q05",
    axis1_input="moved (1 -> 3) - and it moves UPWARD, against the "
                "'largely absent before treatment' framing",
    axis2_stability="floor - 1 at six of seven seeds and 2 at one; 1/4/4/5 across "
                    "metrics; 1 (limma-voom) or 4 (DESeq2) at the sample level; 3 on "
                    "sound inputs; 2 under MAST spec A",
    gate="all three sound-input populations are in signal-bearing pre contrasts "
         "(Endothelial yes, NK yes, TCD8 yes on sound13)",
    verify_numbers_guard="yes - verify_numbers.py:233 check_eq('pre NES with FDR "
                         "< 0.05')",
    guard_expected_value="1  -> 3 on adoption",
    recommendation="restate",
    proposed_wording="and reaches FDR q < 0.05 in three (endothelial, NK and CD8+ T "
                     "cells), none of them the myeloid or epithelial compartments "
                     "the post-treatment signature is carried by,",
    supporting_measurement=S13 + " pre fdr_q: Endothelial 0.0000, NK_cells 0.0474, "
                           "TCD8_cells 0.0156",
    note="This is the row that most needs the author's eye. Adoption makes the "
         "PRE-treatment picture stronger, not weaker, and the three populations "
         "concerned are lymphoid and endothelial - which is why the framing sentence "
         "C13 still holds while this count does not.")

row(claim_id="C10", ledger_id="M09 / R05a", source="main text",
    location="clean.docx para 69 s10 (cont.); edits.py:380",
    verbatim="and monocytes/macrophages are not significant (NES = 1.06, q = 0.47),",
    printed_value="NES 1.06, q 0.47 (table 1.0553 / 0.4721, rank 3 of 45)",
    sound_input_value="NES -1.0011, q 0.8925, rank 21 of 45 - SIGN FLIP",
    sound_input_table=S13 + " (ttest, pre, MoMac); " + CC + " MoMac_pre_nes / _q",
    axis1_input="moved - sign flip, +1.06 -> -1.00, q 0.47 -> 0.89",
    axis2_stability="no signal - on the live table the ES is identical to nine "
                    "decimals (0.362106) at all seven seeds while NES ranges "
                    "0.998-1.098, nominal P 0.087-0.400 and q 0.407-0.592; the "
                    "printed 0.47 is not any seed's answer; zero sets at FDR < 0.25; "
                    "only 4 of 45 sets share the positive ES sign, so the NES "
                    "denominator is estimated from a four-member pool",
    gate="live MoMac pre has_signal=no; sound13 MoMac pre has_signal=weak "
         "(max |NES| 1.415, one set at FDR < 0.25)",
    verify_numbers_guard="yes - verify_numbers.py:243 check('MoMac pre NES', 1.06) "
                         "and :245 check('MoMac pre FDR q', 0.47); also :323 "
                         "check('MoMac pre NES (S10C)', 1.06)",
    guard_expected_value="1.06; 0.47  -> -1.00; 0.89 on adoption, in BOTH places",
    recommendation="restate",
    proposed_wording="and the set was not enriched in monocytes/macrophages before "
                     "treatment (FDR q = 0.89),",
    supporting_measurement=S13 + " MoMac pre NES -1.0011, fdr_q 0.8925, rank 21 of "
                           "45. MAST specification B on the same sound inputs agrees "
                           "on the sign (-1.529, q 0.076): "
                           "14_MAST_Specification/outputs/nfkb_by_specification.csv",
    note="The sentence's CONCLUSION is right on every analysis ever run; its "
         "EVIDENCE is not. Seven runs of this one contrast put it at +1.06, "
         "+1.10, +1.75 (significant), +1.23, +1.68 (significant), -1.00 and -1.53. "
         "The proposed wording quotes only the q, from a contrast that is still only "
         "weakly signal-bearing. If a number is wanted, "
         "'(NES = -1.00, q = 0.89, ranked 21st of 45)' is what the adopted table "
         "gives - but the safest form quotes no NES at all.")

row(claim_id="C11", ledger_id="M10 / R05b", source="main text",
    location="clean.docx para 69 s10 (cont.); edits.py:381-382",
    verbatim="whereas after treatment they show the strongest enrichment of any "
             "population (NES = 2.23, q < 0.001; Fig. S10C).",
    printed_value="NES 2.23, q < 0.001, strongest of 13",
    sound_input_value="NES 2.2283, q 0.0000, rank 1 of 43, strongest of 13",
    sound_input_table=S13 + " (ttest, post, MoMac)",
    axis1_input="unchanged (2.2331 -> 2.2283, within the printed 2 d.p.)",
    axis2_stability="safe - the most robust number in the paper. NES 2.170-2.225 "
                    "over 7 seeds with rank 1 and strongest at all 7; 2.225 / 2.779 "
                    "/ 3.667 / 3.679 under the four metrics, rank 1 and q 0 in all; "
                    "+2.410 rank 1 (limma-voom) and +1.943 rank 2 (DESeq2) at the "
                    "sample level; +1.736 rank 1 under MAST spec A; +2.380 rank 1 "
                    "under MAST spec B",
    gate="MoMac post is the strongest signal-bearing contrast in the study - 18 sets "
         "at FDR < 0.25 on sound13, 19 on live",
    verify_numbers_guard="yes - verify_numbers.py:241 check('MoMac post NES', 2.23); "
                         ":242 ('MoMac post FDR q', 0.0); :324 ('MoMac post NES "
                         "(S10C)'); :328 the idxmax 'strongest' check",
    guard_expected_value="2.23; 0.0  -> unchanged; all four checks still pass",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="every axis agrees; the published ranking metric returns "
                           "the SMALLEST of the four metric values, so the number is "
                           "not being flattered",
    note="Survives adoption untouched, guard included.")

row(claim_id="C12", ledger_id="M11 / R06a", source="main text",
    location="clean.docx para 69 s15; edits.py:387-388",
    verbatim="indicating that epithelium is a minor contributor to the cytokine pool "
             "while showing NF-κB target-gene enrichment itself (NES = 1.86).",
    printed_value="NES 1.86 (table 1.8608)",
    sound_input_value="NES 1.9003, q 0.0000, rank 1 of 42",
    sound_input_table=S13 + " (ttest, post, Epithelial)",
    axis1_input="moved (1.86 -> 1.90)",
    axis2_stability="safe - 1.855-1.884 over 7 seeds with rank 1 at all 7; rank 1 "
                    "under all four metrics; positive and q < 0.001 under both "
                    "sample-level methods; rank 1 under MAST specs A and B",
    gate="epithelial post is signal-bearing on both tables",
    verify_numbers_guard="yes - verify_numbers.py:247 check('Epithelial post NES', "
                         "1.86, tol 0.01)",
    guard_expected_value="1.86  -> 1.90 on adoption (tol=0.01 fails otherwise)",
    recommendation="restate",
    proposed_wording="indicating that epithelium is a minor contributor to the "
                     "cytokine pool while showing NF-κB target-gene enrichment "
                     "itself (NES = 1.90).",
    supporting_measurement=S13 + " Epithelial post NES 1.900318",
    note="The only change is the second decimal. The claim, its rank and its "
         "significance are unaffected; the companion 6-31% cytokine range "
         "(verify_numbers.py:275-279) is not touched by the recompute.")

row(claim_id="C13", ledger_id="(not in ledger)", source="main text",
    location="clean.docx para 69 s10 (opening clause); edits.py:377-378",
    verbatim="Comparison with pre-treatment samples showed that the NF-κB signature "
             "is largely acquired on treatment:",
    printed_value="framing claim, no number",
    sound_input_value="strengthened - MoMac pre goes from +1.06 to -1.00 while MoMac "
                      "post is unmoved at 2.23; positive pre falls 6 -> 5 while "
                      "positive post holds at 12",
    sound_input_table=S13 + "; " + CC,
    axis1_input="unchanged (direction); its supporting counts all move",
    axis2_stability="rests on C08-C11; the contrast between the halves is larger on "
                    "the sound inputs than on the shipped ones",
    gate="8 of 13 post contrasts signal-bearing against 7 of 13 pre on sound13",
    verify_numbers_guard="indirect - through C08, C09, C11 and verify_numbers.py:332",
    guard_expected_value="",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement=CC + ": post 12 positive / pre 5 positive on sound13, "
                           "against 12 / 6 on live",
    note="MISSED BY THE LEDGER as a separate claim. It matters because it is the "
         "sentence the reader takes away, and it is the one part of this passage "
         "that adoption makes stronger. Keep the Abstract's disjoint-timepoint "
         "caveat adjacent to it: the pre and post sample sets share one patient in "
         "the whole cohort, so this is a comparison between different people.")

row(claim_id="C14", ledger_id="(not in ledger)", source="main text",
    location="clean.docx para 69 s17-s19; edits.py:393-396",
    verbatim="Both NF-κB regulons were more active in post-treatment non-responder "
             "monocytes/macrophages (NFKB2, mean AUCell 0.166 versus 0.106, "
             "P = 0.0087; NFKB1, 0.196 versus 0.135, P = 0.017; BACH1, 0.222 versus "
             "0.140, P = 0.017; two-sided Mann-Whitney on sample means, Fig. S11C).",
    printed_value="NFKB2 0.166/0.106 P 0.0087; NFKB1 0.196/0.135 P 0.017; "
                  "BACH1 0.222/0.140 P 0.017",
    sound_input_value="not measured",
    sound_input_table="none - the AUCell matrix is "
                      "Round_5/02_Preparation_for_Panels/SCENIC/aucell_matrix.csv, "
                      "frozen by the author's ruling of 1 September 2026",
    axis1_input="not covered", axis2_stability="not measured",
    gate="not applicable - not a gene-set enrichment",
    verify_numbers_guard="yes - verify_numbers.py:347-352 check the P and both means "
                         "for all three regulons; :353-354 check n_NR = 6, n_R = 5",
    guard_expected_value="0.017/0.196/0.135; 0.0087/0.166/0.106; 0.017/0.222/0.140; "
                         "6; 5  -> unchanged",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="07_R1.8_NFkB_Specificity/outputs/nfkb_regulon_activity.csv",
    note="MISSED BY THE LEDGER, which covers GSEA values only. This is the "
         "orthogonal, motif-anchored evidence the NF-kB claim actually rests on when "
         "the gene-set numbers move, so it belongs in the comparison. Caveat for the "
         "author, not a defect found here: the deposited SCENIC run is a 12k-cell "
         "subset (see the project memory note), and n = 6 versus 5 samples.")

row(claim_id="C15", ledger_id="(not in ledger)", source="main text",
    location="clean.docx para 69 s21; edits.py:396-399",
    verbatim="The pathway's own negative-feedback targets (NFKBIA, TNFAIP3, NFKB2, "
             "RELB, BIRC3, TRAF1), which are transcribed only after nuclear "
             "translocation of the complex, were elevated in the same direction.",
    printed_value="direction only",
    sound_input_value="not measured", sound_input_table="none",
    axis1_input="not covered", axis2_stability="not measured",
    gate="not applicable",
    verify_numbers_guard="no", guard_expected_value="",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="07_R1.8_NFkB_Specificity/outputs/nfkb_feedback_targets.csv",
    note="MISSED BY THE LEDGER. A direction claim with no number; unaffected by "
         "adoption. Unguarded by verify_numbers.py.")

row(claim_id="C16", ledger_id="(not in ledger)", source="main text",
    location="clean.docx para 69 s22; edits.py:399-404",
    verbatim="We emphasize that these analyses measure the transcriptional output "
             "and inferred activity of NF-κB and do not establish biochemical "
             "pathway activation; no protein-level measurements of nuclear p65, "
             "phospho-p65, or IκB degradation were available for this cohort.",
    printed_value="scope limitation, no number",
    sound_input_value="not applicable", sound_input_table="none",
    axis1_input="not covered", axis2_stability="not measured",
    gate="not applicable", verify_numbers_guard="no", guard_expected_value="",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="",
    note="MISSED BY THE LEDGER. It is the sentence that makes every other NF-kB "
         "claim defensible and it must survive any restatement above it unchanged.")

row(claim_id="C17", ledger_id="(not in ledger)", source="main text",
    location="clean.docx para 65 s10",
    verbatim="Similarly, NFKB1 regulon activity was elevated in non-responder "
             "IL-1β+ macrophages (two-sided exact permutation test, P = 0.038; "
             "Fig. 5E, Table S6),",
    printed_value="P = 0.038", sound_input_value="not measured",
    sound_input_table="none - pySCENIC AUCell, frozen",
    axis1_input="not covered", axis2_stability="not measured",
    gate="not applicable",
    verify_numbers_guard="yes - verify_numbers.py:125 check('NFKB1 regulon P', "
                         "sweep_p('NFKB1 regulon activity'), 0.038)",
    guard_expected_value="0.038  -> unchanged",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="the sweep table read by verify_numbers.py:125",
    note="MISSED BY THE LEDGER. A second, earlier NF-kB claim in the Results, "
         "guarded but never entered in ledger.csv.")

row(claim_id="C26", ledger_id="M22 / F01 (the sentence, not the panel)",
    source="main text", location="clean.docx para 63 s4",
    verbatim="TNFα signaling via NF-κB, inflammatory response, interferon gamma "
             "response, and IL-6/JAK/STAT3 signaling were among the most "
             "significantly enriched pathways (Fig. 5A, 5B).",
    printed_value="direction only; Fig. 5A plots TNF-alpha/NF-kB at +2.140, the "
                  "largest positive of nine bars",
    sound_input_value="not covered - Fig. 5A is a MAST prerank of a Round_5 "
                      "differential expression, a third pipeline",
    sound_input_table="Round_5/02_Preparation_for_Panels/GSEA/"
                      "MoMac_mast_prerank_gsea.csv (frozen; reproduces the printed "
                      "panel to 0.0001 NES)",
    axis1_input="not covered", axis2_stability="not measured - Fig. 5A was not "
                "seed-swept; PROVENANCE.csv records reproduces_published = yes",
    gate="not applicable",
    verify_numbers_guard="no", guard_expected_value="",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="PROVENANCE.csv figure 5 panel A, verdict reproduced "
                           "2026-09-01; VERIFICATION_ADDENDUM.md addendum 3",
    note="MISSED BY THE LEDGER as a sentence (ledger.csv carries the panel, M22/F01, "
         "but not the text that quotes it). It is the one Results sentence still "
         "resting on MAST after MAST was removed from the Methods - see the "
         "'held in reserve' section.")

# ==================================================== Results, downstream cascade
row(claim_id="C18", ledger_id="M13 + M14", source="main text",
    location="clean.docx para 71 s16 (ledger para 70)",
    verbatim="Epithelial cells showed enrichment for inflammatory response "
             "(NES = 1.48, P = 0.02) and hypoxia (NES = 1.56, P = 0.006), both "
             "established NF-κB-regulated programs that promote tumor cell survival "
             "and plasticity.",
    printed_value="1.48 / P 0.02; 1.56 / P 0.006",
    sound_input_value="not covered - a third pipeline, whose panels read .raw and are "
                      "therefore not affected by the double normalisation",
    sound_input_table="03_Revised_Panels/Main_Figures/05_Figure_5/05_GSEA_Summary/"
                      "gsea_data/gsea_combined_5types.csv",
    axis1_input="not covered",
    axis2_stability="NES reproduces (1.4780 and 1.5601 at seed 42; 1.452-1.478 and "
                    "1.552-1.601 over 7 seeds). The printed P values are NOT the "
                    "table's: table 0.0115 and 0.0083 against printed 0.02 and "
                    "0.006. FDR q is 0.090-0.116 and 0.077-0.135 - neither survives "
                    "FDR < 0.05",
    gate="all 50 Hallmark sets tested; ranked lists 56,000-57,000 genes",
    verify_numbers_guard="no - verify_numbers.py has no check for any para-71 value",
    guard_expected_value="",
    recommendation="keep verbatim (flagged)",
    proposed_wording="",
    supporting_measurement=P70 + " rows Epithelial / Inflammatory Response and "
                          "Epithelial / Hypoxia",
    note="CLAUDE.md rule 1. The printed figure is the ground truth and these are "
         "pre-existing submitted numbers; a P-value mismatch in the code's table is "
         "not grounds to change the text. FOR THE AUTHOR: check the two P values "
         "against printed Fig. 5N before the next submission, and consider whether "
         "quoting a nominal P for a claim at FDR q ~ 0.09-0.13 is the form wanted.")

row(claim_id="C19", ledger_id="M15", source="main text",
    location="clean.docx para 71 s17 (ledger para 70)",
    verbatim="Fibroblasts exhibited enrichment for inflammatory response "
             "(NES = 1.30, P = 0.02), consistent with NF-κB-driven stromal "
             "activation.",
    printed_value="NES 1.30, P 0.02",
    sound_input_value="not covered",
    sound_input_table="gsea_combined_5types.csv gives NES 1.279612, NOM_pval "
                      "0.019694, FDR_qval 0.931735",
    axis1_input="not covered",
    axis2_stability="NOT REPRODUCIBLE from any table on disk - 1.2796 at seed 42, "
                    "1.2796-1.2980 over 7 seeds; the band does not reach 1.30, so "
                    "the 0.020 gap is not a seed effect. Every CSV with an NES "
                    "column under 03_Revised_Panels/Main_Figures, 02_New_Analyses "
                    "and Round_5/02_Preparation_for_Panels was scanned; none gives "
                    "1.30. FDR q 0.676-0.932",
    gate="all 50 Hallmark sets tested",
    verify_numbers_guard="no", guard_expected_value="",
    recommendation="keep verbatim (flagged)", proposed_wording="",
    supporting_measurement=P70 + " row Fibroblast / Inflammatory Response, "
                          "abs_diff_from_printed 0.0204",
    note="CLAUDE.md rule 1 - recorded as not reproducible, no correction proposed. "
         "Mechanism not found. The far more consequential fact for the author is "
         "the FDR: q = 0.68 to 0.93 across seeds, so this claim rests on a nominal "
         "P alone.")

row(claim_id="C20", ledger_id="M16-M19", source="main text",
    location="clean.docx para 71 s18 (ledger para 70); edits.py:875-878",
    verbatim="Monocytes/macrophages themselves showed the broadest downstream "
             "activation, including inflammatory response (NES = 1.67), EMT "
             "(NES = 1.66), a program whose intermediate states are increasingly "
             "resolved at single-cell level (38,39), hypoxia (NES = 1.53), and "
             "angiogenesis (NES = 1.50; all P ≤ 0.05), reflecting both autocrine "
             "NF-κB amplification and their role as the primary source of "
             "pro-inflammatory signaling.",
    printed_value="1.67 / 1.66 / 1.53 / 1.50, all P <= 0.05",
    sound_input_value="not covered",
    sound_input_table="03_Revised_Panels/Main_Figures/05_Figure_5/05_GSEA_Summary/"
                      "gsea_data/gsea_momac.csv: 1.672274, 1.696945, 1.531753, "
                      "1.500481; NOM_pval 0.0000, 0.0000, 0.004739, 0.023622",
    axis1_input="not covered",
    axis2_stability="three of the four reproduce (1.6723, 1.5318, 1.5005). EMT does "
                    "NOT: the table gives 1.6969 against a printed 1.66, and the "
                    "seed band is 1.6529-1.6969, so the 0.037 gap is not a seed "
                    "effect. 'all P <= 0.05' holds on the table's nominal P. These "
                    "four are the only para-71 values that survive FDR < 0.10",
    gate="all 50 Hallmark sets tested",
    verify_numbers_guard="no", guard_expected_value="",
    recommendation="keep verbatim (flagged)", proposed_wording="",
    supporting_measurement=P70 + " rows MoMac / Inflammatory Response, EMT, Hypoxia, "
                          "Angiogenesis",
    note="CLAUDE.md rule 1 - the EMT mismatch is recorded, not acted on. Note that "
         "MoMac post TNFa/NF-kB has three different printed values in three panels: "
         "+2.140 (Fig. 5A, MAST prerank), +2.233 (Fig. 5H / S9E / Results) and "
         "+1.894 (Fig. 5N, this table). All positive, all at or near rank 1; the "
         "decimals differ because they are three analyses, which the paper does not "
         "tell the reader.")

row(claim_id="C21", ledger_id="M20", source="main text",
    location="clean.docx para 71 s14 (ledger para 70)",
    verbatim="Systematic GSEA across monocytes/macrophages, epithelial cells, "
             "fibroblasts, and dendritic cells revealed that NF-κB-driven pathways "
             "were enriched in all four cell types in post-treatment non-responders, "
             "with cell-type-specific functional manifestations (Fig. 5N).",
    printed_value="all four positive",
    sound_input_value="not covered",
    sound_input_table="gsea_momac.csv and gsea_combined_5types.csv",
    axis1_input="not covered",
    axis2_stability="supported - TNF-alpha/NF-kB NES is +1.894 (MoMac), +1.748 "
                    "(Epithelial), +1.468 (Fibroblast), +1.522 (DC), and "
                    "inflammatory response is positive in all four",
    gate="all 50 sets tested in each",
    verify_numbers_guard="no", guard_expected_value="",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="the two gsea_data tables, read directly; three of the "
                           "four reach FDR q < 0.05 on TNF-alpha/NF-kB (DC q = 0.073)",
    note="Direction claim, supported as written.")

row(claim_id="C22", ledger_id="(not in ledger)", source="main text",
    location="clean.docx para 71 s19; edits.py:882-885",
    verbatim="These convergent yet cell-type-specific pathway enrichments show that "
             "an NF-κB transcriptional signature is detectable across tumor, "
             "stromal, and immune compartments, accompanying coordinated TME "
             "reprogramming that extends well beyond checkpoint ligand upregulation.",
    printed_value="direction only",
    sound_input_value="not covered", sound_input_table="",
    axis1_input="not covered", axis2_stability="supported by C21",
    gate="", verify_numbers_guard="no", guard_expected_value="",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="",
    note="MISSED BY THE LEDGER. Already softened this round from "
         "'macrophage-derived NF-kB signaling propagates ... driving' to "
         "'a signature is detectable ... accompanying'. Nothing in the adoption "
         "touches it.")

# ============================================ Abstract, significance, Discussion
row(claim_id="C23", ledger_id="(not in ledger)", source="main text",
    location="clean.docx para 17 s7 (Abstract)",
    verbatim="These cells show elevated NFKB1 and NFKB2 regulon activity and "
             "accompany a coordinated NF-κB transcriptional signature across tumor "
             "and microenvironment compartments, together with PD-L1 upregulation, "
             "epithelial-mesenchymal transition, and chronic inflammation.",
    printed_value="direction only", sound_input_value="unchanged",
    sound_input_table=CC, axis1_input="unchanged",
    axis2_stability="safe - 12 of 13 positive post is stable on every axis",
    gate="", verify_numbers_guard="no", guard_expected_value="",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="C01, C11, C14",
    note="MISSED BY THE LEDGER. The Abstract quotes no NF-kB number, which is why it "
         "survives adoption untouched.")

row(claim_id="C24", ledger_id="(not in ledger)", source="main text",
    location="clean.docx para 20 s2 (Significance / summary)",
    verbatim="CEACAM5/6⁺ cancer cells are enriched in pre-treatment non-responders "
             "and occupy immune-poor spatial regions, whereas IL-1β⁺ macrophages "
             "are enriched in post-treatment non-responders and accompany an NF-κB "
             "transcriptional signature across the tumor microenvironment.",
    printed_value="direction only", sound_input_value="unchanged",
    sound_input_table=CC, axis1_input="unchanged", axis2_stability="safe",
    gate="", verify_numbers_guard="no", guard_expected_value="",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="C01, C11",
    note="MISSED BY THE LEDGER.")

row(claim_id="C25", ledger_id="(not in ledger)", source="main text",
    location="clean.docx para 76 s7 (Discussion)",
    verbatim="Although the CEACAM5/6 spatial niche includes local myeloid "
             "recruitment, the NF-κB-associated inflammatory program is broader and "
             "is detectable only after treatment.",
    printed_value="'detectable only after treatment'",
    sound_input_value="contradicted in three populations: on sound inputs the set "
                      "reaches FDR q < 0.05 BEFORE treatment in endothelial, NK and "
                      "CD8+ T cells",
    sound_input_table=S13 + " (ttest, pre, fdr_q)",
    axis1_input="moved - the pre-treatment count at q < 0.05 rises from 1 to 3",
    axis2_stability="floor, inherited from C09; the count is 1/4/4/5 across metrics "
                    "and 1 or 4 at the sample level, so 'only' was already the "
                    "weakest reading available",
    gate="all three are signal-bearing pre contrasts on sound13",
    verify_numbers_guard="no", guard_expected_value="",
    recommendation="restate",
    proposed_wording="Although the CEACAM5/6 spatial niche includes local myeloid "
                     "recruitment, the NF-κB-associated inflammatory program is "
                     "broader and is far more pronounced after treatment, where it "
                     "extends to the myeloid and epithelial compartments that carry "
                     "it.",
    supporting_measurement=S13 + " pre fdr_q: Endothelial 0.0000, NK 0.0474, "
                           "TCD8 0.0156; post: MoMac 0.0000, Epithelial 0.0000, "
                           "Fibroblast 0.0000",
    note="MISSED BY THE LEDGER, and it is the one absolute word in the paper that "
         "adoption breaks. 'Only after treatment' becomes false the moment the "
         "pre-treatment q < 0.05 count goes from one to three.")

# ==================================================================== Methods
row(claim_id="C27", ledger_id="(part of M21)", source="main text - Methods",
    location="clean.docx para 118 s3 (ledger para 117)",
    verbatim="Enrichment was assessed against MSigDB Hallmark 2020 gene sets with "
             "1,000 permutations, minimum gene set size of 15, and maximum of 500.",
    printed_value="min 15, max 500, 1,000 permutations",
    sound_input_value="agrees with the code",
    sound_input_table="12_R1.8_DEG_Recompute/scripts/recompute_deg.py:351 "
                      "gp.prerank(min_size=15, max_size=500, seed=42)",
    axis1_input="unchanged", axis2_stability="correct as printed",
    gate="", verify_numbers_guard="no - verify_numbers.py checks values, not settings",
    guard_expected_value="",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="VERIFICATION_ADDENDUM.md addendum 1; the module's "
                           "positive control reproduces the deposited tables "
                           "bit-for-bit at min_size 15 and not at 5",
    note="CORRECTION TO THE LEDGER, already recorded in VERIFICATION_ADDENDUM.md. "
         "FINDINGS.md section 5.3 says 'the printed Methods declare min_size = 5'. "
         "That is true of para 121 and FALSE of this paragraph, which declares 15 "
         "and is right.")

row(claim_id="C28", ledger_id="M21", source="main text - Methods",
    location="clean.docx para 121 s1 and s5 (ledger para 120)",
    verbatim="For NF-κB pathway enrichment analysis across cell types (Fig. 5H, I) "
             "... GSEA was then performed using GSEApy prerank against MSigDB "
             "Hallmark 2020 gene sets (min_size = 5, max_size = 500; 1,000 "
             "permutations; seed = 42).",
    printed_value="min_size = 5, scoped to 'Fig. 5H, I'",
    sound_input_value="the paragraph's three settings - the 5,000-cell subsample, "
                      "the absence of a detection filter and min_size = 5 - are all "
                      "Fig. 5I's. Fig. 5H is made by recompute_deg.py: no subsample, "
                      "10% detection filter, min_size = 15",
    sound_input_table="run_gsea_momac.py:55 and run_gsea_5celltypes.py:65 "
                      "(min_size 5); recompute_deg.py:121 (MIN_PCT 0.10) and :351 "
                      "(min_size 15)",
    axis1_input="not a number - unchanged",
    axis2_stability="wrong in scope, not in settings. The two printed ordinals (C04, "
                    "C05) really are min_size = 15 quantities; at min_size = 5 they "
                    "would be 33 of 50 and 42 of 43",
    gate="", verify_numbers_guard="no", guard_expected_value="",
    recommendation="restate",
    proposed_wording="For the cell-type GSEA shown in Fig. 5I and Fig. 5N, "
                     "differential expression was computed within each cell type "
                     "after random subsampling of populations exceeding 5,000 cells "
                     "(n = 5,000; random_state = 42), and GSEApy prerank was run "
                     "against MSigDB Hallmark 2020 gene sets (min_size = 5, "
                     "max_size = 500; 1,000 permutations; seed = 42). The "
                     "cell-type NF-κB comparison shown in Fig. 5H and Fig. S9E and "
                     "S10C uses the differential expression described above, with a "
                     "10% detection filter, no subsampling, and min_size = 15.",
    supporting_measurement="VERIFICATION_ADDENDUM.md addendum 1, which checked both "
                           "Methods paragraphs and both pipelines against source",
    note="A scoping fix, not a settings fix. No number moves. It matters because "
         "min_size sets the DENOMINATOR of C04 and C05 while leaving the NES "
         "untouched to six decimals, so a reader taking min_size = 5 at face value "
         "cannot reconstruct either ordinal.")

row(claim_id="C29", ledger_id="(not in ledger)", source="main text - Methods",
    location="clean.docx para 115 area; edits.py:1113-1126 (ROUND 10, point "
             "'R1.8, ED.4')",
    verbatim="Differential expression between responders and non-responders was "
             "performed with a Welch t-test across cells with overestimated variance "
             "(scanpy rank_genes_groups, method = 't-test_overestim_var') ... For "
             "cell populations exceeding 5,000 cells, random subsampling was "
             "performed (n = 5,000; random_state = 42).",
    printed_value="Welch t-test; 10% detection filter; BH; 5,000-cell subsample",
    sound_input_value="the first three describe recompute_deg.py correctly; the "
                      "subsampling sentence does not - recompute_deg.py has no "
                      "subsample",
    sound_input_table="recompute_deg.py:121 (MIN_PCT = 0.10) and the absence of any "
                      "subsample; run_gsea_momac.py:38-40 and "
                      "run_gsea_5celltypes.py:37-49 (where the subsample lives)",
    axis1_input="not a number - unchanged",
    axis2_stability="the same conflation as C28: one paragraph describing two "
                    "pipelines as one",
    gate="", verify_numbers_guard="no", guard_expected_value="",
    recommendation="restate",
    proposed_wording="Differential expression between responders and non-responders "
                     "was performed with a Welch t-test across cells with "
                     "overestimated variance (scanpy rank_genes_groups, method = "
                     "'t-test_overestim_var'), applied separately to pre-treatment "
                     "and post-treatment stomach samples. Genes detected in fewer "
                     "than 10% of the cells entering a comparison were excluded, and "
                     "P values were adjusted using the Benjamini-Hochberg method. "
                     "(delete the subsampling sentence here; it belongs in the "
                     "Fig. 5I / 5N paragraph - see C28.)",
    supporting_measurement="VERIFICATION_ADDENDUM.md addendum 2, section 'A second "
                           "inaccuracy in the same replacement paragraph'; the note "
                           "on edits.py:1152-1157 already flags it as open for the "
                           "author",
    note="MISSED BY THE LEDGER as a claim. edits.py's own note records it as still "
         "open. Fixing C28 and C29 together removes the conflation once.")

# ==================================================================== Legends
def legend(cid, lid, loc, text, rec, note, guard="no", axis1="unchanged",
           axis2="inherits", proposed="", support=""):
    row(claim_id=cid, ledger_id=lid, source="figure legend", location=loc,
        verbatim=text, printed_value="panel description", sound_input_value="",
        sound_input_table="", axis1_input=axis1, axis2_stability=axis2,
        gate="", verify_numbers_guard=guard, guard_expected_value="",
        recommendation=rec, proposed_wording=proposed,
        supporting_measurement=support, note=note)


legend("C30", "M22 / F01", "clean.docx para 242 s1 (Figure 5 legend)",
       "(A) Top 9 Hallmark pathways enriched in Monocytes/Macrophages "
       "(post-treatment NR versus R), ranked by normalized enrichment score.",
       "keep verbatim",
       "PROVENANCE.csv: figure 5, printed panel A, build dir 05_A, "
       "reproduces_published = yes, reproduced 2026-09-01 to 0.0001 NES against "
       "Round_5/.../MoMac_mast_prerank_gsea.csv. Frozen; a Version A ship.",
       axis1="not covered", axis2="reproduces the printed panel",
       support="PROVENANCE.csv; VERIFICATION_ADDENDUM.md addendum 3")

legend("C31", "M23 / F02", "clean.docx para 242 s12 (Figure 5 legend)",
       "(H) Radar plot comparing NF-κB pathway NES between pre- and post-treatment "
       "across cell types.",
       "keep verbatim",
       "PROVENANCE.csv: figure 5, printed panel H, build dir 05_F (the directory "
       "letter is NOT the panel letter). The panel reads "
       "nfkb_per_celltype.csv, so adoption redraws it; the legend sentence itself is "
       "unaffected. Author's ruling of 1 September 2026: the pre-treatment ring "
       "stays and is not greyed or annotated for the signal gate.",
       axis1="the plotted values move; the legend does not",
       axis2="the post half is sound; the pre half plots five contrasts with no "
             "gene-set signal",
       support="signal.csv; VERIFICATION_ADDENDUM.md, author's ruling")

legend("C32", "F05", "clean.docx para 242 s13 (Figure 5 legend)",
       "(I) GSEA enrichment curves for TNF-α Signaling via NF-κB in fibroblasts and "
       "epithelial cells (post-treatment R vs. NR).",
       "keep verbatim",
       "Fig. 5I is the panel-script pipeline (min_size = 5, 5,000-cell subsample, no "
       "detection filter), not the recompute, so adoption does not touch it. "
       "PROVENANCE.csv records reproduces_published = yes.",
       axis1="not covered", axis2="separate pipeline")

legend("C33", "F06", "clean.docx para 242 s22 (Figure 5 legend)",
       "(N) Hallmark pathway enrichment (GSEA) in monocytes/macrophages, epithelial "
       "cells, fibroblasts, and DC cells (post-treatment NR versus R).",
       "keep verbatim",
       "The source of C18-C21. Not touched by adoption.",
       axis1="not covered", axis2="separate pipeline")

legend("C34", "F03", "clean.docx para 257 s7 (Figure S9 legend)",
       "(E) Hallmark TNFα signaling via NF-κB normalized enrichment score, FDR and "
       "rank among all Hallmark gene sets tested, for each of the 13 cell types after "
       "treatment.",
       "keep verbatim",
       "The panel that carries C01-C06. Adoption redraws every value in it; the "
       "legend sentence describes the panel correctly either way.",
       axis1="the plotted values move; the legend does not",
       axis2="inherits C01-C06")

legend("C35", "F04", "clean.docx para 258 s5 (Figure S10 legend)",
       "(C) NF-κB normalized enrichment score before and after treatment for each "
       "cell type.",
       "keep verbatim",
       "The panel that carries C08-C11. Adoption redraws it. verify_numbers.py:321-333 "
       "checks four values on its source table.",
       guard="indirect - verify_numbers.py:323-333",
       axis1="the plotted values move; the legend does not",
       axis2="inherits C08-C11")

legend("C36", "(not in ledger)", "clean.docx para 258 s7 (Figure S10 legend)",
       "(E) Hallmark pathways most enriched in each adaptive lineage at each "
       "timepoint.",
       "keep verbatim",
       "MISSED BY THE LEDGER. The legend is a description, not a claim, and stays "
       "true after the repoint the author has ruled for - it is the SENTENCE C07 "
       "that does not.",
       axis1="the plotted values move almost completely (2 of 40 top-5 slots "
             "survive); the legend does not",
       axis2="see C07")

legend("C37", "(not in ledger)", "clean.docx para 259 s5 (Figure S11 legend)",
       "(C) NFKB1, NFKB2 and BACH1 regulon activity from motif-anchored pySCENIC, "
       "with the NF-κB negative-feedback target genes.",
       "keep verbatim",
       "MISSED BY THE LEDGER. The panel that carries C14 and C15. pySCENIC output is "
       "frozen by the author's ruling of 1 September 2026, so adoption does not "
       "touch it.",
       guard="indirect - verify_numbers.py:345-354",
       axis1="not covered", axis2="not measured")

# ============================================================== Response letter
row(claim_id="C38", ledger_id="R01 + R02 + R03", source="response letter",
    location="Response_to_Reviewers_CIR260753ET_v3_clean.docx para 85 s3 "
             "(ledger para 84)",
    verbatim="The Hallmark TNFα signaling via NF-κB program is enriched toward "
             "post-treatment non-response in 12 of 13 populations, with four reaching "
             "FDR q < 0.05, and it is the top-ranked Hallmark set in five populations "
             "while ranking 31st of 49 in pericytes and 37th of 38 in B cells, so the "
             "enrichment is graded rather than uniform; the manuscript therefore no "
             "longer states that “every population” shows NF-κB activation "
             "(Fig. S9E).",
    printed_value="12 of 13; four; five; 31 of 49; 37 of 38",
    sound_input_value="12 of 13; three; three; 30 of 49; 38 of 38",
    sound_input_table=S13 + "; " + CC + " row sound13",
    axis1_input="moved in four of its five numbers (only 12 of 13 holds)",
    axis2_stability="inherits C01 (safe), C02 (floor), C03 (unstable), C04 and C05 "
                    "(ordinals)",
    gate="see C01-C05",
    verify_numbers_guard="yes - the same checks as C01-C05 "
                         "(verify_numbers.py:228-233, 257-267)",
    guard_expected_value="see C01-C05",
    recommendation="restate",
    proposed_wording="The Hallmark TNFα signaling via NF-κB program is enriched "
                     "toward post-treatment non-response in 12 of 13 populations, "
                     "three of them at FDR q < 0.05 and a fourth at q = 0.05, and it "
                     "is among the three highest-ranked Hallmark sets in seven "
                     "populations while ranking in the lower third in pericytes and "
                     "last in B cells, so the enrichment is graded rather than "
                     "uniform; the manuscript therefore no longer states that "
                     "“every population” shows NF-κB activation (Fig. S9E).",
    supporting_measurement="the same tables as C01-C05; the letter must move with the "
                           "main text or verify_numbers.py's cross-check on the two "
                           "will diverge",
    note="The letter and the main text must be edited in the same pass. The letter "
         "is edited only through apply_consistency_fixes.py.")

row(claim_id="C39", ledger_id="R06", source="response letter",
    location="Response letter para 87 s2 (ledger para 86)",
    verbatim="These cytokines are detectable in post-treatment epithelial cells but "
             "are substantially lower than in the myeloid compartment (6-31% of the "
             "corresponding monocyte/macrophage levels), while epithelial cells "
             "themselves carry a clear NF-κB target-gene signature (NES = 1.86, "
             "FDR q < 0.001, the top-ranked Hallmark set in that compartment; "
             "Fig. S9E).",
    printed_value="6-31%; NES 1.86; q < 0.001; rank 1",
    sound_input_value="6-31% unchanged; NES 1.9003; q 0.0000; rank 1 of 42",
    sound_input_table=S13 + " (ttest, post, Epithelial)",
    axis1_input="moved in one number (1.86 -> 1.90); the rank, the q and the "
                "cytokine range are unchanged",
    axis2_stability="safe - see C12",
    gate="epithelial post is signal-bearing on both tables",
    verify_numbers_guard="yes - verify_numbers.py:247 (the NES) and :275-279 "
                         "(the 6-31% range)",
    guard_expected_value="1.86 -> 1.90; the 5.5-6.5 / 30.0-32.0 range is unchanged",
    recommendation="restate",
    proposed_wording="... while epithelial cells themselves carry a clear NF-κB "
                     "target-gene signature (NES = 1.90, FDR q < 0.001, the "
                     "top-ranked Hallmark set in that compartment; Fig. S9E).",
    supporting_measurement=S13 + " Epithelial post NES 1.900318, fdr_q 0.0, rank 1",
    note="A one-decimal edit. The rank-1 half of the claim is one of the most stable "
         "statements in the paper - rank 1 at all 7 seeds, under all four metrics, "
         "on the sound inputs and under both corrected MAST specifications.")

row(claim_id="C40", ledger_id="R04 + R05", source="response letter",
    location="Response letter para 99 s4 (ledger para 98)",
    verbatim="Consistently, the Hallmark TNFα/NF-κB set is positively enriched in 6 "
             "of 13 populations before treatment and in 12 of 13 after, and "
             "monocytes/macrophages move from not significant before treatment "
             "(NES = 1.06, q = 0.47) to the strongest enrichment of any population "
             "after it (NES = 2.23, q < 0.001) (Fig. S10A-C).",
    printed_value="6 of 13; 12 of 13; 1.06 / 0.47; 2.23 / q < 0.001",
    sound_input_value="5 of 13; 12 of 13; -1.00 / 0.89; 2.2283 / q 0.0000",
    sound_input_table=S13 + "; " + CC + " row sound13",
    axis1_input="moved in two of its four numbers (12 of 13 and 2.23 hold)",
    axis2_stability="inherits C08 (unstable), C01 (safe), C10 (no signal), C11 (safe)",
    gate="see C08-C11",
    verify_numbers_guard="yes - verify_numbers.py:321-333 checks all four against "
                         "nfkb_pre_vs_post.csv",
    guard_expected_value="6 and 12 (:332); 1.06 (:323); 2.23 (:324)  -> 5 and 12; "
                         "-1.00; 2.23",
    recommendation="restate",
    proposed_wording="Consistently, the Hallmark TNFα/NF-κB set is positively "
                     "enriched in five of the thirteen populations before treatment "
                     "and in twelve of thirteen after, and monocytes/macrophages "
                     "move from no enrichment before treatment (FDR q = 0.89) to the "
                     "strongest enrichment of any population after it (NES = 2.23, "
                     "q < 0.001) (Fig. S10A-C).",
    supporting_measurement="see C08, C10, C11",
    note="The letter's own next sentence - that the two timepoints are predominantly "
         "unpaired - is the right caveat and should stay exactly where it is.")

row(claim_id="C41", ledger_id="R07", source="response letter",
    location="Response letter para 108 s4 (ledger para 107)",
    verbatim="The most informative additional result is in adaptive T cells: "
             "post-treatment responder CD4+ and CD8+ T cells show stronger "
             "interferon-γ/inflammatory programs, whereas non-responder T cells are "
             "enriched for proliferation-associated programs.",
    printed_value="direction only",
    sound_input_value="does not survive - see C07",
    sound_input_table=AUDIT,
    axis1_input="moved - almost completely",
    axis2_stability="unstable to the test as well as to the input (MAST vs t-test on "
                    "the same sound matrix share 0-1 of 5 top-5 slots)",
    gate="not scored", verify_numbers_guard="no", guard_expected_value="",
    recommendation="drop",
    proposed_wording="The additional lineage-level Hallmark results are provided as a "
                     "resource; we do not draw a pathway-level distinction between "
                     "responders and non-responders from them, because the "
                     "lineage-level rankings are not stable across "
                     "differential-expression methods.",
    supporting_measurement=AUDIT,
    note="This sentence was written to a reviewer as 'the most informative "
         "additional result'. It is the least stable claim in the set. It must go "
         "with C07.")

row(claim_id="C42", ledger_id="(not in ledger)", source="response letter",
    location="Response letter para 85 s5",
    verbatim="As an orthogonal transcription-factor-specific analysis, NFKB1 and "
             "NFKB2 pySCENIC regulon activity is higher in post-treatment "
             "non-responder monocyte/macrophage cells.",
    printed_value="direction only",
    sound_input_value="not measured",
    sound_input_table="none - pySCENIC AUCell, frozen",
    axis1_input="not covered", axis2_stability="not measured",
    gate="not applicable",
    verify_numbers_guard="yes - verify_numbers.py:345-354 (the values behind it)",
    guard_expected_value="unchanged",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="07_R1.8_NFkB_Specificity/outputs/nfkb_regulon_activity.csv",
    note="MISSED BY THE LEDGER. This is the letter's answer to the reviewer's "
         "'enrichment is not activation' point, and it does not depend on any number "
         "adoption moves - which makes it more valuable after adoption, not less.")

row(claim_id="C43", ledger_id="(not in ledger)", source="response letter",
    location="Response letter para 85 s6",
    verbatim="These results are described as NF-κB transcriptional/regulon activity, "
             "not biochemical activation.",
    printed_value="scope limitation", sound_input_value="not applicable",
    sound_input_table="", axis1_input="not covered", axis2_stability="not measured",
    gate="", verify_numbers_guard="no", guard_expected_value="",
    recommendation="keep verbatim", proposed_wording="",
    supporting_measurement="",
    note="MISSED BY THE LEDGER. Pairs with C16 in the main text.")

legend("C44", "(not in ledger)", "Response letter para 89 s3 (S9E legend, in the "
       "letter)",
       "Hallmark TNFα signaling via NF-κB across the 13 cell types after treatment, "
       "with the FDR and the rank of the set among all Hallmark sets tested in that "
       "cell type.",
       "keep verbatim",
       "MISSED BY THE LEDGER. Duplicate of C34, carried in the letter; must not "
       "diverge from it.")

legend("C45", "(not in ledger)", "Response letter para 105 s3 (S10C legend, in the "
       "letter)",
       "NF-κB normalized enrichment score before and after treatment across cell "
       "types.",
       "keep verbatim",
       "MISSED BY THE LEDGER. Duplicate of C35, carried in the letter.")


def main():
    CSV_OUT.parent.mkdir(parents=True, exist_ok=True)
    with CSV_OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=F)
        w.writeheader()
        w.writerows(R)
    from collections import Counter
    c = Counter(r["recommendation"] for r in R)
    print("rows:", len(R))
    for k, v in sorted(c.items()):
        print(f"  {k}: {v}")
    print("ledger-missed:", sum(1 for r in R if "not in ledger" in r["ledger_id"]))
    print("axis1 not covered:", sum(1 for r in R
                                    if r["axis1_input"].startswith("not covered")))
    print("axis1 moved:", sum(1 for r in R if r["axis1_input"].startswith("moved")))
    print("axis1 unchanged:", sum(1 for r in R
                                  if r["axis1_input"].startswith("unchanged")))
    print("guarded:", sum(1 for r in R
                          if r["verify_numbers_guard"].startswith("yes")))
    print("WROTE", CSV_OUT)


if __name__ == "__main__":
    main()
