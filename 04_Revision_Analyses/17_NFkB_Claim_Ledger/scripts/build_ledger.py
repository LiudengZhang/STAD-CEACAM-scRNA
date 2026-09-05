"""
ledger.csv - one row per printed GSEA-derived claim in the shipped paper.

The claims are taken from the authoritative sources, not from reading prose:

  main text     04_Manuscript_R1/01_Main_Text/edits.py, and the paragraphs it
                produces in Manuscript_R1_clean.docx (read only)
  what is
  checked       04_Manuscript_R1/verify_numbers.py - the expression each claim
                is checked by names the table and the column it reads
  letter        04_Manuscript_R1/05_Response_to_Reviewers/
                apply_consistency_fixes.py, and the shipped v3_clean.docx
  panels        03_Final_Panels/PROVENANCE.csv rows attributed to NF-kB
                analyses: Figure 5A, 5H, 5I, 5N, S9E, S10C

`value_in_source` is read back out of the named table by this script, so the
"quoted vs source" column is a measurement and not a transcription.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from gsea_common import ROOT, OUT

NF = ROOT / "04_Revision_Analyses/07_R1.8_NFkB_Specificity/outputs/nfkb_per_celltype.csv"
PP = ROOT / "04_Revision_Analyses/08_R2.1_PreTx_Inflammatory/outputs/nfkb_pre_vs_post.csv"
MOMAC_A = (ROOT.parent / "Round_5" / "02_Preparation_for_Panels" / "GSEA"
           / "MoMac_mast_prerank_gsea.csv")
G5 = (ROOT / "03_Final_Panels/05_Figure_5/05_GSEA_Summary/"
      "gsea_data/gsea_combined_5types.csv")
GM = (ROOT / "03_Final_Panels/05_Figure_5/05_GSEA_Summary/"
      "gsea_data/gsea_momac.csv")
AD = ROOT / "04_Revision_Analyses/09_R2.2_Adaptive_Immune/outputs/adaptive_hallmark_top.csv"

nf = pd.read_csv(NF)
nf = nf[nf["method"] == "ttest"]
post = nf[nf["phase"] == "post"].set_index("cell_type")
pre = nf[nf["phase"] == "pre"].set_index("cell_type")
g5 = pd.read_csv(G5)
gm = pd.read_csv(GM)


def rel(p):
    return str(Path(p).resolve()).replace(str(ROOT) + "/", "")


def g5v(ct, pw, col="NES"):
    r = g5[(g5["CellType"] == ct) & (g5["Pathway"] == pw)]
    return None if r.empty else float(r.iloc[0][col])


def gmv(pw, col="NES"):
    r = gm[gm["Pathway"] == pw]
    return None if r.empty else float(r.iloc[0][col])


C = []


def add(**kw):
    kw.setdefault("cell_type", "")
    kw.setdefault("contrast", "")
    kw.setdefault("DE_set_it_is_counted_over", "")
    kw.setdefault("source_table_path", "")
    kw.setdefault("source_column", "")
    kw.setdefault("quoted_value", "")
    kw.setdefault("value_in_source", "")
    kw.setdefault("checked_by_verify_numbers", "")
    kw.setdefault("note", "")
    C.append(kw)


MT = "main text, Results, NF-kB paragraph (clean.docx para 68)"
MT70 = "main text, Results, downstream cascade (clean.docx para 70)"
LET = "response letter v3_clean"
LIVE = ("live (12_R1.8_DEG_Recompute/outputs/gsea, read by "
        "07_.../nfkb_specificity.py load_gsea)")

# ---------------------------------------------------------------- main text
add(claim_id="M01", where_printed=MT, claim_type="count",
    verbatim_text="positive enrichment of the Hallmark TNFα signaling via NF-κB "
                  "gene set in non-responders in 12 of 13 populations after treatment",
    cell_type="all 13", contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes > 0, method=ttest, phase=post",
    quoted_value="12 of 13", value_in_source=f"{int((post['nes']>0).sum())} of {len(post)}",
    checked_by_verify_numbers="cell types with positive post NES == 12")

add(claim_id="M02", where_printed=MT, claim_type="count",
    verbatim_text="reaching FDR q < 0.05 in four (monocytes/macrophages, dendritic "
                  "cells, epithelial cells and fibroblasts; Fig. S9E)",
    cell_type="MoMac, DC, Epithelial, Fibroblast", contrast="post",
    DE_set_it_is_counted_over=LIVE, source_table_path=rel(NF),
    source_column="nes > 0 & fdr_q < 0.05",
    quoted_value="4",
    value_in_source=str(int(((post['nes']>0)&(post['fdr_q']<0.05)).sum())),
    checked_by_verify_numbers="post NES with FDR < 0.05 == 4")

top = sorted(post.index[post["rank"] == 1])
add(claim_id="M03", where_printed=MT, claim_type="ordinal count",
    verbatim_text="it was the top-ranked Hallmark gene set in those four populations "
                  "and in mast cells",
    cell_type="; ".join(top), contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="rank == 1",
    quoted_value="5", value_in_source=f"{len(top)} ({', '.join(top)})",
    checked_by_verify_numbers="populations where TNFa/NF-kB is the top-ranked "
                              "Hallmark set == 5, and the named set")

add(claim_id="M04", where_printed=MT, claim_type="ordinal",
    verbatim_text="but ranked 31st of 49 sets in pericytes",
    cell_type="Pericyte", contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="rank / n_sets",
    quoted_value="31 of 49",
    value_in_source=f"{int(post.loc['Pericyte','rank'])} of {int(post.loc['Pericyte','n_sets'])}",
    checked_by_verify_numbers="Pericyte rank == 31; Hallmark sets tested == 49")

add(claim_id="M05", where_printed=MT, claim_type="ordinal",
    verbatim_text="and 37th of 38 in B cells",
    cell_type="B_cells", contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="rank / n_sets",
    quoted_value="37 of 38",
    value_in_source=f"{int(post.loc['B_cells','rank'])} of {int(post.loc['B_cells','n_sets'])}",
    checked_by_verify_numbers="B cells rank == 37; Hallmark sets tested == 38")

add(claim_id="M06", where_printed=MT, claim_type="count",
    verbatim_text="the one population with a negative score",
    cell_type="B_cells", contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes < 0",
    quoted_value="1", value_in_source=str(int((post["nes"] < 0).sum())),
    checked_by_verify_numbers="B cells are claimed to be the only population with a "
                              "negative post-treatment NES")

add(claim_id="M07", where_printed=MT, claim_type="count",
    verbatim_text="before treatment the gene set is positively enriched in only 6 of "
                  "13 populations",
    cell_type="all 13", contrast="pre", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes > 0, phase=pre",
    quoted_value="6 of 13", value_in_source=f"{int((pre['nes']>0).sum())} of {len(pre)}",
    checked_by_verify_numbers="cell types with positive pre NES == 6")

add(claim_id="M08", where_printed=MT, claim_type="count",
    verbatim_text="and reaches q < 0.05 in one",
    cell_type="Endothelial_cells", contrast="pre", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes > 0 & fdr_q < 0.05, phase=pre",
    quoted_value="1",
    value_in_source=str(int(((pre['nes']>0)&(pre['fdr_q']<0.05)).sum())),
    checked_by_verify_numbers="pre NES with FDR < 0.05 == 1")

add(claim_id="M09", where_printed=MT, claim_type="point value",
    verbatim_text="monocytes/macrophages are not significant (NES = 1.06, q = 0.47)",
    cell_type="MoMac", contrast="pre", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes, fdr_q",
    quoted_value="NES 1.06 / q 0.47",
    value_in_source=f"NES {pre.loc['MoMac','nes']:.4f} / q {pre.loc['MoMac','fdr_q']:.4f}",
    checked_by_verify_numbers="MoMac pre NES == 1.06 (tol 0.01); MoMac pre FDR q == 0.47")

add(claim_id="M10", where_printed=MT, claim_type="point value + superlative",
    verbatim_text="whereas after treatment they show the strongest enrichment of any "
                  "population (NES = 2.23, q < 0.001; Fig. S10C)",
    cell_type="MoMac", contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes, fdr_q; max over cell types",
    quoted_value="NES 2.23 / q < 0.001, rank 1 of 13 populations",
    value_in_source=f"NES {post.loc['MoMac','nes']:.4f} / q {post.loc['MoMac','fdr_q']:.4f}; "
                    f"largest post NES is {post['nes'].idxmax()}",
    checked_by_verify_numbers="MoMac post NES == 2.23; MoMac post FDR q == 0.0; "
                              "MoMac has the strongest post NES")

add(claim_id="M11", where_printed=MT, claim_type="point value",
    verbatim_text="while showing NF-κB target-gene enrichment itself (NES = 1.86)",
    cell_type="Epithelial", contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes",
    quoted_value="1.86", value_in_source=f"{post.loc['Epithelial','nes']:.4f}",
    checked_by_verify_numbers="Epithelial post NES == 1.86")

add(claim_id="M12", where_printed=MT, claim_type="qualitative (direction)",
    verbatim_text="while pathway-level analysis separated the groups, with "
                  "interferon-γ programs enriched in responders and proliferation "
                  "programs in non-responders (Fig. S10D, S10E)",
    cell_type="CD4+ T, CD8+ T, NK, B, plasma", contrast="post",
    DE_set_it_is_counted_over="09_R2.2_Adaptive_Immune (its own prerank)",
    source_table_path=rel(AD), source_column="NES by lineage and set",
    quoted_value="direction only, no number printed",
    value_in_source="see adaptive_hallmark_top.csv",
    checked_by_verify_numbers="not checked as a GSEA value; only the BH claim on "
                              "adaptive_state_tests.csv is checked")

# ------------------------------------------------ main text, para 70, Fig 5N/5I
P70 = ("separate GSEA pipeline: differential expression computed inside the panel "
       "script (05_G/create_panel_g_gsea_2types.py, 05_GSEA_Summary/"
       "run_gsea_5celltypes.py, run_gsea_momac.py), prerank at min_size=5")
for cid, ct, pw, q, qp, src in [
    ("M13", "Epithelial", "Inflammatory Response", 1.48, 0.02, "g5"),
    ("M14", "Epithelial", "Hypoxia", 1.56, 0.006, "g5"),
    ("M15", "Fibroblast", "Inflammatory Response", 1.30, 0.02, "g5"),
    ("M16", "MoMac", "Inflammatory Response", 1.67, None, "gm"),
    ("M17", "MoMac", "Epithelial Mesenchymal Transition", 1.66, None, "gm"),
    ("M18", "MoMac", "Hypoxia", 1.53, None, "gm"),
    ("M19", "MoMac", "Angiogenesis", 1.50, None, "gm"),
]:
    v = g5v(ct, pw) if src == "g5" else gmv(pw)
    p = g5v(ct, pw, "NOM_pval") if src == "g5" else gmv(pw, "NOM_pval")
    add(claim_id=cid, where_printed=MT70, claim_type="point value",
        verbatim_text=f"{ct} {'showed' if src=='g5' else ''} enrichment for "
                      f"{pw} (NES = {q:.2f}"
                      + (f", P = {qp}" if qp is not None else "") + ")",
        cell_type=ct, contrast="post",
        DE_set_it_is_counted_over="neither `live` nor `sound13` - " + P70,
        source_table_path=rel(G5 if src == "g5" else GM),
        source_column="NES / NOM_pval",
        quoted_value=f"NES {q:.2f}" + (f", P {qp}" if qp is not None else ""),
        value_in_source=f"NES {v:.4f}, NOM P {p:.4f}" if v is not None else "not found",
        checked_by_verify_numbers="NOT CHECKED - verify_numbers.py has no check for "
                                  "any para-70 NES",
        note="")

add(claim_id="M20", where_printed=MT70, claim_type="qualitative (direction)",
    verbatim_text="Systematic GSEA across monocytes/macrophages, epithelial cells, "
                  "fibroblasts, and dendritic cells revealed that NF-κB-driven "
                  "pathways were enriched in all four cell types in post-treatment "
                  "non-responders (Fig. 5N)",
    cell_type="MoMac, Epithelial, Fibroblast, DC", contrast="post",
    DE_set_it_is_counted_over="neither - " + P70,
    source_table_path=rel(G5) + "; " + rel(GM), source_column="NES",
    quoted_value="all four positive", value_in_source="see stability.csv",
    checked_by_verify_numbers="NOT CHECKED")

add(claim_id="M21", where_printed="main text, Methods (clean.docx para 120)",
    claim_type="method statement",
    verbatim_text="GSEA was then performed using GSEApy prerank against MSigDB "
                  "Hallmark 2020 gene sets (min_size = 5, max_size = 500; 1,000 "
                  "permutations; seed = 42)",
    cell_type="all 13", contrast="pre and post", DE_set_it_is_counted_over=LIVE,
    source_table_path="04_Revision_Analyses/12_R1.8_DEG_Recompute/scripts/recompute_deg.py:351",
    source_column="gp.prerank(..., min_size=15, max_size=500, seed=42)",
    quoted_value="min_size = 5", value_in_source="min_size = 15",
    checked_by_verify_numbers="NOT CHECKED - verify_numbers.py checks values, not "
                              "the Methods' settings",
    note="The declared min_size is not the one the printed numbers were produced at. "
         "See ordinal_denominator.csv: NES is unchanged, but every printed ordinal "
         "and denominator is a min_size=15 quantity.")

add(claim_id="M22", where_printed="main text, Figure 5 legend",
    claim_type="panel description",
    verbatim_text="(A) Top 9 Hallmark pathways enriched in Monocytes/Macrophages",
    cell_type="MoMac", contrast="post",
    DE_set_it_is_counted_over="neither - MAST prerank, Round_5 GSEA/MoMac_mast_prerank_gsea.csv",
    source_table_path=rel(MOMAC_A), source_column="NES",
    quoted_value="9 bars", value_in_source="35 sets in the table, top 9 by |NES| plotted",
    checked_by_verify_numbers="NOT CHECKED")

add(claim_id="M23", where_printed="main text, Figure 5 legend",
    claim_type="panel description",
    verbatim_text="(H) Radar plot comparing NF-κB pathway NES between pre- and "
                  "post-treatment across cell types.",
    cell_type="all 13", contrast="pre and post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes",
    quoted_value="26 NES values", value_in_source="26 rows, method=ttest",
    checked_by_verify_numbers="indirectly - the panel reads the checked table")

# ------------------------------------------------------------- response letter
add(claim_id="R01", where_printed=LET + ", para 84", claim_type="count",
    verbatim_text="enriched toward post-treatment non-response in 12 of 13 "
                  "populations, with four reaching FDR q < 0.05",
    cell_type="all 13", contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes > 0; nes > 0 & fdr_q < 0.05",
    quoted_value="12 of 13; four",
    value_in_source=f"{int((post['nes']>0).sum())} of {len(post)}; "
                    f"{int(((post['nes']>0)&(post['fdr_q']<0.05)).sum())}",
    checked_by_verify_numbers="same checks as M01/M02")

add(claim_id="R02", where_printed=LET + ", para 84", claim_type="ordinal count",
    verbatim_text="and it is the top-ranked Hallmark set in five populations",
    cell_type="; ".join(top), contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="rank == 1",
    quoted_value="5", value_in_source=str(len(top)),
    checked_by_verify_numbers="same check as M03")

add(claim_id="R03", where_printed=LET + ", para 84", claim_type="ordinal",
    verbatim_text="while ranking 31st of 49 in pericytes and 37th of 38 in B cells, "
                  "so the enrichment is graded rather than uniform",
    cell_type="Pericyte; B_cells", contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="rank / n_sets",
    quoted_value="31 of 49; 37 of 38",
    value_in_source=f"{int(post.loc['Pericyte','rank'])} of {int(post.loc['Pericyte','n_sets'])}; "
                    f"{int(post.loc['B_cells','rank'])} of {int(post.loc['B_cells','n_sets'])}",
    checked_by_verify_numbers="same checks as M04/M05")

add(claim_id="R04", where_printed=LET + ", para 98", claim_type="count",
    verbatim_text="the Hallmark TNFα/NF-κB set is positively enriched in 6 of 13 "
                  "populations before treatment and in 12 of 13 after",
    cell_type="all 13", contrast="pre and post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(PP), source_column="pre_nes > 0; post_nes > 0",
    quoted_value="6; 12",
    value_in_source=f"{int((pd.read_csv(PP)['pre_nes']>0).sum())}; "
                    f"{int((pd.read_csv(PP)['post_nes']>0).sum())}",
    checked_by_verify_numbers="S10C shows 6 of 13 positive before and 12 of 13 after")

add(claim_id="R05", where_printed=LET + ", para 98", claim_type="point value",
    verbatim_text="monocytes/macrophages move from not significant before treatment "
                  "(NES = 1.06, q = 0.47) to the strongest enrichment of any "
                  "population after it (NES = 2.23, q < 0.001)",
    cell_type="MoMac", contrast="pre and post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(PP), source_column="pre_nes, pre_fdr, post_nes, post_fdr",
    quoted_value="1.06 / 0.47 ; 2.23 / <0.001",
    value_in_source=f"{pre.loc['MoMac','nes']:.4f} / {pre.loc['MoMac','fdr_q']:.4f} ; "
                    f"{post.loc['MoMac','nes']:.4f} / {post.loc['MoMac','fdr_q']:.4f}",
    checked_by_verify_numbers="MoMac pre NES (S10C) == 1.06; MoMac post NES (S10C) == 2.23")

add(claim_id="R06", where_printed=LET + ", para 86", claim_type="point value + ordinal",
    verbatim_text="while epithelial cells themselves carry a clear NF-κB target-gene "
                  "signature (NES = 1.86, FDR q < 0.001, the top-ranked Hallmark set "
                  "in that compartment; Fig. S9E)",
    cell_type="Epithelial", contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes, fdr_q, rank",
    quoted_value="1.86 / q<0.001 / rank 1",
    value_in_source=f"{post.loc['Epithelial','nes']:.4f} / "
                    f"{post.loc['Epithelial','fdr_q']:.4f} / "
                    f"rank {int(post.loc['Epithelial','rank'])} of "
                    f"{int(post.loc['Epithelial','n_sets'])}",
    checked_by_verify_numbers="Epithelial post NES == 1.86 (the rank is checked only "
                              "through the top-ranked set membership check)")

add(claim_id="R07", where_printed=LET + ", para 107", claim_type="qualitative (direction)",
    verbatim_text="post-treatment responder CD4+ and CD8+ T cells show stronger "
                  "interferon-γ/inflammatory programs, whereas non-responder T cells "
                  "are enriched for proliferation-associated programs",
    cell_type="CD4+ T, CD8+ T", contrast="post",
    DE_set_it_is_counted_over="09_R2.2_Adaptive_Immune (its own prerank)",
    source_table_path=rel(AD), source_column="NES by lineage and set",
    quoted_value="direction only", value_in_source="see adaptive_hallmark_top.csv",
    checked_by_verify_numbers="NOT CHECKED as a GSEA value")

# ------------------------------------------------------------------- panels
add(claim_id="F01", where_printed="figure panel, Figure 5A (printed letter A; "
                                  "PROVENANCE build dir 05_A)",
    claim_type="plotted values",
    verbatim_text="nine horizontal bars, NES; TNF-alpha Signaling via NF-kB +2.140 "
                  "the largest positive",
    cell_type="MoMac", contrast="post",
    DE_set_it_is_counted_over="neither - MAST prerank on the Round_5 preparation",
    source_table_path=rel(MOMAC_A), source_column="NES",
    quoted_value="+2.140, +1.860, +1.750, +1.748; -1.488, -1.803, -1.875, -2.223, -2.450",
    value_in_source="reproduces the printed panel to 0.0001 NES "
                    "(create_panel_a_momac_enrichment.py docstring, 2026-09-01)",
    checked_by_verify_numbers="NOT CHECKED by verify_numbers.py; checked by "
                              "PROVENANCE.csv reproduces_published=yes",
    note="This panel's NES scale is NOT the same analysis as Fig 5H / S9E: TNFa is "
         "+2.140 here (MAST prerank, Round_5) and +2.233 there (t-test, module 12).")

add(claim_id="F02", where_printed="figure panel, Figure 5H (printed letter H; "
                                  "PROVENANCE build dir 05_F)",
    claim_type="plotted values",
    verbatim_text="radar, 13 cell types x pre/post NES",
    cell_type="all 13", contrast="pre and post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes (method=ttest)",
    quoted_value="26 NES values", value_in_source="26 rows",
    checked_by_verify_numbers="indirect - the table is the checked one")

add(claim_id="F03", where_printed="figure panel, Figure S9E",
    claim_type="plotted values",
    verbatim_text="Hallmark TNFα signaling via NF-κB across the 13 cell types after "
                  "treatment, with the FDR and the rank of the set among all Hallmark "
                  "sets tested in that cell type",
    cell_type="all 13", contrast="post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(NF), source_column="nes, fdr_q, rank, n_sets",
    quoted_value="13 NES + 13 ranks + 13 denominators",
    value_in_source="13 rows",
    checked_by_verify_numbers="M01-M06 are checks on exactly these values")

add(claim_id="F04", where_printed="figure panel, Figure S10C",
    claim_type="plotted values",
    verbatim_text="pre- versus post-treatment NF-κB NES per cell type",
    cell_type="all 13", contrast="pre and post", DE_set_it_is_counted_over=LIVE,
    source_table_path=rel(PP), source_column="pre_nes, post_nes, pre_fdr, post_fdr",
    quoted_value="26 NES values", value_in_source="13 rows x 2 phases",
    checked_by_verify_numbers="MoMac pre/post NES (S10C); 6 of 13 / 12 of 13")

add(claim_id="F05", where_printed="figure panel, Figure 5I (PROVENANCE build dir 05_G)",
    claim_type="plotted values",
    verbatim_text="GSEA running enrichment for epithelial cells and fibroblasts",
    cell_type="Epithelial, Fibroblast", contrast="post",
    DE_set_it_is_counted_over="neither - " + P70,
    source_table_path="03_Final_Panels/05_Figure_5/05_G/"
                      "gsea_{Epithelial,Fibroblast}/gseapy.gene_set.prerank.report.csv",
    source_column="ES / NES", quoted_value="curves, no printed number",
    value_in_source="prerank reports on disk",
    checked_by_verify_numbers="NOT CHECKED")

add(claim_id="F06", where_printed="figure panel, Figure 5N "
                                  "(PROVENANCE build dir 05_GSEA_Summary)",
    claim_type="plotted values",
    verbatim_text="four Hallmark dotplots: monocytes/macrophages, epithelial cells, "
                  "fibroblasts, dendritic cells",
    cell_type="MoMac, Epithelial, Fibroblast, DC", contrast="post",
    DE_set_it_is_counted_over="neither - " + P70,
    source_table_path=rel(GM) + "; " + rel(G5), source_column="NES, NOM_pval, FDR_qval",
    quoted_value="NES per set per cell type",
    value_in_source="50 and 4x50 rows",
    checked_by_verify_numbers="NOT CHECKED")

d = pd.DataFrame(C)
cols = ["claim_id", "where_printed", "claim_type", "verbatim_text", "cell_type",
        "contrast", "DE_set_it_is_counted_over", "source_table_path",
        "source_column", "quoted_value", "value_in_source",
        "checked_by_verify_numbers", "note"]
d = d[cols]
d.to_csv(OUT / "ledger.csv", index=False)
print(d[["claim_id", "claim_type", "quoted_value", "value_in_source"]].to_string(index=False))
print(f"\n{len(d)} claims written to {OUT/'ledger.csv'}")
sys.exit(0)
