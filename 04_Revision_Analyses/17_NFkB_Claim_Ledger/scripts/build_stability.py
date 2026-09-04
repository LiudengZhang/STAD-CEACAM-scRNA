"""
stability.csv - one row per ledger claim, five independent axes.

  (a) seed        this module's own 7-seed sweep at fixed data, metric and
                  settings  -> outputs/signal_runs.csv, counted_claims_by_seed.csv
  (b) metric      the four ranking metrics of 16_GSEA_Metric_Sensitivity
                  -> its counted_claims_by_metric_live.csv, nfkb_by_metric_wide_live.csv
  (c) replication the change of unit in 15_Pseudobulk_Sample_Level
                  -> its counted_claims.csv, nfkb_comparison.csv (limma-voom, DESeq2)
  (d) MAST        14_MAST_Specification -> nfkb_by_specification.csv. STILL RUNNING
                  at the time of writing: specification A (~ condition + cngeneson)
                  and T (the Welch comparator) are complete for all 26 contrasts;
                  A0 is partial and B (the mixed model) and dep (the deposited
                  call) are not yet available, so those are left empty rather
                  than guessed.
  (e) sound input 13_R1.8_Neutrophil_Rebuilt_Recompute -> counted_claims.csv
                  (the `sound13` row) and nfkb_per_celltype_sound13.csv

Cells read: agrees / disagrees / partly / not applicable / untestable / pending.
"""
import sys
import numpy as np
import pandas as pd
from gsea_common import ROOT, OUT

A = pd.read_csv(OUT / "counted_claims_by_seed.csv")
A = A[A["degset"] == "live"]
Apost = A[A["phase"] == "post"]
Apre = A[A["phase"] == "pre"]
SIG = pd.read_csv(OUT / "signal.csv", comment="#")
SIGl = SIG[SIG["degset"] == "live"].set_index(["cell_type", "contrast"])

M16 = ROOT / "02_New_Analyses/16_GSEA_Metric_Sensitivity/outputs"
B = pd.read_csv(M16 / "counted_claims_by_metric_live.csv").set_index("metric")
BW = pd.read_csv(M16 / "nfkb_by_metric_wide_live.csv").set_index(["cell", "phase"])

M15 = ROOT / "02_New_Analyses/15_Pseudobulk_Sample_Level/outputs"
C = pd.read_csv(M15 / "counted_claims.csv").set_index("label")
CW = pd.read_csv(M15 / "nfkb_comparison.csv").set_index(["cell_type", "phase"])

M14 = ROOT / "02_New_Analyses/14_MAST_Specification/outputs/nfkb_by_specification.csv"
D = pd.read_csv(M14)
DA = D[D["spec"] == "A"].set_index(["cell", "phase"])

M13 = ROOT / "02_New_Analyses/13_R1.8_Neutrophil_Rebuilt_Recompute/outputs"
E = pd.read_csv(M13 / "counted_claims.csv").set_index("label")
EW = pd.read_csv(M13 / "nfkb_per_celltype_sound13.csv")
EW = EW[EW["method"] == "ttest"].set_index(["cell_type", "phase"])


def rng(s):
    s = list(s)
    return f"{min(s)}" if min(s) == max(s) else f"{min(s)}-{max(s)}"


def d_counts(phase, pred):
    """MAST spec A counted over the 13 contrasts of one phase."""
    g = DA.xs(phase, level="phase")
    return pred(g), len(g)


R = []


def row(cid, quoted, a, av, b, bv, c, cv, dd, dv, e, ev, sig, verdict):
    R.append(dict(claim_id=cid, quoted_value=quoted,
                  a_seed=a, a_seed_value=av,
                  b_ranking_metric=b, b_metric_value=bv,
                  c_replication_unit=c, c_replication_value=cv,
                  d_MAST_specification=dd, d_MAST_value=dv,
                  e_sound_input_recompute=e, e_sound_value=ev,
                  contrast_has_signal=sig, verdict=verdict))


gA = DA.xs("post", level="phase")
gApre = DA.xs("pre", level="phase")
mast_post_pos = int((gA["nfkb_nes"] > 0).sum())
mast_post_q05 = int(((gA["nfkb_nes"] > 0) & (gA["nfkb_fdr_q"] < 0.05)).sum())
mast_post_r1 = int((gA["nfkb_rank"] == 1).sum())
mast_pre_pos = int((gApre["nfkb_nes"] > 0).sum())
mast_pre_q05 = int(((gApre["nfkb_nes"] > 0) & (gApre["nfkb_fdr_q"] < 0.05)).sum())

# ---- M01 / R01a : 12 of 13 positive after treatment
row("M01", "12 of 13",
    "agrees", f"12 at all 7 seeds",
    "agrees", "12 under all four metrics (16_)",
    "agrees (stronger)", "13 of 13 under limma-voom and DESeq2 (15_)",
    "disagrees", f"{mast_post_pos} of 13 under MAST spec A (~condition+cngeneson)",
    "agrees", "12 of 13 on sound13 (13_)",
    "8 of 13 post contrasts have signal", "SAFE")

# ---- M02 : q<0.05 in four
row("M02", "4",
    "agrees", "4 at all 7 seeds",
    "disagrees", "published 4, signed_p 6, logFC 8, t-stat 9 (16_)",
    "disagrees (higher)", "8 under both limma-voom and DESeq2 (15_)",
    "disagrees (lower)", f"{mast_post_q05} under MAST spec A",
    "disagrees (lower)", "3 on sound13 (13_)",
    "n/a - a count over 13 contrasts", "FLOOR, not an estimate")

# ---- M03 / R02 : top-ranked in five
row("M03", "5",
    "DISAGREES", "6,5,4,6,5,5,5 across seeds 1/7/13/42/101/777/2026; "
                 "6 at the published seed 42 in a fresh run",
    "disagrees", "6 under all four metrics (16_)",
    "partly", "4 (limma-voom) / 5 (DESeq2) (15_)",
    "disagrees", f"{mast_post_r1} under MAST spec A",
    "disagrees", "3 on sound13 (13_)",
    "n/a - a count over 13 contrasts", "UNSTABLE - do not print as an integer")

# ---- M04 : Pericyte 31st of 49
pr = Apost["pericyte_rank"]
row("M04", "31 of 49",
    "partly", f"rank {rng(pr)} of 49 across 7 seeds",
    "disagrees", "rank 32 / 22 / 12 / 12 across the four metrics (16_)",
    "disagrees", "23 of 49 (limma-voom), 3 of 49 (DESeq2) (15_)",
    "disagrees", f"rank {int(DA.loc[('Pericyte','post'),'nfkb_rank'])} of "
                 f"{int(DA.loc[('Pericyte','post'),'n_sets'])} (spec A)",
    "partly", "30 of 49 on sound13 (13_)",
    SIGl.loc[("Pericyte", "post"), "has_signal"],
    "ORDINAL - denominator is a min_size artefact; at the min_size=5 the "
    "printed Methods declare it is 33 of 50")

# ---- M05 : B cells 37th of 38
br = Apost["bcells_rank"]
row("M05", "37 of 38",
    "agrees", f"rank {rng(br)} of 38 at all 7 seeds",
    "partly", "rank 37 (published, signed_p) / 38 (logFC, t-stat) of 38 (16_)",
    "DISAGREES", "11 of 44 (limma-voom), 14 of 44 (DESeq2); both positive (15_)",
    "partly", f"rank {int(DA.loc[('B_cells','post'),'nfkb_rank'])} of "
              f"{int(DA.loc[('B_cells','post'),'n_sets'])} (spec A), NES "
              f"{DA.loc[('B_cells','post'),'nfkb_nes']:.3f}",
    "partly", "38 of 38 on sound13 (13_)",
    SIGl.loc[("B_cells", "post"), "has_signal"],
    "ORDINAL - stable to seed and metric, not to replication unit; at min_size=5 "
    "it is 42 of 43")

# ---- M06 : the one population with a negative score
row("M06", "1",
    "agrees", "1 at all 7 seeds",
    "agrees", "B cells is the only negative under all four metrics (16_)",
    "DISAGREES", "0 of 13 negative under both sample-level methods (15_)",
    "disagrees", f"{int((gA['nfkb_nes']<0).sum())} negative under MAST spec A",
    "agrees", "1 on sound13 (13_)",
    SIGl.loc[("B_cells", "post"), "has_signal"],
    "SAFE within the per-cell t-test; not a property of the data")

# ---- M07 : 6 of 13 positive before treatment
row("M07", "6 of 13",
    "agrees", "6 at all 7 seeds",
    "disagrees", "published 6, the other three metrics 7 (one degenerate) (16_)",
    "disagrees", "8 of 13 (limma-voom), 9 of 13 (DESeq2) (15_)",
    "disagrees", f"{mast_pre_pos} of 13 under MAST spec A",
    "disagrees", "5 of 13 on sound13 (13_)",
    "5 of 13 pre contrasts have signal", "UNSTABLE across designs")

# ---- M08 : q<0.05 in one before treatment
row("M08", "1",
    "partly", f"{rng(Apre['n_positive_q05'])} across 7 seeds",
    "disagrees", "published 1, signed_p 4, logFC 4, t-stat 5 (16_)",
    "disagrees", "1 (limma-voom), 4 (DESeq2) (15_)",
    "disagrees", f"{mast_pre_q05} under MAST spec A",
    "disagrees", "3 on sound13 (13_)",
    "5 of 13 pre contrasts have signal", "FLOOR, not an estimate")

# ---- M09 / R05a : MoMac pre NES 1.06 q 0.47
s = SIGl.loc[("MoMac", "pre")]
row("M09", "NES 1.06 / q 0.47",
    "DISAGREES", f"NES {rng([f'{v:.3f}' for v in Apre['momac_nes']])}, "
                 f"q {rng([f'{v:.3f}' for v in Apre['momac_q']])}, rank 3 at all "
                 f"7 seeds",
    "DISAGREES", "1.098/q0.407 (published), 1.190/0.293, 1.000 degenerate/0.628, "
                 "1.749/q0.014 (t-stat: significant) (16_)",
    "DISAGREES", "+1.228 q0.498 (limma-voom), +1.676 q0.031 (DESeq2) (15_)",
    "DISAGREES", f"{DA.loc[('MoMac','pre'),'nfkb_nes']:.3f} q "
                 f"{DA.loc[('MoMac','pre'),'nfkb_fdr_q']:.3f} (spec A) - sign flip",
    "DISAGREES", "-1.001 q0.892 rank 21 of 45 on sound13 (13_) - sign flip",
    f"{s['has_signal']} (max |NES| {s['max_abs_nes_any_set']}, "
    f"{s['n_sets_fdr25']} sets at FDR<0.25, only {s['n_sets_same_es_sign']} of "
    f"{s['n_hallmark_sets_tested']} sets share the positive ES sign)",
    "NO SIGNAL - must not be quoted in any direction")

# ---- M10 / R05b : MoMac post 2.23, strongest
s = SIGl.loc[("MoMac", "post")]
row("M10", "NES 2.23 / q<0.001 / strongest",
    "agrees", f"NES {rng([f'{v:.3f}' for v in Apost['momac_nes']])}, q 0.0000, "
              f"rank 1 and strongest of the 13 at all 7 seeds",
    "agrees", "2.225 / 2.779 / 3.667 / 3.679, rank 1, q 0 under all four (16_)",
    "agrees", "+2.410 rank 1 (limma-voom), +1.943 rank 2 (DESeq2), both q<0.001 (15_)",
    "agrees (weaker)", f"{DA.loc[('MoMac','post'),'nfkb_nes']:.3f} q "
                       f"{DA.loc[('MoMac','post'),'nfkb_fdr_q']:.3f} rank "
                       f"{int(DA.loc[('MoMac','post'),'nfkb_rank'])} (spec A)",
    "agrees", "2.228 q0.000 rank 1 of 43 on sound13 (13_)",
    f"{s['has_signal']} (max |NES| {s['max_abs_nes_any_set']}, "
    f"{s['n_sets_fdr25']} sets at FDR<0.25)",
    "SAFE - the most robust number in the analysis")

# ---- M11 / R06 : Epithelial post NES 1.86
s = SIGl.loc[("Epithelial", "post")]
row("M11", "NES 1.86 (rank 1, q<0.001 in the letter)",
    "agrees", f"NES {rng([f'{v:.3f}' for v in Apost['epithelial_nes']])}, rank 1 "
              f"at all 7 seeds",
    "agrees", f"{BW.loc[('Epithelial','post'),'published_nes']:.3f} / "
              f"{BW.loc[('Epithelial','post'),'signed_p_nes']:.3f} / "
              f"{BW.loc[('Epithelial','post'),'logfc_nes']:.3f} / "
              f"{BW.loc[('Epithelial','post'),'tstat_nes']:.3f}, rank 1 in all (16_)",
    "agrees", f"+{CW.loc[('Epithelial','post'),'limma_nes']:.3f} rank "
              f"{int(CW.loc[('Epithelial','post'),'limma_rank'])}, "
              f"+{CW.loc[('Epithelial','post'),'deseq2_nes']:.3f} rank "
              f"{int(CW.loc[('Epithelial','post'),'deseq2_rank'])} (15_)",
    "agrees", f"{DA.loc[('Epithelial','post'),'nfkb_nes']:.3f} q "
              f"{DA.loc[('Epithelial','post'),'nfkb_fdr_q']:.3f} rank "
              f"{int(DA.loc[('Epithelial','post'),'nfkb_rank'])} (spec A)",
    "agrees", f"{EW.loc[('Epithelial','post'),'nes']:.3f} rank "
              f"{int(EW.loc[('Epithelial','post'),'rank'])} on sound13 (13_)",
    f"{s['has_signal']}", "SAFE")

row("M12", "direction only",
    "untestable", "not measured here - a different module's prerank",
    "untestable", "16_ did not cover the adaptive-lineage GSEA",
    "untestable", "15_ did not cover it",
    "untestable", "14_ did not cover it",
    "untestable", "13_ did not cover it",
    "not measured",
    "UNTESTED - no module has examined the Fig. S10D/E enrichment")

# ---- para 70 / Fig 5N-5I
p70 = OUT / "para70_summary.csv"
p70d = pd.read_csv(p70) if p70.exists() else None
P70KEY = {"M13": ("Epithelial", "Inflammatory Response"),
          "M14": ("Epithelial", "Hypoxia"),
          "M15": ("Fibroblast", "Inflammatory Response"),
          "M16": ("MoMac", "Inflammatory Response"),
          "M17": ("MoMac", "Epithelial Mesenchymal Transition"),
          "M18": ("MoMac", "Hypoxia"),
          "M19": ("MoMac", "Angiogenesis")}
for cid, (cell, pw) in P70KEY.items():
    q = {"M13": "1.48", "M14": "1.56", "M15": "1.30", "M16": "1.67",
         "M17": "1.66", "M18": "1.53", "M19": "1.50"}[cid]
    if p70d is not None:
        r = p70d[(p70d["cell"] == cell) & (p70d["pathway"] == pw)]
        av = (f"NES {r.iloc[0]['nes_seed42']:.4f} at seed 42, "
              f"{r.iloc[0]['nes_min']:.4f}-{r.iloc[0]['nes_max']:.4f} over 7 seeds; "
              f"printed {q}" if len(r) else "not reproduced")
        a = ("agrees" if len(r) and abs(r.iloc[0]["nes_seed42"] - float(q)) < 0.005
             else "DISAGREES")
    else:
        a, av = "pending", "para70_check.py had not finished when this was written"
    row(cid, q, a, av,
        "untestable", "16_ swept only the NF-kB per-cell-type pipeline",
        "untestable", "15_ swept only the NF-kB per-cell-type pipeline",
        "untestable", "14_ swept only the NF-kB per-cell-type pipeline",
        "not applicable", "a different pipeline: DE computed inside the panel "
                          "script, no 10% detection filter, min_size=5",
        "not measured by signal.csv (different ranked list)",
        "SEPARATE PIPELINE - see FINDINGS")

row("M20", "all four positive",
    "pending" if p70d is None else "see para70_all_terms.csv", "",
    "untestable", "", "untestable", "", "untestable", "",
    "not applicable", "", "not measured", "SEPARATE PIPELINE")

row("M21", "Methods declare min_size = 5",
    "not applicable", "settings, not a seed question",
    "agrees with 16_", "16_ measured max |dNES| = 0.000 across min_size 5/10/25",
    "not applicable", "", "not applicable", "",
    "not applicable", "",
    "n/a",
    "WRONG AS PRINTED - the numbers were produced at min_size = 15 "
    "(recompute_deg.py:351). NES is unaffected; every printed ordinal and "
    "denominator is not. See ordinal_denominator.csv")

row("M22", "Fig 5A, nine bars",
    "untestable", "the MAST prerank table's seed sweep is not in scope here",
    "untestable", "", "untestable", "", "untestable", "",
    "not applicable", "Round_5 MAST prerank, a third pipeline again",
    "n/a", "REPRODUCES the printed panel to 0.0001 NES (PROVENANCE, 2026-09-01)")

row("M23", "Fig 5H, 26 NES values",
    "agrees", "the 26 values move by at most 0.12 NES over 7 seeds",
    "agrees", "sign agrees on 25 of 26 under all four metrics (16_)",
    "partly", "sign agrees on 12-13 of 13 post, 5-8 of 13 pre (15_)",
    "partly", "spec A flips sign on 4 post contrasts", "partly",
    "sound13 flips sign on 3 contrasts incl. MoMac pre and B cells post",
    "8 of 13 post and 5 of 13 pre contrasts have signal",
    "The post half is sound; the pre half plots values from contrasts with no "
    "signal")

# response-letter rows mirror the main text
for cid, src, note in [("R01", "M01+M02", "same values"),
                       ("R02", "M03", "same value"),
                       ("R03", "M04+M05", "same values"),
                       ("R04", "M07+M01", "same values"),
                       ("R05", "M09+M10", "same values"),
                       ("R06", "M11", "same value, plus the rank-1 claim"),
                       ("R07", "M12", "same qualitative claim")]:
    row(cid, f"as {src}", "see " + src, note, "see " + src, "", "see " + src, "",
        "see " + src, "", "see " + src, "", "see " + src,
        "INHERITS the verdict of " + src)

for cid, note in [("F01", "= M22"), ("F02", "= M23"), ("F03", "= M01-M06"),
                  ("F04", "= M07-M10"), ("F05", "= M13-M15"),
                  ("F06", "= M13-M19")]:
    row(cid, note, "see " + note[2:], "", "see " + note[2:], "",
        "see " + note[2:], "", "see " + note[2:], "", "see " + note[2:], "",
        "see " + note[2:], "INHERITS " + note[2:])

d = pd.DataFrame(R)
d.to_csv(OUT / "stability.csv", index=False)
pd.set_option("display.width", 250)
print(d[["claim_id", "quoted_value", "a_seed", "b_ranking_metric",
         "c_replication_unit", "d_MAST_specification",
         "e_sound_input_recompute", "verdict"]].to_string(index=False))
sys.exit(0)
