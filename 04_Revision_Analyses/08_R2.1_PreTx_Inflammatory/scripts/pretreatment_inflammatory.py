"""
WP9 / Reviewer 2 point R2.1.

  "...better clarify whether the inflammatory programs observed in the tumor
   microenvironment of post-treatment non-responders are already detectable in
   pre-treatment responders or non-responders... whether IL-1b-driven
   inflammation represents a purely acquired resistance mechanism or whether
   elements of this program are already present before treatment."

Three measurements, all comparing the pre-treatment contrast against the
post-treatment one on the same scale:

  1. Abundance of the IL-1b+ inflammatory state as a fraction of the
     monocyte/macrophage compartment, in all four groups.
  2. The IL-1b+ state's own transcriptional signature scored in every
     monocyte/macrophage cell, by group.
  3. TNFa/NF-kB Hallmark enrichment pre versus post, per cell type, from the
     module 12 recompute on .raw - the same producer as Fig. S9E, Welch t-test
     branch. Its NES is already non-responder-relative.

Inputs : submission-tree/01_Raw_Inputs/01_H5AD/MoMac.h5ad
         12_R1.8_DEG_Recompute/outputs/gsea/*_{pre,post}_ttest_hallmark.csv
         13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/nfkb_per_celltype_sound13.csv
             - the REPORT's section 3 only. See ADOPTED_GSEA below.
Outputs: pretx_state_abundance.csv, pretx_signature_scores.csv,
         il1b_signature_genes.csv, nfkb_pre_vs_post.csv,
         pretx_inflammatory_report.txt, README.txt

This module computes these results and does not draw them. The three panels it
carried belonged to a supplementary figure the article does not print, so
drawing them wrote figures that no assembled figure reads and no page carries;
the five tables above are the whole output, and they are what the Results, the
response letter and verify_numbers.py read.

The tables this module writes are the live run, nfkb_pre_vs_post.csv included.
Section 3 of pretx_inflammatory_report.txt is rendered instead from the adopted
sound-input table, through load_adopted_shift(); sections 1 and 2 are this
module's own live run.

That distinction now also travels with the tables. write_outputs_readme() puts
outputs/README.txt beside them saying which of the two NF-kB runs each file is,
because outputs/ is the directory that ships and a reader who opens
nfkb_pre_vs_post.csv there and compares it with Supplementary Table S10 needs
to be told, in that directory, that the two are different runs and neither is
a corrected copy of the other.

--rewrite-report-from-adopted rebuilds the report alone, reading sections 1 and 2
back from this module's own deposited CSVs, so the report can be regenerated
without re-running the analysis and putting new contents behind cited
filenames.
"""

from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD  # noqa: E402
from shared.sample_ids import sample_id_map, to_study_ids  # noqa: E402

RECOMPUTE_GSEA = (Path(__file__).resolve().parents[2]
                  / "12_R1.8_DEG_Recompute" / "outputs" / "gsea")

# ---------------------------------------------------------------- the adopted table
# Section 3 of the report is rendered from the SOUND-INPUT recompute, not from
# this module's own nfkb_pre_vs_post.csv, so that it agrees with the Results and
# the response letter, which take every NF-kB number from that same table.
# RECOMPUTE_GSEA above is the "live" run, computed on the doubly-normalised .X
# of 00_Data_Audit/FINDINGS.md sections 1 and 7.
#
# The analysis below is untouched: nfkb_pre_vs_post.csv is still the live run,
# byte for byte, and one file name still means one set of contents. Six
# quantities differ between the two tables, and only these six:
#     post FDR q < 0.05 count   4  -> 3     rank-1 count      5  -> 3
#     pericyte ordinal         31  -> 30    B-cell ordinal   37  -> 38
#     pre positive count        6  -> 5     pre q < 0.05      1  -> 3
ADOPTED_GSEA = (Path(__file__).resolve().parents[2]
                / "13_R1.8_Neutrophil_Rebuilt_Recompute" / "outputs"
                / "nfkb_per_celltype_sound13.csv")

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

MIN_CELLS = 20
TARGET = "C3_Mac_Inflam_IL1B"
GROUPS = ["Pre-R", "Pre-NR", "Post-R", "Post-NR"]


def four_group_label(row):
    if row["Treatment phase"] == "Pre":
        m = {"Responsed": "Pre-R", "No-response": "Pre-NR"}
        return m.get(row["stomach_pre_grouping"])
    m = {"Responsed": "Post-R", "No-response": "Post-NR"}
    return m.get(row["stomach_post_grouping"])


def mw(a, b):
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(p), float(2.0 * u / (len(a) * len(b)) - 1.0)


def main():
    ad = sc.read_h5ad(MOMAC_H5AD)
    # Resolved before subsetting, so the crosswalk is checked against every
    # specimen in the object rather than only the stomach ones.
    ids = sample_id_map(ad.obs)
    ad = ad[ad.obs["Sample site"] == "Stomach"].copy()
    ad.obs["group4"] = ad.obs.apply(four_group_label, axis=1)
    ad = ad[ad.obs["group4"].isin(GROUPS)].copy()
    ad.obs["minor_cell_state"] = ad.obs["minor_cell_state"].astype(str)

    # ------------------------------------------------- 1. state abundance
    df = ad.obs[["sample", "group4", "minor_cell_state"]].copy()
    df["sample"] = df["sample"].astype(str)
    keep = df.groupby("sample").size()
    df = df[df["sample"].isin(keep[keep >= MIN_CELLS].index)]
    frac = (df.groupby(["sample", "group4"], observed=True)["minor_cell_state"]
            .apply(lambda s: (s == TARGET).mean()).rename("fraction").reset_index())
    # Specimens are named the way Supplementary Table 1 names them. The relabel
    # is in place, after the grouping, so no row moves and no value changes.
    frac["sample"] = to_study_ids(frac["sample"], ids)
    frac.to_csv(OUT / "pretx_state_abundance.csv", index=False)

    # -------------------------------------------- 2. signature score
    # The signature is the state's own markers, derived from this object so the
    # score is internally consistent rather than imported from elsewhere.
    sc.tl.rank_genes_groups(ad, "minor_cell_state", groups=[TARGET],
                            method="wilcoxon", n_genes=50)
    sig = list(sc.get.rank_genes_groups_df(ad, group=TARGET)["names"][:50])
    sc.tl.score_genes(ad, sig, score_name="il1b_signature", random_state=0)

    sdf = ad.obs[["sample", "group4", "il1b_signature"]].copy()
    sdf["sample"] = sdf["sample"].astype(str)
    sdf = sdf[sdf["sample"].isin(keep[keep >= MIN_CELLS].index)]
    score = (sdf.groupby(["sample", "group4"], observed=True)["il1b_signature"]
             .mean().rename("score").reset_index())
    score["sample"] = to_study_ids(score["sample"], ids)
    score.to_csv(OUT / "pretx_signature_scores.csv", index=False)
    pd.Series(sig, name="gene").to_csv(OUT / "il1b_signature_genes.csv", index=False)

    # ----------------------------------------------- 3. NF-kB pre vs post
    # From 12_R1.8_DEG_Recompute, the same producer as Fig. S9E, so every
    # NF-kB number added in this revision comes from one analysis. The prepared
    # GSEA/{pre,post} tables are not used: they descend from MAST runs on a
    # matrix that was never log1p CP10K (00_Data_Audit/FINDINGS.md section 7).
    # The Welch t-test branch is used because that is the test Fig. 5H reports;
    # module 12 ranks non-responders against responders, so no sign flip.
    rows = []
    for phase in ("pre", "post"):
        for f in sorted(RECOMPUTE_GSEA.glob(f"*_{phase}_ttest_hallmark.csv")):
            d = pd.read_csv(f)
            hit = d[d["Term"].str.contains("NF-kB", case=False, regex=False)]
            if not len(hit):
                continue
            r = hit.iloc[0]
            rows.append(dict(phase=phase,
                             cell_type=f.name.replace(f"_{phase}_ttest_hallmark.csv", ""),
                             nes_NRvsR=float(r["NES"]),
                             fdr_q=float(r["FDR q-val"])))
    if not rows:
        raise SystemExit(
            f"no NF-kB rows under {RECOMPUTE_GSEA} - run "
            "12_R1.8_DEG_Recompute/scripts/recompute_deg.py first")
    long = pd.DataFrame(rows)
    nfkb = long.pivot(index="cell_type", columns="phase", values="nes_NRvsR")
    nfkb.columns = [f"{c}_nes" for c in nfkb.columns]
    fdr = long.pivot(index="cell_type", columns="phase", values="fdr_q")
    fdr.columns = [f"{c}_fdr" for c in fdr.columns]
    nfkb = nfkb.join(fdr).reset_index()
    nfkb.to_csv(OUT / "nfkb_pre_vs_post.csv", index=False)

    # -------------------------------------------------------------- report
    # Section 3 is rendered from the ADOPTED sound-input table.
    # nfkb_pre_vs_post.csv, written just above, remains the live run and is not
    # what the report quotes; the report's own SOURCES block says so.
    report = render_report(frac, score, load_adopted_shift(),
                           sources=ADOPTED_SOURCES)
    (OUT / "pretx_inflammatory_report.txt").write_text(report, encoding="utf-8")
    write_outputs_readme()
    print(report)


# ------------------------------------------------------- the outputs/ note
# nfkb_pre_vs_post.csv disagrees with Table S10, on purpose, and a reader who
# opens outputs/ has no way of knowing that. Everything above explains it - the
# module docstring, the ADOPTED_GSEA comment, the report's own SOURCES block -
# and none of it is in the directory the CSV is in. That directory is what
# ships: 06_Code/code/04_Revision_Analyses/08_R2.1_PreTx_Inflammatory/outputs/
# is where the code reads these tables from, where the archive README.txt sends
# a reader looking for the numbers behind a claim, and where the response
# letter's other two analysis-output citations land them. A note anywhere else
# is a note they do not reach.
#
# It is written from this module's own constants rather than typed into the
# directory as a file of its own, so that the paths it names cannot drift from
# the paths the code actually reads. main() writes it with the tables.
README_NAME = "README.txt"

OUTPUTS_README = """\
Outputs of 08_R2.1_PreTx_Inflammatory (reviewer point R2.1)

  pretx_state_abundance.csv    IL-1b+ inflammatory MoMac abundance, by group
  pretx_signature_scores.csv   that state's signature, scored per cell, by group
  il1b_signature_genes.csv     the genes that signature is built from
  nfkb_pre_vs_post.csv         TNFa/NF-kB Hallmark NES and FDR, pre and post,
                               per cell type - READ THE NOTE BELOW
  pretx_inflammatory_report.txt  the three sections written up, each with the
                               table it was rendered from named in its SOURCES

nfkb_pre_vs_post.csv IS NOT THE TABLE THE PAPER QUOTES.

  Two NF-kB enrichment runs exist and both are kept. They differ because they
  are computed on different matrices, not because one is a corrected copy of
  the other:

    live      this file. Computed by this module from
              {live}
              which was produced on the .X of the input objects. Eight of those
              carry a double normalisation dated 2026-07-30 (see the data
              audit, FINDINGS.md sections 1 and 7).

    adopted   {adopted_name}, written by
              13_R1.8_Neutrophil_Rebuilt_Recompute. Same analysis, recomputed
              on sound input. This is the one Supplementary Table S10, the
              Results, the response letter and section 3 of the report beside
              this file all quote.

  So the numbers here will not match Supplementary Table S10, and are not meant
  to. Measured 2026-09-10 over the 52 paired values ST10's per-cell columns and
  this file have in common - 13 cell types x {{pre, post}} x {{NES, FDR}} - 43
  of the 52 differ, and MoMac pre-treatment differs in sign: here NES +1.055
  with FDR 0.472, in ST10 NES -1.001 with FDR 0.893.

  This file is kept, unaltered, because it is what this module computed and
  because one file name has to keep meaning one set of contents. Regenerating
  it from the adopted table would put a second set of contents behind a name
  that has already been cited, which is the defect this project has already
  been bitten by twice. If you want the paper's numbers, read Table S10 or
  {adopted_name}.
"""


def write_outputs_readme():
    """Write outputs/README.txt, naming the two NF-kB lineages beside the CSV."""
    text = OUTPUTS_README.format(
        live=RECOMPUTE_GSEA.name + "/  (12_R1.8_DEG_Recompute)",
        adopted_name=ADOPTED_GSEA.name)
    (OUT / README_NAME).write_text(text, encoding="utf-8")
    return OUT / README_NAME


SECTION_SOURCE = {
    1: "   Source: this module's own live run, "
       "outputs/pretx_state_abundance.csv",
    2: "   Source: this module's own live run, "
       "outputs/pretx_signature_scores.csv",
    3: "   Source: the ADOPTED table 13_R1.8_Neutrophil_Rebuilt_Recompute/\n"
       "           outputs/nfkb_per_celltype_sound13.csv (method = ttest),\n"
       "           NOT this module's own outputs/nfkb_pre_vs_post.csv",
}

ADOPTED_SOURCES = [
    "   Sections 1 and 2 are this module's own live run, read from its own",
    "   outputs/pretx_state_abundance.csv and outputs/pretx_signature_scores.csv.",
    "",
    "   Section 3, and the NF-kB sentence in the CONCLUSION that repeats it, are",
    "   read from the ADOPTED sound-input table",
    "",
    "       13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/"
    "nfkb_per_celltype_sound13.csv",
    "",
    "   and NOT from this module's own",
    "   outputs/nfkb_pre_vs_post.csv. That file is the live run on the doubly-",
    "   normalised matrix of 00_Data_Audit/FINDINGS.md sections 1 and 7. It is",
    "   unchanged and still sits beside this report; it is simply not the table",
    "   the counts below are taken from.",
    "",
    "   The adopted table is the one Supplementary Figure S10 panel C is drawn",
    "   from, so this report and that panel report the same numbers.",
]


def render_report(frac, score, n, sources=None):
    """
    The report text, from frames that are already prepared.

    `n` is the per-cell-type NF-kB frame for section 3, already restricted to
    cell types that have both phases and already sorted by post NES descending
    - that is, exactly what load_adopted_shift() returns.

    `sources` is the list of lines naming the table behind each section. With
    sources=None nothing is inserted and the output is the original layout
    unchanged, which is what makes this refactor checkable against the report
    that shipped.
    """
    L = ["IS THE INFLAMMATORY PROGRAMME ALREADY PRESENT BEFORE TREATMENT?",
         "Reviewer 2 point R2.1", "=" * 92, ""]
    if sources:
        L.append("SOURCES")
        L.append("-" * 92)
        L.extend(sources)
        L.append("")

    L.append("1. ABUNDANCE OF THE IL-1b+ INFLAMMATORY STATE "
             "(% of monocytes/macrophages)")
    L.append("-" * 92)
    if sources:
        L.append(SECTION_SOURCE[1])
        L.append("")
    for g in GROUPS:
        v = frac.loc[frac["group4"] == g, "fraction"].values
        L.append(f"   {g:<9} n={len(v)}  mean = {v.mean() * 100:5.2f}%  "
                 f"median = {np.median(v) * 100:5.2f}%")
    L.append("")
    for a, b, lab in (("Pre-NR", "Pre-R", "pre-treatment NR vs R"),
                      ("Post-NR", "Post-R", "post-treatment NR vs R"),
                      ("Post-NR", "Pre-NR", "post vs pre within NR"),
                      ("Post-R", "Pre-R", "post vs pre within R")):
        p, r = mw(frac.loc[frac["group4"] == a, "fraction"].values,
                  frac.loc[frac["group4"] == b, "fraction"].values)
        L.append(f"   {lab:<26} P = {p:.4f}   r = {r:+.3f}")
    L.append("")

    L.append("2. IL-1b+ STATE SIGNATURE SCORE ACROSS ALL MONOCYTES/MACROPHAGES")
    L.append("-" * 92)
    if sources:
        L.append(SECTION_SOURCE[2])
        L.append("")
    L.append(f"   Signature: top 50 markers of {TARGET} within this object")
    for g in GROUPS:
        v = score.loc[score["group4"] == g, "score"].values
        L.append(f"   {g:<9} n={len(v)}  mean score = {v.mean():+.4f}")
    L.append("")
    for a, b, lab in (("Pre-NR", "Pre-R", "pre-treatment NR vs R"),
                      ("Post-NR", "Post-R", "post-treatment NR vs R")):
        p, r = mw(score.loc[score["group4"] == a, "score"].values,
                  score.loc[score["group4"] == b, "score"].values)
        L.append(f"   {lab:<26} P = {p:.4f}   r = {r:+.3f}")
    L.append("")

    L.append("3. TNFa/NF-kB HALLMARK ENRICHMENT, PRE VERSUS POST "
             "(positive = enriched in non-responders)")
    L.append("-" * 92)
    if sources:
        L.append(SECTION_SOURCE[3])
        L.append("")
    L.append(f"   {'cell type':<22}{'pre NES':>10}{'pre q':>9}"
             f"{'post NES':>10}{'post q':>9}{'change':>10}")
    for _, r in n.iterrows():
        L.append(f"   {r['cell_type']:<22}{r['pre_nes']:>10.3f}{r['pre_fdr']:>9.3f}"
                 f"{r['post_nes']:>10.3f}{r['post_fdr']:>9.3f}"
                 f"{r['post_nes'] - r['pre_nes']:>+10.3f}")
    L.append("")
    n_pre = int((n["pre_nes"] > 0).sum())
    n_post = int((n["post_nes"] > 0).sum())
    top = n.iloc[0]
    q_pre = int(((n["pre_nes"] > 0) & (n["pre_fdr"] < 0.05)).sum())
    q_post = int(((n["post_nes"] > 0) & (n["post_fdr"] < 0.05)).sum())
    L.append(f"   Positive enrichment in {n_pre}/{len(n)} cell types before treatment "
             f"({q_pre} at q < 0.05)")
    L.append(f"   and in {n_post}/{len(n)} after ({q_post} at q < 0.05).")
    L.append(f"   The strongest after treatment is {top['cell_type']} "
             f"(NES {top['post_nes']:+.3f}, q = {top['post_fdr']:.3f}); before")
    L.append(f"   treatment it is NES {top['pre_nes']:+.3f} (q = {top['pre_fdr']:.3f}).")
    L.append("")
    L.append("CONCLUSION")
    L.append("-" * 92)
    pre_p, _ = mw(frac.loc[frac["group4"] == "Pre-NR", "fraction"].values,
                  frac.loc[frac["group4"] == "Pre-R", "fraction"].values)
    post_p, _ = mw(frac.loc[frac["group4"] == "Post-NR", "fraction"].values,
                   frac.loc[frac["group4"] == "Post-R", "fraction"].values)
    # The abundance difference does not reach significance at either timepoint.
    # Saying the separation "appears" after treatment overstates P = 0.33; what
    # changes is the effect size, from r = +0.13 to r = +0.40.
    L.append(f"   The IL-1b+ state does not distinguish responders from non-responders")
    L.append(f"   before treatment (P = {pre_p:.3f}); after treatment it separates in")
    L.append(f"   the direction of non-response without reaching significance at these")
    L.append(f"   group sizes (P = {post_p:.3f}). NF-kB enrichment follows the same")
    L.append(f"   timing: positive in {n_pre} of {len(n)} cell types before treatment and")
    L.append(f"   {q_pre} at q < 0.05, against {n_post} of {len(n)} and {q_post} after, with")
    L.append("   monocytes/macrophages the strongest")
    L.append("   population after treatment. The programme is therefore best described")
    L.append("   as emerging on treatment rather than as a pre-existing feature of the")
    L.append("   non-responder microenvironment.")

    return "\n".join(L)


def load_adopted_shift():
    """
    Pre and post NES and FDR per cell type from the adopted sound table, in the
    shape, column names and sort order section 3 of the report expects.

    The sort is by post NES, descending; it is the order the section's rows are
    printed in, and it is recomputed rather than carried over, because the
    adopted table may order two cell types differently.
    """
    if not ADOPTED_GSEA.exists():
        raise SystemExit(
            f"the adopted NF-kB table is missing: {ADOPTED_GSEA} - it is "
            "written by 13_R1.8_Neutrophil_Rebuilt_Recompute")
    a = pd.read_csv(ADOPTED_GSEA)
    a = a[a["method"] == "ttest"]
    nes = a.pivot(index="cell_type", columns="phase", values="nes")
    nes.columns = [f"{c}_nes" for c in nes.columns]
    fdr = a.pivot(index="cell_type", columns="phase", values="fdr_q")
    fdr.columns = [f"{c}_fdr" for c in fdr.columns]
    out = nes.join(fdr).reset_index()
    return (out.dropna(subset=["pre_nes", "post_nes"])
            .sort_values("post_nes", ascending=False))


def rewrite_report_from_adopted():
    """
    Rebuild ONLY pretx_inflammatory_report.txt, from tables already on disk.

    main() rereads MoMac.h5ad and rewrites five outputs. Putting a second,
    different set of contents behind an already-cited filename is the
    nfkb_rankings_13types_mast.csv defect this project already carries one
    instance of, so the report is regenerated from the tables, without
    re-running the analysis.

    Sections 1 and 2 are read back from this module's own deposited CSVs, which
    are the frames main() computed those sections from - same columns, same
    rows, so the same numbers. Section 3 comes from the adopted table. No other
    file in outputs/ is opened for writing.
    """
    frac = pd.read_csv(OUT / "pretx_state_abundance.csv")
    score = pd.read_csv(OUT / "pretx_signature_scores.csv")
    report = render_report(frac, score, load_adopted_shift(),
                           sources=ADOPTED_SOURCES)
    (OUT / "pretx_inflammatory_report.txt").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--rewrite-report-from-adopted", action="store_true",
        help="rewrite pretx_inflammatory_report.txt only, from tables already "
             "on disk: sections 1 and 2 from this module's own deposited CSVs, "
             "section 3 from the adopted sound table. No analysis table is "
             "rewritten and MoMac.h5ad is never opened.")
    ap.add_argument(
        "--write-outputs-readme", action="store_true",
        help="write outputs/README.txt only, from this module's own path "
             "constants. No table is rewritten, no analysis runs and "
             "MoMac.h5ad is never opened.")
    args = ap.parse_args()
    if args.rewrite_report_from_adopted:
        rewrite_report_from_adopted()
    elif args.write_outputs_readme:
        print(f"wrote {write_outputs_readme()}")
    else:
        main()
