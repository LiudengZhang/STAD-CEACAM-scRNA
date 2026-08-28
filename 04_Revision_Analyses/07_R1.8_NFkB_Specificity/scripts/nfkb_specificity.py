"""
WP8 / Reviewer 1 point R1.8.

  "The statement that 'every population exhibited elevated NF-kB pathway
   activity' is unusually broad and may reflect a common inflammatory state,
   sample composition, or an analytical artifact. Transcriptomic enrichment of
   the Hallmark TNFa/NF-kB gene set does not necessarily demonstrate biochemical
   activation... Moreover, IL-1 expression and other inflammatory genes should
   also be evaluated in the epithelial compartment."

Three things are produced.

  1. The quantitative version of "every population": NES, nominal P and FDR for
     HALLMARK TNFa signalling via NF-kB in each of the 13 major cell types, so
     the claim can be stated with numbers instead of the word "every".
  2. A specificity control: where TNFa/NF-kB ranks among all 50 Hallmark sets in
     each cell type, and which other sets outrank it. If generic inflammatory
     sets rank alongside it, that is reported rather than hidden.
  3. Epithelial inflammatory cytokine expression (IL1A, IL1B, IL6, TNF, CXCL8)
     by treatment phase and response, to test whether epithelium contributes to
     the inflammatory programme or only responds to it.

No protein-level evidence (nuclear p65, phospho-p65, IkB degradation) exists for
this cohort and none is generated here; that limitation is stated in the
manuscript rather than worked around.

Inputs : Round_5/02_Preparation_for_Panels/GSEA/{pre,post}/*_gsea_hallmark.csv
         Round_5/01_Raw_Inputs/01_H5AD/Epithelial.h5ad, MoMac.h5ad
Outputs: nfkb_per_celltype.csv, hallmark_specificity.csv,
         epithelial_cytokines.csv, nfkb_specificity_report.txt,
         panels S9_E, S9_F (the regulon panel is S11_C)
"""

from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import PREPARATION, EPITHELIAL_H5AD, MOMAC_H5AD, REVISED_PANELS  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S9 = REVISED_PANELS / "Supplementary_New" / "S9_Mechanism_Specificity"
GSEA = PREPARATION / "GSEA"

SCALE, CM, DPI = 4, 1 / 2.54, 300
MIN_CELLS = 20
NFKB_TERM = "TNF-alpha Signaling via NF-kB"
CYTOKINES = ["IL1A", "IL1B", "IL6", "TNF", "CXCL8"]
COLOR_R, COLOR_NR = "#2166AC", "#B2182B"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def _find_nfkb_term(terms):
    """The Hallmark term is spelled differently across GSEApy versions."""
    for t in terms:
        s = str(t).lower().replace("_", " ").replace("-", "")
        if "nf" in s and "kb" in s and ("tnf" in s or "tnfa" in s):
            return t
    return None


def load_gsea():
    """
    Read the per-cell-type Hallmark tables and express every NES on a single
    "non-responder versus responder" axis, so a positive value always means
    enriched in non-responders.

    Sign convention, established empirically rather than from the comment in
    run_gsea_from_mast.py, which is wrong. In DEG/post/MoMac_mast_deg.csv the
    genes IL1B, IL1A, IL6, TNF, CXCL8 and NFKBIA all carry POSITIVE
    logfoldchanges, and all six are independently confirmed here to be higher in
    post-treatment non-responders. MAST logfoldchanges is therefore already
    NR-relative-to-R. run_gsea_from_mast.py negates it before ranking, so in the
    saved tables a NEGATIVE NES means enriched in non-responders. The sign is
    flipped back here.
    """
    rows = []
    for phase in ("pre", "post"):
        for f in sorted((GSEA / phase).glob("*_gsea_hallmark.csv")):
            cell = f.name.replace("_gsea_hallmark.csv", "")
            d = pd.read_csv(f)
            term = _find_nfkb_term(d["Term"])
            if term is None:
                continue
            d["NES_NRvsR"] = -d["NES"]
            # Rank by enrichment in non-responders, most enriched first.
            d = d.sort_values("NES_NRvsR", ascending=False).reset_index(drop=True)
            d["rank"] = np.arange(1, len(d) + 1)
            r = d[d["Term"] == term].iloc[0]
            outranking = d[d["rank"] < r["rank"]]["Term"].tolist()
            rows.append(dict(
                phase=phase, cell_type=cell, term=term,
                nes=float(r["NES_NRvsR"]), nes_as_stored=float(r["NES"]),
                nom_p=float(r["NOM p-val"]),
                fdr_q=float(r["FDR q-val"]), rank=int(r["rank"]),
                n_sets=len(d),
                sets_outranking_it="; ".join(outranking[:5]),
            ))
    return pd.DataFrame(rows)


def load_ttest_rankings():
    """
    The NES values that Figure 5H actually plots. The Methods describe Figure
    5H/I as a Welch t-test analysis, and nfkb_rankings_13types_mast.csv is the
    file the radar panel reads - but that file is byte-identical to
    nfkb_rankings_13types_ttest_reproduce.csv, i.e. despite its name it holds
    the t-test values, not the MAST ones. The MAST values live in
    nfkb_rankings_13types_mast_BACKUP_mast_method.csv. The filename is a
    reproducibility hazard and is corrected in the code release.
    """
    f = GSEA / "nfkb_rankings_13types_ttest_reproduce.csv"
    d = pd.read_csv(f)
    return d[["cell_type", "pre_nes", "post_nes"]].rename(
        columns={"pre_nes": "ttest_pre_nes", "post_nes": "ttest_post_nes"})


def sample_means(path, gene, phase, group_col):
    ad = sc.read_h5ad(path)
    ad = ad[ad.obs["Sample site"] == "Stomach"]
    ad = ad[ad.obs["Treatment phase"] == phase]
    ad = ad[ad.obs[group_col].isin(["Responsed", "No-response"])].copy()
    src = ad.raw if (ad.raw is not None and gene in ad.raw.var_names) else ad
    if gene not in src.var_names:
        return None, None
    i = list(src.var_names).index(gene)
    x = src.X[:, i]
    x = x.toarray().flatten() if hasattr(x, "toarray") else np.asarray(x).flatten()
    df = pd.DataFrame({"sample": ad.obs["sample"].astype(str).values,
                       "group": ad.obs[group_col].astype(str).values, "value": x})
    keep = df.groupby("sample").size()
    df = df[df["sample"].isin(keep[keep >= MIN_CELLS].index)]
    agg = df.groupby(["sample", "group"], observed=True)["value"].mean().reset_index()
    return (agg.loc[agg["group"] == "No-response", "value"].values,
            agg.loc[agg["group"] == "Responsed", "value"].values)


def main():
    # ------------------------------------------- 1 & 2: NF-kB per cell type
    g = load_gsea()
    g.to_csv(OUT / "nfkb_per_celltype.csv", index=False)

    spec = g[["phase", "cell_type", "rank", "n_sets", "nes", "fdr_q",
              "sets_outranking_it"]].copy()
    spec.to_csv(OUT / "hallmark_specificity.csv", index=False)

    L = ["NF-kB SPECIFICITY AND THE EPITHELIAL COMPARTMENT - Reviewer 1 point R1.8",
         "=" * 96, ""]
    L.append(f"Hallmark term used: {g['term'].iloc[0]}")
    L.append("")
    for phase in ("pre", "post"):
        sub = g[g["phase"] == phase].sort_values("nes", ascending=False)
        n_pos = int((sub["nes"] > 0).sum())
        n_fdr25 = int(((sub["nes"] > 0) & (sub["fdr_q"] < 0.25)).sum())
        n_fdr05 = int(((sub["nes"] > 0) & (sub["fdr_q"] < 0.05)).sum())
        L.append(f"[{phase.upper()}-TREATMENT, non-responders vs responders]")
        L.append("-" * 96)
        L.append(f"  {'cell type':<22}{'NES':>8}{'nom P':>10}{'FDR q':>10}"
                 f"{'rank of 50':>12}")
        for _, r in sub.iterrows():
            L.append(f"  {r['cell_type']:<22}{r['nes']:>8.3f}{r['nom_p']:>10.4f}"
                     f"{r['fdr_q']:>10.4f}{r['rank']:>8} / {r['n_sets']}")
        L.append("")
        L.append(f"  Positive NES in {n_pos}/{len(sub)} cell types; "
                 f"FDR q < 0.25 in {n_fdr25}; FDR q < 0.05 in {n_fdr05}.")
        L.append("")
    L.append("CONCORDANCE WITH THE INDEPENDENT DE METHOD USED FOR FIGURE 5H")
    L.append("-" * 96)
    tt = load_ttest_rankings()
    cmp_df = (g[g["phase"] == "post"][["cell_type", "nes", "fdr_q"]]
              .rename(columns={"nes": "mast_post_nes"})
              .merge(tt[["cell_type", "ttest_post_nes"]], on="cell_type", how="inner"))
    cmp_df["same_direction"] = (
        np.sign(cmp_df["mast_post_nes"]) == np.sign(cmp_df["ttest_post_nes"]))
    cmp_df.to_csv(OUT / "nfkb_method_concordance.csv", index=False)
    L.append("  Figure 5H is drawn from a Welch t-test analysis; the tables above are")
    L.append("  from the independent MAST hurdle model. Both are expressed as")
    L.append("  non-responder versus responder.")
    L.append(f"  {'cell type':<22}{'MAST NES':>10}{'t-test NES':>12}{'agree':>8}")
    for _, r in cmp_df.sort_values("mast_post_nes", ascending=False).iterrows():
        L.append(f"  {r['cell_type']:<22}{r['mast_post_nes']:>10.3f}"
                 f"{r['ttest_post_nes']:>12.3f}{'yes' if r['same_direction'] else 'NO':>8}")
    agree = int(cmp_df["same_direction"].sum())
    rho = cmp_df["mast_post_nes"].corr(cmp_df["ttest_post_nes"], method="spearman")
    L.append(f"  Direction agrees in {agree}/{len(cmp_df)} cell types; "
             f"Spearman rho = {rho:.3f}.")
    L.append("")
    L.append("  NOTE FOR THE CODE RELEASE: nfkb_rankings_13types_mast.csv is")
    L.append("  byte-identical to nfkb_rankings_13types_ttest_reproduce.csv, i.e. the")
    L.append("  file named '_mast' contains the t-test values. The MAST values are in")
    L.append("  nfkb_rankings_13types_mast_BACKUP_mast_method.csv. Figure 5H reads the")
    L.append("  misnamed file, which is consistent with the Methods (Figure 5H/I is")
    L.append("  described as a Welch t-test analysis) but the filename must be")
    L.append("  corrected before the code is deposited.")
    L.append("")
    L.append("SPECIFICITY: WHAT OUTRANKS TNFa/NF-kB")
    L.append("-" * 96)
    post = g[g["phase"] == "post"].sort_values("rank")
    for _, r in post.iterrows():
        if r["rank"] == 1:
            L.append(f"  {r['cell_type']:<22} top-ranked Hallmark set")
        else:
            L.append(f"  {r['cell_type']:<22} rank {r['rank']}, outranked by: "
                     f"{r['sets_outranking_it']}")
    L.append("")
    L.append("  Revised wording: instead of 'every population exhibited elevated NF-kB")
    L.append("  pathway activity', the manuscript now reports the number of cell types")
    L.append("  with a positive NES and the FDR in each, and describes the finding as a")
    L.append("  coordinated inflammatory transcriptional programme in which TNFa/NF-kB")
    L.append("  is among the most consistently enriched components. Transcriptomic")
    L.append("  enrichment is not equated with biochemical pathway activation; no")
    L.append("  nuclear p65, phospho-p65 or IkB data exist for this cohort.")
    L.append("")

    # -------------------------------------------- 3: epithelial cytokines
    rows = []
    for compartment, path in (("Epithelial", EPITHELIAL_H5AD),
                              ("Monocytes/Macrophages", MOMAC_H5AD)):
        for phase, col in (("Pre", "stomach_pre_grouping"),
                           ("Post", "stomach_post_grouping")):
            for gene in CYTOKINES:
                nr, r = sample_means(path, gene, phase, col)
                if nr is None or len(nr) < 2 or len(r) < 2:
                    continue
                u, p = stats.mannwhitneyu(nr, r, alternative="two-sided")
                rows.append(dict(
                    compartment=compartment, phase=phase, gene=gene,
                    n_NR=len(nr), n_R=len(r),
                    mean_NR=float(nr.mean()), mean_R=float(r.mean()),
                    rank_biserial_r=float(2.0 * u / (len(nr) * len(r)) - 1.0),
                    p_two_tailed=float(p)))
    cyto = pd.DataFrame(rows)
    cyto.to_csv(OUT / "epithelial_cytokines.csv", index=False)

    L.append("EPITHELIAL VERSUS MYELOID INFLAMMATORY CYTOKINE EXPRESSION")
    L.append("-" * 96)
    L.append(f"  {'compartment':<24}{'phase':<7}{'gene':<8}"
             f"{'mean NR':>10}{'mean R':>10}{'P':>10}")
    for _, r in cyto.iterrows():
        L.append(f"  {r['compartment']:<24}{r['phase']:<7}{r['gene']:<8}"
                 f"{r['mean_NR']:>10.4f}{r['mean_R']:>10.4f}{r['p_two_tailed']:>10.4f}")
    L.append("")
    epi_post = cyto[(cyto["compartment"] == "Epithelial") & (cyto["phase"] == "Post")]
    mm_post = cyto[(cyto["compartment"] == "Monocytes/Macrophages")
                   & (cyto["phase"] == "Post")]
    if len(epi_post) and len(mm_post):
        ratio = (epi_post.set_index("gene")["mean_NR"]
                 / mm_post.set_index("gene")["mean_NR"]).dropna()
        L.append("  Post-treatment non-responder epithelial expression as a fraction of")
        L.append("  monocyte/macrophage expression of the same gene:")
        for gene, v in ratio.items():
            L.append(f"     {gene:<8}{v * 100:6.1f}%")

    report = "\n".join(L)
    (OUT / "nfkb_specificity_report.txt").write_text(report, encoding="utf-8")
    print(report)

    _panel_nfkb(g)
    _panel_cytokines(cyto)


def _panel_nfkb(g):
    d = (S9 / "S9_E"); d.mkdir(parents=True, exist_ok=True)
    post = g[g["phase"] == "post"].sort_values("nes")
    fig, ax = plt.subplots(figsize=(7.5 * SCALE * CM, 5.0 * SCALE * CM))
    y = np.arange(len(post))
    colors = ["#B2182B" if q < 0.25 else "#cccccc" for q in post["fdr_q"]]
    ax.barh(y, post["nes"], color=colors, edgecolor="#444444", linewidth=0.5,
            height=0.65)
    # Annotations always sit to the right of the bar's far end, so the ones on
    # negative bars do not run into the cell-type labels on the axis.
    for i, (_, r) in enumerate(post.iterrows()):
        ax.text(max(r["nes"], 0.0) + 0.05, i,
                f"q={r['fdr_q']:.2f}  #{r['rank']}/{r['n_sets']}",
                va="center", ha="left",
                fontsize=4.5 * SCALE, color="#333333")
    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace("_", " ") for c in post["cell_type"]],
                       fontsize=5.5 * SCALE)
    ax.set_xlabel("NES, TNF$\\alpha$ signalling via NF-$\\kappa$B\n"
                  "(post-treatment, non-responders vs responders)",
                  fontsize=6 * SCALE)
    ax.set_xlim(min(0, post["nes"].min()) - 0.3, post["nes"].max() + 0.9)
    ax.tick_params(axis="x", labelsize=5.5 * SCALE, width=0.8, length=3)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="#444444",
                             linewidth=0.5, label=l)
               for c, l in (("#B2182B", "FDR q < 0.25"), ("#cccccc", "FDR q >= 0.25"))]
    ax.legend(handles=handles, loc="lower right", frameon=False,
              fontsize=5.5 * SCALE)
    fig.subplots_adjust(left=0.24, right=0.98, top=0.97, bottom=0.20)
    stem = d / "S9_E_nfkb_per_celltype"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


def _panel_cytokines(cyto):
    d = (S9 / "S9_F"); d.mkdir(parents=True, exist_ok=True)
    post = cyto[cyto["phase"] == "Post"]
    fig, ax = plt.subplots(figsize=(7.0 * SCALE * CM, 4.0 * SCALE * CM))
    genes = CYTOKINES
    x = np.arange(len(genes))
    w = 0.38
    for k, (comp, color) in enumerate((("Epithelial", "#4d4d4d"),
                                       ("Monocytes/Macrophages", "#B2182B"))):
        sub = post[post["compartment"] == comp].set_index("gene")
        vals = [sub["mean_NR"].get(gene, np.nan) for gene in genes]
        ax.bar(x + (k - 0.5) * w, vals, width=w, color=color, alpha=0.85,
               edgecolor="#333333", linewidth=0.5, label=comp)
    ax.set_xticks(x); ax.set_xticklabels(genes, fontsize=6 * SCALE, style="italic")
    ax.set_ylabel("Mean expression,\npost-treatment non-responders",
                  fontsize=6 * SCALE)
    ax.set_yscale("symlog", linthresh=0.01)
    ax.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=5.5 * SCALE, loc="upper left")
    fig.subplots_adjust(left=0.19, right=0.98, top=0.95, bottom=0.16)
    stem = d / "S9_F_epithelial_vs_myeloid_cytokines"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
