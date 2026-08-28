"""
WP6c / Reviewer 1 point R1.6 - the properly powered test, n = 407.

The Visium cohort has 10 sections and the TIGER cohort 45 tumours. TCGA-STAD
gives 407 tumours with matched deconvolved epithelial expression, which is the
largest independent test available for the question the reviewer is really
asking: is CEACAM5/6 expression associated with less immune infiltrate once
tumour content is accounted for?

Method, chosen to stay consistent with the manuscript rather than to introduce a
new one: ESTIMATE (Yoshihara 2013), the same algorithm already used to derive
tumour purity for the TIGER cohort, applied to the cleaned tumour-only count
matrix. ESTIMATE ranks genes within each sample, so it is insensitive to the
count-versus-FPKM distinction. Epithelial CEACAM5/6 comes from the existing
BayesPrism deconvolution, so the exposure is epithelium-specific and cannot be
driven by CEACAM expression in infiltrating cells.

The critical covariate is tumour purity. Epithelial expression and immune score
are both functions of purity, so an unadjusted correlation between them is
uninterpretable; purity is therefore included in every model.

Inputs : BayesPrism_TCGA/tcga_bulk_counts_tumor_only.tsv (407 tumours)
         BayesPrism_TCGA/tcga_bayesprism_epithelial_expression.tsv
         one raw augmented_star_gene_counts.tsv, for the Ensembl-to-symbol map
Outputs: tcga_estimate_scores.csv, tcga_immune_exclusion.csv,
         tcga_immune_report.txt, panel S11_B
"""

from pathlib import Path
import subprocess
import sys
import tempfile
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import PREPARATION, RAW_INPUTS, REVIEWER_MATERIALS, REVISED_PANELS  # noqa: E402

ABSOLUTE_PURITY = (REVIEWER_MATERIALS / "public_data"
                   / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt")

warnings.filterwarnings("ignore")

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S11 = REVISED_PANELS / "Supplementary_New" / "S11_Affirmative_Analyses"

TCGA = PREPARATION / "BayesPrism_TCGA"
RAWDIR = (RAW_INPUTS / "02_External" / "Bulk" / "TCGA_STAD" / "02_Raw_Data"
          / "Gene_Expression" / "FPKM" / "raw_files")

SCALE, CM, DPI = 4, 1 / 2.54, 300

R_SCRIPT = """
suppressMessages(library(estimate))
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; outdir <- args[2]
gct <- file.path(outdir, "input.gct")
filtered <- file.path(outdir, "filtered.gct")
scores <- file.path(outdir, "estimate_scores.gct")
filterCommonGenes(input.f = infile, output.f = filtered, id = "GeneSymbol")
estimateScore(input.ds = filtered, output.ds = scores, platform = "illumina")
cat("done\\n")
"""


def ensembl_to_symbol():
    """Build the mapping from any one raw GDC count file, which carries both."""
    f = next(RAWDIR.glob("*.augmented_star_gene_counts.tsv"), None)
    if f is None:
        raise SystemExit(f"No raw GDC count file under {RAWDIR}")
    d = pd.read_csv(f, sep="\t", comment="#", low_memory=False)
    d = d[d["gene_id"].astype(str).str.startswith("ENSG")]
    return dict(zip(d["gene_id"], d["gene_name"]))


def run_estimate(expr_symbols):
    """expr_symbols: genes x samples, index = HGNC symbol."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        infile = td / "expr.txt"
        expr_symbols.to_csv(infile, sep="\t")
        rs = td / "run.R"
        rs.write_text(R_SCRIPT)
        p = subprocess.run(
            ["conda", "run", "-n", "r_demo", "Rscript", str(rs), str(infile), str(td)],
            capture_output=True, text=True)
        out = td / "estimate_scores.gct"
        if not out.exists():
            raise SystemExit(
                "ESTIMATE did not produce output.\n"
                f"stdout:\n{p.stdout[-2000:]}\nstderr:\n{p.stderr[-2000:]}")
        g = pd.read_csv(out, sep="\t", skiprows=2, index_col=0)
        g = g.drop(columns=[c for c in g.columns if c.lower() == "description"])
        return g.T


def main():
    L = ["TCGA-STAD IMMUNE INFILTRATION versus EPITHELIAL CEACAM5/6",
         "Reviewer 1 point R1.6, n = 407", "=" * 92, ""]

    # ------------------------------------------------------------- expression
    bulk = pd.read_csv(TCGA / "tcga_bulk_counts_tumor_only.tsv", sep="\t",
                       index_col=0)
    bulk = bulk[~bulk.index.astype(str).str.startswith("N_")]
    mapping = ensembl_to_symbol()
    bulk.index = [mapping.get(i, None) for i in bulk.index]
    bulk = bulk[bulk.index.notna()]
    bulk = bulk.groupby(level=0).max()  # one row per symbol
    L.append(f"Bulk matrix after symbol mapping: "
             f"{bulk.shape[0]:,} genes x {bulk.shape[1]} tumours")

    est = run_estimate(bulk)
    est.index = [str(i).replace(".", "-") for i in est.index]
    est.to_csv(OUT / "tcga_estimate_scores.csv")
    L.append(f"ESTIMATE scores computed for {len(est)} tumours "
             f"({', '.join(est.columns[:4])})")
    L.append("")

    # --------------------------------------------------- deconvolved CEACAM
    epi = pd.read_csv(TCGA / "tcga_bayesprism_epithelial_expression.tsv",
                      sep="\t", index_col=0)
    epi.index = [str(i)[:12] for i in epi.index]
    est.index = [str(i)[:12] for i in est.index]

    common = epi.index.intersection(est.index)
    d = pd.DataFrame(index=common)
    d["CEACAM"] = np.log2(epi.loc[common, ["CEACAM5", "CEACAM6"]].mean(axis=1) + 1)
    d["CEACAM5"] = np.log2(epi.loc[common, "CEACAM5"] + 1)
    d["CEACAM6"] = np.log2(epi.loc[common, "CEACAM6"] + 1)
    for c in est.columns:
        d[c] = est.loc[common, c]
    # Purity must come from OUTSIDE the expression data. ESTIMATE derives its
    # purity from ESTIMATEScore = ImmuneScore + StromalScore, so adjusting an
    # ESTIMATE immune score for ESTIMATE purity regresses the outcome on a
    # deterministic function of itself. TCGA PanCanAtlas ABSOLUTE purity is
    # called from DNA copy number and is independent of RNA, so it is used here.
    absolute = pd.read_csv(ABSOLUTE_PURITY, sep="\t")
    absolute["patient"] = absolute["array"].astype(str).str[:12]
    ab = (absolute.dropna(subset=["purity"])
          .groupby("patient")["purity"].mean())
    d["ABSOLUTE_purity"] = ab.reindex(d.index)
    purity_col = "ABSOLUTE_purity"

    n_with = int(d[purity_col].notna().sum())
    L.append(f"Tumours with both deconvolved epithelium and ESTIMATE: {len(d)}")
    L.append(f"Of these, with an ABSOLUTE purity call: {n_with}")
    L.append("Purity covariate: TCGA PanCanAtlas ABSOLUTE (DNA-based), NOT the")
    L.append("ESTIMATE-derived purity, which would be circular against ImmuneScore.")
    L.append("")

    rows = []
    for exposure in ("CEACAM", "CEACAM5", "CEACAM6"):
        for target in ("ImmuneScore", "StromalScore"):
            if target not in d.columns:
                continue
            sub = d.dropna(subset=[exposure, target, purity_col])
            rho, p_raw = stats.spearmanr(sub[exposure], sub[target])
            m = smf.ols(f"{target} ~ {exposure} + {purity_col}", sub).fit()
            ci = m.conf_int().loc[exposure]
            rows.append(dict(
                exposure=exposure, outcome=target, n=len(sub),
                spearman_rho=float(rho), spearman_p=float(p_raw),
                beta_purity_adjusted=float(m.params[exposure]),
                p_purity_adjusted=float(m.pvalues[exposure]),
                ci_lo=float(ci[0]), ci_hi=float(ci[1])))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "tcga_immune_exclusion.csv", index=False)

    L.append("ASSOCIATION WITH INFILTRATION (ESTIMATE), purity-adjusted")
    L.append("-" * 92)
    L.append(f"  {'exposure':<10}{'outcome':<14}{'n':>6}{'rho':>9}{'P raw':>11}"
             f"{'beta adj':>11}{'P adj':>11}")
    for _, r in t.iterrows():
        flag = "  *" if r["p_purity_adjusted"] < 0.05 else ""
        L.append(f"  {r['exposure']:<10}{r['outcome']:<14}{r['n']:>6}"
                 f"{r['spearman_rho']:>9.3f}{r['spearman_p']:>11.3g}"
                 f"{r['beta_purity_adjusted']:>11.2f}"
                 f"{r['p_purity_adjusted']:>11.3g}{flag}")
    L.append("")

    imm = t[(t["exposure"] == "CEACAM") & (t["outcome"] == "ImmuneScore")]
    if len(imm):
        r = imm.iloc[0]
        direction = ("lower" if r["beta_purity_adjusted"] < 0 else "higher")
        L.append("INTERPRETATION")
        L.append("-" * 92)
        L.append(f"  At n = {r['n']}, tumours with higher deconvolved epithelial")
        L.append(f"  CEACAM5/6 show {direction} ESTIMATE immune score after adjustment")
        L.append(f"  for tumour purity (beta = {r['beta_purity_adjusted']:.1f}, "
                 f"P = {r['p_purity_adjusted']:.3g}).")
        L.append("  This is a whole-tumour measure. It is compatible with the spatial")
        L.append("  finding only if exclusion is local rather than global; the two")
        L.append("  analyses answer different questions and are reported as such.")

    report = "\n".join(L)
    (OUT / "tcga_immune_report.txt").write_text(report, encoding="utf-8")
    print(report)
    _panel(d, t, purity_col)


def _panel(d, t, purity_col):
    dd = (S11 / "S11_B"); dd.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(8.0 * SCALE * CM, 4.0 * SCALE * CM))

    ax = axes[0]
    if "ImmuneScore" in d.columns:
        ax.scatter(d["CEACAM"], d["ImmuneScore"], s=4 * SCALE, c="#4d4d4d",
                   alpha=0.45, edgecolors="none")
        z = np.polyfit(d["CEACAM"].dropna(),
                       d.loc[d["CEACAM"].notna(), "ImmuneScore"], 1)
        xs = np.linspace(d["CEACAM"].min(), d["CEACAM"].max(), 50)
        ax.plot(xs, np.polyval(z, xs), color="#B2182B", linewidth=1.2)
        ax.set_xlabel("Epithelial CEACAM5/6\n(log2, deconvolved)", fontsize=6 * SCALE)
        ax.set_ylabel("ESTIMATE immune score", fontsize=6 * SCALE)
        row = t[(t["exposure"] == "CEACAM") & (t["outcome"] == "ImmuneScore")]
        if len(row):
            ax.set_title(f"TCGA-STAD, n = {int(row['n'].iloc[0])}\n"
                         f"purity-adjusted P = {row['p_purity_adjusted'].iloc[0]:.3g}",
                         fontsize=6 * SCALE)

    ax = axes[1]
    sub = t[t["outcome"] == "ImmuneScore"]
    y = np.arange(len(sub))
    cols = ["#B2182B" if v < 0 else "#2166AC" for v in sub["beta_purity_adjusted"]]
    ax.barh(y, sub["beta_purity_adjusted"], color=cols, edgecolor="#333",
            linewidth=0.4, height=0.6)
    ax.errorbar(sub["beta_purity_adjusted"], y,
                xerr=[sub["beta_purity_adjusted"] - sub["ci_lo"],
                      sub["ci_hi"] - sub["beta_purity_adjusted"]],
                fmt="none", ecolor="#333", elinewidth=0.8, capsize=2)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["exposure"], fontsize=5.5 * SCALE)
    ax.axvline(0, color="#666", linewidth=0.8)
    ax.set_xlabel("Immune score change\nper log2 CEACAM\n(purity-adjusted)",
                  fontsize=5.5 * SCALE)

    for ax in axes:
        ax.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    fig.subplots_adjust(left=0.14, right=0.98, top=0.84, bottom=0.32, wspace=0.5)
    stem = dd / "S11_B_tcga_immune_exclusion"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
