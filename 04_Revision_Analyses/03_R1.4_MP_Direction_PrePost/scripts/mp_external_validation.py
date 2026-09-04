"""
WP4b / Reviewer 1 point R1.4 - external validation of the metaprograms.

The corrected in-cohort result is that MP4 is lowest in pre-treatment
responders, but the direct responder-versus-non-responder contrast is null at
4 versus 4. Rather than leave the metaprogram claim resting on an underpowered
comparison, the signature is scored in an independent, properly sized cohort:
45 anti-PD-1-treated gastric tumours (PRJEB25780/TIGER) with known response and
BayesPrism-deconvolved epithelial expression, so the score is computed in the
same compartment the metaprograms were derived from.

Scoring is deliberately simple - the mean of within-cohort z-scores of the
signature genes - because a more elaborate score would be harder to defend than
the signal it recovers.

Inputs : ST3_gene_signatures.csv (MP gene lists as published)
         BayesPrism/tiger_bayesprism_epithelial_expression.tsv
         PRJEB25780 metadata (response, ESTIMATE purity)
Outputs: mp_external_validation.csv, mp_external_report.txt, panel S8_D
"""

from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MANUSCRIPT, PREPARATION, TIGER_META, REVISED_PANELS  # noqa: E402

warnings.filterwarnings("ignore")

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S8 = REVISED_PANELS / "Supplementary_New" / "S8_CEACAM_Metaprogram"

SCALE, CM, DPI = 4, 1 / 2.54, 300
COLOR_R, COLOR_NR = "#2166AC", "#B2182B"


def load_signatures():
    st3 = pd.read_csv(MANUSCRIPT / "04_Tables" / "ST3_gene_signatures.csv")
    sigs = {}
    for _, r in st3.iterrows():
        name = str(r["Signature"])
        if "MP" not in name.upper():
            continue
        genes = [g.strip() for g in str(r["Genes"]).split(";") if g.strip()]
        sigs[name] = genes
    return sigs


def main():
    L = ["METAPROGRAM EXTERNAL VALIDATION - Reviewer 1 point R1.4", "=" * 92, ""]

    sigs = load_signatures()
    if not sigs:
        raise SystemExit("No MP signatures found in ST3_gene_signatures.csv")

    expr = pd.read_csv(PREPARATION / "BayesPrism"
                       / "tiger_bayesprism_epithelial_expression.tsv",
                       sep="\t", index_col=0)
    meta = pd.read_csv(TIGER_META, sep="\t").set_index("sample_id")
    common = expr.index.intersection(meta.index)
    expr, meta = expr.loc[common], meta.loc[common]
    meta = meta[meta["response_NR"].isin(["R", "N"])]
    expr = expr.loc[meta.index]

    lx = np.log2(expr + 1)
    z = (lx - lx.mean()) / lx.std(ddof=1).replace(0, np.nan)

    L.append(f"Cohort: PRJEB25780 (TIGER), {len(meta)} anti-PD-1-treated tumours "
             f"({(meta['response_NR'] == 'N').sum()} NR, "
             f"{(meta['response_NR'] == 'R').sum()} R)")
    L.append("Compartment: BayesPrism-deconvolved epithelium")
    L.append("")

    rows, scores = [], {}
    for name, genes in sorted(sigs.items()):
        present = [g for g in genes if g in z.columns]
        if len(present) < 5:
            L.append(f"  {name}: only {len(present)} of {len(genes)} genes present, "
                     "skipped")
            continue
        s = z[present].mean(axis=1)
        scores[name] = s
        nr = s[meta["response_NR"] == "N"]
        r = s[meta["response_NR"] == "R"]
        u, p = stats.mannwhitneyu(nr, r, alternative="two-sided")
        eff = 2.0 * u / (len(nr) * len(r)) - 1.0
        df = pd.DataFrame({"score": s,
                           "is_NR": (meta["response_NR"] == "N").astype(int),
                           "purity": meta["TumorPurity"]}).dropna()
        m = smf.ols("score ~ is_NR + purity", df).fit()
        rows.append(dict(
            program=name, n_genes_used=len(present), n_genes_total=len(genes),
            n_NR=len(nr), n_R=len(r), mean_NR=float(nr.mean()),
            mean_R=float(r.mean()), rank_biserial_r=float(eff),
            p_two_tailed=float(p),
            beta_purity_adjusted=float(m.params["is_NR"]),
            p_purity_adjusted=float(m.pvalues["is_NR"])))

    v = pd.DataFrame(rows)
    v.to_csv(OUT / "mp_external_validation.csv", index=False)
    # The per-sample scores as well as the summary. _panel() draws these 45
    # points per programme, and until 2026-09-01 they existed only inside this
    # function, so the panel could not be redrawn without recomputing the score
    # - which would have put the scoring code in two places. One file, no
    # number changed.
    pd.DataFrame(scores).rename_axis("sample").to_csv(
        OUT / "mp_external_sample_scores.csv")

    L.append(f"  {'program':<26}{'genes':>9}{'mean NR':>10}{'mean R':>10}"
             f"{'r':>8}{'P':>10}{'P adj':>10}")
    for _, x in v.iterrows():
        flag = "  *" if x["p_two_tailed"] < 0.05 else ""
        L.append(f"  {x['program']:<26}"
                 f"{x['n_genes_used']:>5}/{x['n_genes_total']:<3}"
                 f"{x['mean_NR']:>10.3f}{x['mean_R']:>10.3f}"
                 f"{x['rank_biserial_r']:>8.3f}{x['p_two_tailed']:>10.4f}"
                 f"{x['p_purity_adjusted']:>10.4f}{flag}")
    L.append("")
    mp4 = v[v["program"].str.upper().str.contains("MP4")]
    if len(mp4):
        x = mp4.iloc[0]
        direction = "higher" if x["mean_NR"] > x["mean_R"] else "lower"
        L.append(f"  MP4 is {direction} in non-responders in this independent")
        L.append(f"  cohort (two-sided P = {x['p_two_tailed']:.4f}, purity-adjusted "
                 f"P = {x['p_purity_adjusted']:.4f}).")
        L.append("  The corrected in-cohort result predicts 'higher in non-responders'.")

    report = "\n".join(L)
    (OUT / "mp_external_report.txt").write_text(report, encoding="utf-8")
    print(report)
    if scores:
        _panel(v, scores, meta)


def _panel(v, scores, meta):
    d = (S8 / "S8_D"); d.mkdir(parents=True, exist_ok=True)
    keys = sorted(scores)
    fig, axes = plt.subplots(
        1, len(keys), figsize=(max(2.2 * len(keys), 4.0) * SCALE * CM,
                               4.0 * SCALE * CM))
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(0)
    for ax, name in zip(axes, keys):
        s = scores[name]
        r = s[meta["response_NR"] == "R"].values
        nr = s[meta["response_NR"] == "N"].values
        bp = ax.boxplot([r, nr], positions=[0, 1], widths=0.55, patch_artist=True,
                        showfliers=False, boxprops=dict(linewidth=0.8),
                        whiskerprops=dict(linewidth=0.8),
                        capprops=dict(linewidth=0.8),
                        medianprops=dict(color="black", linewidth=1.2))
        bp["boxes"][0].set_facecolor(COLOR_R); bp["boxes"][0].set_alpha(0.55)
        bp["boxes"][1].set_facecolor(COLOR_NR); bp["boxes"][1].set_alpha(0.55)
        for i, (vals, c) in enumerate(((r, COLOR_R), (nr, COLOR_NR))):
            ax.scatter(i + rng.uniform(-0.1, 0.1, len(vals)), vals, s=7 * SCALE,
                       c=c, zorder=3, edgecolors="white", linewidths=0.3 * SCALE)
        p = v.loc[v["program"] == name, "p_two_tailed"]
        short = name.replace(" (stomach)", "").replace("NMF ", "")
        ax.set_title(f"{short}\nP = {p.iloc[0]:.3f}" if len(p) else short,
                     fontsize=5.5 * SCALE)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["R", "NR"], fontsize=5.5 * SCALE)
        ax.tick_params(axis="both", labelsize=5 * SCALE, width=0.8, length=3)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    axes[0].set_ylabel("Signature score\n(PRJEB25780 epithelium)",
                       fontsize=5.5 * SCALE)
    fig.subplots_adjust(left=0.14, right=0.98, top=0.80, bottom=0.12, wspace=0.45)
    stem = d / "S8_D_metaprogram_external_validation"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
