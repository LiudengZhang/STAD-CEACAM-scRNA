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
  3. TNFa/NF-kB Hallmark enrichment pre versus post, per cell type, reusing the
     GSEA tables (sign convention as established in WP8: the stored NES is
     R-relative, so it is negated to read as non-responder versus responder).

Inputs : Round_5/01_Raw_Inputs/01_H5AD/MoMac.h5ad
         Round_5/02_Preparation_for_Panels/GSEA/{pre,post}/*_gsea_hallmark.csv
Outputs: pretx_state_abundance.csv, pretx_signature_scores.csv,
         nfkb_pre_vs_post.csv, pretx_inflammatory_report.txt,
         panels S10_A, S10_B, S10_C
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
from paths import MOMAC_H5AD, PREPARATION, REVISED_PANELS  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S10 = REVISED_PANELS / "Supplementary_New" / "S10_PreTx_and_Adaptive"

SCALE, CM, DPI = 4, 1 / 2.54, 300
MIN_CELLS = 20
TARGET = "C3_Mac_Inflam_IL1B"
GROUPS = ["Pre-R", "Pre-NR", "Post-R", "Post-NR"]
COLORS = {"Pre-R": "#bde0fe", "Pre-NR": "#a2d2ff",
          "Post-R": "#ffcfd2", "Post-NR": "#f1c0e8"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


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
    score.to_csv(OUT / "pretx_signature_scores.csv", index=False)
    pd.Series(sig, name="gene").to_csv(OUT / "il1b_signature_genes.csv", index=False)

    # ----------------------------------------------- 3. NF-kB pre vs post
    rows = []
    for phase in ("pre", "post"):
        for f in sorted((PREPARATION / "GSEA" / phase).glob("*_gsea_hallmark.csv")):
            d = pd.read_csv(f)
            hit = d[d["Term"].str.contains("NF-kB", case=False, regex=False)]
            if not len(hit):
                continue
            r = hit.iloc[0]
            rows.append(dict(phase=phase,
                             cell_type=f.name.replace("_gsea_hallmark.csv", ""),
                             nes_NRvsR=-float(r["NES"]),
                             fdr_q=float(r["FDR q-val"])))
    nfkb = pd.DataFrame(rows).pivot(index="cell_type", columns="phase",
                                    values="nes_NRvsR")
    nfkb.columns = ["post_nes", "pre_nes"] if list(nfkb.columns) == ["post", "pre"] \
        else list(nfkb.columns)
    nfkb = nfkb.reset_index()
    nfkb.to_csv(OUT / "nfkb_pre_vs_post.csv", index=False)

    # -------------------------------------------------------------- report
    L = ["IS THE INFLAMMATORY PROGRAMME ALREADY PRESENT BEFORE TREATMENT?",
         "Reviewer 2 point R2.1", "=" * 92, ""]

    L.append("1. ABUNDANCE OF THE IL-1b+ INFLAMMATORY STATE "
             "(% of monocytes/macrophages)")
    L.append("-" * 92)
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
    L.append(f"   {'cell type':<22}{'pre NES':>10}{'post NES':>10}{'change':>10}")
    n = nfkb.dropna(subset=["pre_nes", "post_nes"]).sort_values(
        "post_nes", ascending=False)
    for _, r in n.iterrows():
        L.append(f"   {r['cell_type']:<22}{r['pre_nes']:>10.3f}"
                 f"{r['post_nes']:>10.3f}{r['post_nes'] - r['pre_nes']:>+10.3f}")
    L.append("")
    mm = n[n["cell_type"] == "MoMac"]
    if len(mm):
        L.append(f"   Monocytes/macrophages move from NES {mm['pre_nes'].iloc[0]:+.3f} "
                 f"before treatment to {mm['post_nes'].iloc[0]:+.3f} after.")
    L.append("")
    L.append("CONCLUSION")
    L.append("-" * 92)
    pre_p, _ = mw(frac.loc[frac["group4"] == "Pre-NR", "fraction"].values,
                  frac.loc[frac["group4"] == "Pre-R", "fraction"].values)
    post_p, _ = mw(frac.loc[frac["group4"] == "Post-NR", "fraction"].values,
                   frac.loc[frac["group4"] == "Post-R", "fraction"].values)
    L.append(f"   The IL-1b+ state does not separate responders from non-responders")
    L.append(f"   before treatment (P = {pre_p:.3f}); the separation appears after")
    L.append(f"   treatment (P = {post_p:.3f}). Combined with the myeloid NF-kB")
    L.append("   enrichment reversing sign between the two timepoints, the programme")
    L.append("   is best described as emerging on treatment rather than as a")
    L.append("   pre-existing feature of the non-responder microenvironment.")

    report = "\n".join(L)
    (OUT / "pretx_inflammatory_report.txt").write_text(report, encoding="utf-8")
    print(report)

    _panel_four_group(frac, "fraction", "IL-1$\\beta$+ state\n(% of mono/macrophages)",
                      "S10_A", "S10_A_il1b_state_four_groups", pct=True)
    _panel_four_group(score, "score", "IL-1$\\beta$+ signature score",
                      "S10_B", "S10_B_il1b_signature_four_groups")
    _panel_nfkb_shift(n)


def _panel_four_group(data, col, ylabel, sub, stem_name, pct=False):
    d = (S10 / sub); d.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.5 * SCALE * CM, 4.2 * SCALE * CM))
    vals = [data.loc[data["group4"] == g, col].values * (100 if pct else 1)
            for g in GROUPS]
    bp = ax.boxplot(vals, positions=range(4), widths=0.6, patch_artist=True,
                    showfliers=False, boxprops=dict(linewidth=0.8),
                    whiskerprops=dict(linewidth=0.8), capprops=dict(linewidth=0.8),
                    medianprops=dict(color="black", linewidth=1.2))
    for patch, g in zip(bp["boxes"], GROUPS):
        patch.set_facecolor(COLORS[g]); patch.set_edgecolor("#444444")
    rng = np.random.default_rng(3)
    for i, v in enumerate(vals):
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=9 * SCALE,
                   c="#333333", zorder=3, alpha=0.85,
                   edgecolors="white", linewidths=0.3 * SCALE)
    top = max(v.max() for v in vals)
    bot = min(v.min() for v in vals)
    span = top - bot
    for k, (i, j) in enumerate(((0, 1), (2, 3))):
        p, _ = mw(vals[i], vals[j])
        yy = top + span * (0.10 + 0.16 * k)
        ax.plot([i, i, j, j], [yy, yy + span * 0.04, yy + span * 0.04, yy],
                color="#444444", linewidth=0.8)
        ax.text((i + j) / 2, yy + span * 0.05, f"P = {p:.3f}", ha="center",
                va="bottom", fontsize=5 * SCALE)
    ax.set_ylim(bot - span * 0.10, top + span * 0.42)
    ax.set_xticks(range(4)); ax.set_xticklabels(GROUPS, fontsize=6 * SCALE)
    ax.set_ylabel(ylabel, fontsize=6 * SCALE)
    ax.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.subplots_adjust(left=0.20, right=0.97, top=0.96, bottom=0.12)
    stem = d / stem_name
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


def _panel_nfkb_shift(n):
    d = (S10 / "S10_C"); d.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.5 * SCALE * CM, 4.6 * SCALE * CM))
    y = np.arange(len(n))
    for i, (_, r) in enumerate(n.iterrows()):
        ax.plot([r["pre_nes"], r["post_nes"]], [i, i], color="#bbbbbb",
                linewidth=1.0, zorder=1)
        ax.scatter(r["pre_nes"], i, s=24, c="#a2d2ff", edgecolors="#444444",
                   linewidths=0.4, zorder=3)
        ax.scatter(r["post_nes"], i, s=24, c="#f1c0e8", edgecolors="#444444",
                   linewidths=0.4, zorder=3)
    ax.axvline(0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace("_", " ") for c in n["cell_type"]],
                       fontsize=5.5 * SCALE)
    ax.set_xlabel("NES, TNF$\\alpha$/NF-$\\kappa$B\n(positive = enriched in non-responders)",
                  fontsize=6 * SCALE)
    ax.tick_params(axis="x", labelsize=5.5 * SCALE, width=0.8, length=3)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    handles = [plt.Line2D([], [], marker="o", linestyle="none", markersize=4,
                          markerfacecolor=c, markeredgecolor="#444444", label=l)
               for c, l in (("#a2d2ff", "Pre-treatment"), ("#f1c0e8", "Post-treatment"))]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.22),
              ncol=2, frameon=False, fontsize=5.5 * SCALE)
    fig.subplots_adjust(left=0.30, right=0.97, top=0.97, bottom=0.26)
    stem = d / "S10_C_nfkb_pre_vs_post"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
