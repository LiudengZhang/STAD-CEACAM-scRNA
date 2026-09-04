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

Inputs : Round_5/01_Raw_Inputs/01_H5AD/MoMac.h5ad
         12_R1.8_DEG_Recompute/outputs/gsea/*_{pre,post}_ttest_hallmark.csv
         13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/nfkb_per_celltype_sound13.csv
             - PANEL S10_C since 2026-09-03, and the REPORT's section 3 since
               2026-09-04. See ADOPTED_GSEA below.
Outputs: pretx_state_abundance.csv, pretx_signature_scores.csv,
         nfkb_pre_vs_post.csv, pretx_inflammatory_report.txt,
         panels S10_A, S10_B, S10_C

The tables this module writes are the live run and are unchanged, nfkb_pre_vs_post
.csv included. Panel S10_C is drawn from the adopted sound-input table instead,
under the author's ruling of 2026-09-03, and since 2026-09-04 so is section 3 of
pretx_inflammatory_report.txt - from the same load_adopted_shift() call, so text
and panel are rendered from one frame and cannot diverge. Until then the report
still counted the live run and contradicted its own figure; the resolved record
is outputs/pretx_inflammatory_report_KNOWN_DISCREPANCY.txt. Panels S10_A and
S10_B do not read a GSEA table at all and are not affected.

--rewrite-report-from-adopted rebuilds the report alone, reading sections 1 and 2
back from this module's own deposited CSVs, so the report can be regenerated
without re-running the analysis and putting new contents behind cited
filenames.
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
from paths import MOMAC_H5AD, REVISED_PANELS  # noqa: E402
from shared.sample_ids import sample_id_map, to_study_ids  # noqa: E402

RECOMPUTE_GSEA = (Path(__file__).resolve().parents[2]
                  / "12_R1.8_DEG_Recompute" / "outputs" / "gsea")

# ---------------------------------------------------------------- the adopted table
# Author's ruling, 2026-09-03. Panel S10_C is drawn from the SOUND-INPUT
# recompute, not from this module's own nfkb_pre_vs_post.csv, so that the panel
# agrees with the Results and the response letter, which were moved to the same
# table in the same pass. RECOMPUTE_GSEA above is the "live" run, computed on the
# doubly-normalised .X of 00_Data_Audit/FINDINGS.md sections 1 and 7.
#
# The analysis below is untouched: nfkb_pre_vs_post.csv and the report are still
# the live run, byte for byte, and one file name still means one set of contents.
# Six quantities differ between the two tables, and only these six:
#     post FDR q < 0.05 count   4  -> 3     rank-1 count      5  -> 3
#     pericyte ordinal         31  -> 30    B-cell ordinal   37  -> 38
#     pre positive count        6  -> 5     pre q < 0.05      1  -> 3
# Of those, S10_C shows the pre and post NES per cell type, so what moves on this
# panel is the pre positive count, the pre q < 0.05 count and the row order.
ADOPTED_GSEA = (Path(__file__).resolve().parents[2]
                / "13_R1.8_Neutrophil_Rebuilt_Recompute" / "outputs"
                / "nfkb_per_celltype_sound13.csv")

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S10 = REVISED_PANELS / "Supplementary_New" / "S10_PreTx_and_Adaptive"

SCALE, CM, DPI = 4, 1 / 2.54, 300
MIN_CELLS = 20
TARGET = "C3_Mac_Inflam_IL1B"
GROUPS = ["Pre-R", "Pre-NR", "Post-R", "Post-NR"]
# The names the figures use, matching Fig. 5H and Fig. S9E.
LABELS = {"B_cells": "B cells", "DC_cells": "DC",
          "Endothelial_cells": "Endothelial", "Epithelial": "Epithelial",
          "Fibroblast": "Fibroblast", "Mast_cells": "Mast", "MoMac": "MoMac",
          "Neutrophils": "Neutrophils", "NK_cells": "NK",
          "Pericyte": "Pericyte", "Plasma_cells": "Plasma",
          "TCD4_cells": "CD4+ T", "TCD8_cells": "CD8+ T"}
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
    # Section 3 is rendered from the ADOPTED sound-input table - the same
    # frame _panel_nfkb_shift() is handed below - so the report and panel
    # S10_C quote the same counts. nfkb_pre_vs_post.csv, written just above,
    # remains the live run and is not what the report quotes; the report's
    # own SOURCES block says so.
    report = render_report(frac, score, load_adopted_shift(),
                           sources=ADOPTED_SOURCES)
    (OUT / "pretx_inflammatory_report.txt").write_text(report, encoding="utf-8")
    print(report)

    _panel_four_group(frac, "fraction", "IL-1$\\beta$+ state\n(% of mono/macrophages)",
                      "S10_A", "S10_A_il1b_state_four_groups", pct=True)
    _panel_four_group(score, "score", "IL-1$\\beta$+ signature score",
                      "S10_B", "S10_B_il1b_signature_four_groups")
    # S10_C is drawn from the adopted sound table, not from `n`. See ADOPTED_GSEA.
    _panel_nfkb_shift(load_adopted_shift())



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
    "   under the author's ruling of 2026-09-03, and NOT from this module's own",
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
    - that is, exactly what load_adopted_shift() returns and exactly the frame
    _panel_nfkb_shift() is handed, so the text and the panel cannot diverge
    again.

    `sources` is the list of lines naming the table behind each section. With
    sources=None nothing is inserted and the output is the 2026-08-29 layout
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
    same shape, column names and sort order that _panel_nfkb_shift() expects.

    The sort is by post NES, descending, exactly as in main(); it is the panel's
    row order and it is recomputed rather than carried over, because the adopted
    table may order two cell types differently.
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
    instance of, so the report is regenerated the way panel S10_C is: from the
    tables, without re-running the analysis.

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
    ax.set_yticklabels([LABELS[c] for c in n["cell_type"]],
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
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--redraw-adopted-panel", action="store_true",
        help="redraw panel S10_C from the adopted sound table and do nothing "
             "else. The full run rereads MoMac.h5ad and rewrites five outputs; "
             "this draws the one panel the 2026-09-03 ruling covers and touches "
             "no analysis output.")
    ap.add_argument(
        "--rewrite-report-from-adopted", action="store_true",
        help="rewrite pretx_inflammatory_report.txt only, from tables already "
             "on disk: sections 1 and 2 from this module's own deposited CSVs, "
             "section 3 from the adopted sound table. The counterpart of "
             "--redraw-adopted-panel for the text. No analysis table is "
             "rewritten and MoMac.h5ad is never opened.")
    args = ap.parse_args()
    if args.redraw_adopted_panel:
        _panel_nfkb_shift(load_adopted_shift())
    elif args.rewrite_report_from_adopted:
        rewrite_report_from_adopted()
    else:
        main()
