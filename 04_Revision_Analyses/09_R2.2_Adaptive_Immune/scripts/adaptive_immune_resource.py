"""
WP10 / Reviewer 2 point R2.2.

  "...the addition of supplementary analyses comparing transcriptional programs
   across other immune cell populations - particularly adaptive immune cells -
   would further enhance the value of this resource for the broader cancer
   immunology community."

Produces the resource the reviewer asks for, at two levels:

  1. Composition. For CD4+ T, CD8+ T, NK and B cells, the sample-level fraction
     of every annotated minor state across the four treatment/response groups,
     with two-sided tests at both timepoints and Benjamini-Hochberg correction
     within each lineage.
  2. Programme. The Hallmark gene sets most strongly enriched in non-responders
     in each adaptive lineage, pre and post treatment, from the module 12
     recompute. Module 12 ranks non-responders against responders, so a
     positive NES already means enriched in non-responders and no sign flip is
     applied here.

Inputs : Round_5/01_Raw_Inputs/01_H5AD/{TCD4,TCD8,NK_cells,B_cells}.h5ad
         12_R1.8_DEG_Recompute/outputs/gsea/*_{pre,post}_ttest_hallmark.csv
Outputs: adaptive_state_fractions.csv, adaptive_state_tests.csv,
         adaptive_hallmark_top.csv, adaptive_immune_report.txt,
         panels S10_D and S10_E
"""

from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import (TCD4_H5AD, TCD8_H5AD, NK_CELLS_H5AD, B_CELLS_H5AD,
                   REVISED_PANELS)  # noqa: E402
from shared.sample_ids import sample_id_map, to_study_ids  # noqa: E402

RECOMPUTE_GSEA = (Path(__file__).resolve().parents[2]
                  / "12_R1.8_DEG_Recompute" / "outputs" / "gsea")

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S10 = REVISED_PANELS / "Supplementary_New" / "S10_PreTx_and_Adaptive"

SCALE, CM, DPI = 4, 1 / 2.54, 300
MIN_CELLS = 20
GROUPS = ["Pre-R", "Pre-NR", "Post-R", "Post-NR"]
COLORS = {"Pre-R": "#bde0fe", "Pre-NR": "#a2d2ff",
          "Post-R": "#ffcfd2", "Post-NR": "#f1c0e8"}

LINEAGES = {
    "CD4+ T cells": (TCD4_H5AD, "TCD4_cells"),
    "CD8+ T cells": (TCD8_H5AD, "TCD8_cells"),
    "NK cells": (NK_CELLS_H5AD, "NK_cells"),
    "B cells": (B_CELLS_H5AD, "B_cells"),
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def four_group_label(row):
    if row["Treatment phase"] == "Pre":
        return {"Responsed": "Pre-R", "No-response": "Pre-NR"}.get(
            row["stomach_pre_grouping"])
    return {"Responsed": "Post-R", "No-response": "Post-NR"}.get(
        row["stomach_post_grouping"])


def main():
    frac_rows, test_rows = [], []

    for lineage, (path, _) in LINEAGES.items():
        ad = sc.read_h5ad(path)
        # Each lineage object is relabelled from its own crosswalk. They are not
        # pooled into one: NK_cells.h5ad spells eight lymph-node specimens
        # differently from the other objects, and while none of those survive
        # the stomach filter below, a merged crosswalk would refuse to build.
        ad = ad[ad.obs["Sample site"] == "Stomach"].copy()
        ids = sample_id_map(ad.obs)
        ad.obs["group4"] = ad.obs.apply(four_group_label, axis=1)
        ad = ad[ad.obs["group4"].isin(GROUPS)].copy()
        obs = ad.obs[["sample", "group4", "minor_cell_state"]].copy()
        obs["sample"] = obs["sample"].astype(str)
        obs["minor_cell_state"] = obs["minor_cell_state"].astype(str)

        keep = obs.groupby("sample").size()
        obs = obs[obs["sample"].isin(keep[keep >= MIN_CELLS].index)]

        tab = (obs.groupby(["sample", "group4"], observed=True)["minor_cell_state"]
               .value_counts(normalize=True).rename("fraction").reset_index())
        tab["lineage"] = lineage
        frac_rows.append((tab, ids))

        states = sorted(obs["minor_cell_state"].unique())
        for state in states:
            sub = tab[tab["minor_cell_state"] == state]
            # Samples where the state is absent contribute a zero, otherwise the
            # comparison silently drops them and inflates the mean.
            filled = {}
            for g in GROUPS:
                samples = obs.loc[obs["group4"] == g, "sample"].unique()
                s = sub.set_index("sample")["fraction"]
                filled[g] = np.array([s.get(x, 0.0) for x in samples])
            for a, b, phase in (("Pre-NR", "Pre-R", "Pre"),
                                ("Post-NR", "Post-R", "Post")):
                if len(filled[a]) < 2 or len(filled[b]) < 2:
                    continue
                u, p = stats.mannwhitneyu(filled[a], filled[b],
                                          alternative="two-sided")
                test_rows.append(dict(
                    lineage=lineage, minor_cell_state=state, phase=phase,
                    n_NR=len(filled[a]), n_R=len(filled[b]),
                    mean_NR=float(filled[a].mean()), mean_R=float(filled[b].mean()),
                    rank_biserial_r=float(2.0 * u / (len(filled[a]) * len(filled[b])) - 1),
                    p_two_tailed=float(p)))

    # Specimens are named the way Supplementary Table 1 names them. The relabel
    # is in place, after every grouping and test, so no row moves and no value
    # changes.
    fractions = pd.concat(
        [t.assign(sample=to_study_ids(t["sample"], m)) for t, m in frac_rows],
        ignore_index=True)
    fractions.to_csv(OUT / "adaptive_state_fractions.csv", index=False)

    tests = pd.DataFrame(test_rows)
    tests["p_BH_within_lineage"] = np.nan
    for (lin, ph), grp in tests.groupby(["lineage", "phase"]):
        tests.loc[grp.index, "p_BH_within_lineage"] = multipletests(
            grp["p_two_tailed"], method="fdr_bh")[1]
    tests = tests.sort_values(["lineage", "phase", "p_two_tailed"])
    tests.to_csv(OUT / "adaptive_state_tests.csv", index=False)

    # ------------------------------------------------- Hallmark programmes
    hall = []
    for lineage, (_, gsea_name) in LINEAGES.items():
        for phase in ("pre", "post"):
            # The prepared GSEA/{pre,post} tables are not used: they descend
            # from MAST runs on a matrix that was never log1p CP10K
            # (00_Data_Audit/FINDINGS.md section 7). Repointed to module 12 on
            # the author's ruling of 2026-09-03; the Welch t-test branch is the
            # author's choice, consistent with the standing ruling that MAST is
            # dropped and with pretreatment_inflammatory.py, which draws the
            # other panels of this figure. Module 12 ranks non-responders
            # against responders, so no sign flip.
            f = RECOMPUTE_GSEA / f"{gsea_name}_{phase}_ttest_hallmark.csv"
            if not f.exists():
                continue
            d = pd.read_csv(f)
            d["NES_NRvsR"] = d["NES"]
            d = d.sort_values("NES_NRvsR", ascending=False)
            for _, r in d.head(5).iterrows():
                hall.append(dict(lineage=lineage, phase=phase, direction="NR-enriched",
                                 term=r["Term"], nes=float(r["NES_NRvsR"]),
                                 fdr_q=float(r["FDR q-val"])))
            for _, r in d.tail(5).iloc[::-1].iterrows():
                hall.append(dict(lineage=lineage, phase=phase, direction="R-enriched",
                                 term=r["Term"], nes=float(r["NES_NRvsR"]),
                                 fdr_q=float(r["FDR q-val"])))
    hallmark = pd.DataFrame(hall)
    hallmark.to_csv(OUT / "adaptive_hallmark_top.csv", index=False)

    # -------------------------------------------------------------- report
    L = ["ADAPTIVE IMMUNE RESOURCE - Reviewer 2 point R2.2", "=" * 96, ""]
    L.append("1. MINOR CELL STATE COMPOSITION")
    L.append("-" * 96)
    for lineage in LINEAGES:
        sub = tests[tests["lineage"] == lineage]
        L.append(f"  [{lineage}]  {sub['minor_cell_state'].nunique()} annotated states")
        for phase in ("Pre", "Post"):
            s = sub[sub["phase"] == phase].sort_values("p_two_tailed")
            if not len(s):
                continue
            sig = s[s["p_two_tailed"] < 0.05]
            L.append(f"     {phase}-treatment NR vs R: "
                     f"{len(sig)} of {len(s)} states with P < 0.05 "
                     f"(none survive BH within lineage)"
                     if not (s["p_BH_within_lineage"] < 0.05).any() else
                     f"     {phase}-treatment NR vs R: {len(sig)} of {len(s)} "
                     f"states with P < 0.05")
            for _, r in s.head(3).iterrows():
                L.append(f"        {r['minor_cell_state']:<32} "
                         f"NR {r['mean_NR']*100:5.1f}%  R {r['mean_R']*100:5.1f}%  "
                         f"P = {r['p_two_tailed']:.3f}  "
                         f"BH = {r['p_BH_within_lineage']:.3f}")
        L.append("")

    L.append("2. HALLMARK PROGRAMMES (positive NES = enriched in non-responders)")
    L.append("-" * 96)
    for lineage in LINEAGES:
        for phase in ("pre", "post"):
            sub = hallmark[(hallmark["lineage"] == lineage)
                           & (hallmark["phase"] == phase)]
            if not len(sub):
                continue
            L.append(f"  [{lineage}, {phase}-treatment]")
            for direction in ("NR-enriched", "R-enriched"):
                d = sub[sub["direction"] == direction].head(3)
                items = "; ".join(
                    f"{r['term']} (NES {r['nes']:+.2f}, q={r['fdr_q']:.2f})"
                    for _, r in d.iterrows())
                L.append(f"     top {direction}: {items}")
        L.append("")

    report = "\n".join(L)
    (OUT / "adaptive_immune_report.txt").write_text(report, encoding="utf-8")
    print(report)

    _panel_composition(fractions, tests)
    _panel_hallmark(hallmark)


def _panel_composition(fractions, tests):
    d = (S10 / "S10_D"); d.mkdir(parents=True, exist_ok=True)
    lineages = list(LINEAGES)
    fig, axes = plt.subplots(1, len(lineages),
                             figsize=(12.0 * SCALE * CM, 5.0 * SCALE * CM))
    for ax, lineage in zip(axes, lineages):
        sub = fractions[fractions["lineage"] == lineage]
        states = sorted(sub["minor_cell_state"].unique())
        means = (sub.groupby(["group4", "minor_cell_state"], observed=True)["fraction"]
                 .mean().unstack(fill_value=0).reindex(GROUPS).fillna(0))
        means = means.reindex(columns=states, fill_value=0)
        bottom = np.zeros(len(GROUPS))
        cmap = plt.get_cmap("tab20")
        for k, state in enumerate(states):
            ax.bar(range(len(GROUPS)), means[state], bottom=bottom, width=0.7,
                   color=cmap(k % 20), edgecolor="white", linewidth=0.4,
                   label=state.replace("_", " "))
            bottom += means[state].values
        ax.set_xticks(range(len(GROUPS)))
        ax.set_xticklabels(GROUPS, fontsize=5.5 * SCALE, rotation=45, ha="right")
        ax.set_title(lineage, fontsize=6.5 * SCALE)
        ax.set_ylim(0, 1)
        ax.tick_params(axis="both", labelsize=5 * SCALE, width=0.8, length=3)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.legend(loc="upper left", bbox_to_anchor=(0, -0.38), frameon=False,
                  fontsize=3.6 * SCALE, ncol=1, handlelength=1.0,
                  handletextpad=0.4, labelspacing=0.25)
    axes[0].set_ylabel("Fraction of lineage", fontsize=6 * SCALE)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.94, bottom=0.44, wspace=0.28)
    stem = d / "S10_D_adaptive_composition"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


def _panel_hallmark(hallmark):
    d = (S10 / "S10_E"); d.mkdir(parents=True, exist_ok=True)
    post = hallmark[hallmark["phase"] == "post"]
    lineages = list(LINEAGES)
    fig, axes = plt.subplots(1, len(lineages),
                             figsize=(15.0 * SCALE * CM, 4.4 * SCALE * CM))
    for ax, lineage in zip(axes, lineages):
        sub = post[post["lineage"] == lineage].drop_duplicates("term")
        sub = pd.concat([sub[sub["direction"] == "NR-enriched"].head(4),
                         sub[sub["direction"] == "R-enriched"].head(4)])
        sub = sub.sort_values("nes")
        y = np.arange(len(sub))
        colors = ["#B2182B" if v > 0 else "#2166AC" for v in sub["nes"]]
        ax.barh(y, sub["nes"], color=colors, edgecolor="#444444", linewidth=0.4,
                height=0.7)
        ax.axvline(0, color="#666666", linewidth=0.8)
        ax.set_yticks(y)
        # Hallmark names run to 34 characters; anything shorter than the
        # longest one here would cut "Interferon Gamma Response", which the
        # response letter quotes by name.
        ax.set_yticklabels([t if len(t) < 34 else t[:31] + "..."
                            for t in sub["term"]], fontsize=4.5 * SCALE)
        ax.set_title(lineage, fontsize=6.5 * SCALE)
        ax.tick_params(axis="x", labelsize=5 * SCALE, width=0.8, length=3)
        ax.tick_params(axis="y", length=0)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
    fig.supxlabel("NES, post-treatment (positive = enriched in non-responders)",
                  fontsize=6 * SCALE, y=0.03)
    fig.subplots_adjust(left=0.155, right=0.995, top=0.90, bottom=0.18,
                        wspace=1.55)
    stem = d / "S10_E_adaptive_hallmark"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
