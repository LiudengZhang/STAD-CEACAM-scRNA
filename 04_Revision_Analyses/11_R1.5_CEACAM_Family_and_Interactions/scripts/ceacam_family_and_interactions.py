"""
Reviewer 1 point R1.5, follow-up requested by the corresponding author.

R1.5 asks us to justify treating CEACAM5 and CEACAM6 as one epithelial state.
The submitted answer showed that single-positive populations exist and that the
double-positive fraction carries the association. Three further questions are
answered here, all on data already in hand.

  A. Within the CEACAM family, which members are expressed at all in
     pre-treatment gastric tumour epithelium? If CEACAM5 and CEACAM6 dominate
     the family, "CEACAM5/6" names the expressed part of the family rather than
     an arbitrary pair.
  B. How tightly are the two co-expressed, per cell and per sample?
  C. Do they engage the same partners and the same receiver populations?
     Answered with CellPhoneDB v5 statistical analysis, with the pre-treatment
     epithelium split into CEACAM5-only, CEACAM6-only, double-positive and
     double-negative senders, so the two single-positive populations can be
     compared directly against each other.

Inputs : Round_5/01_Raw_Inputs/01_H5AD/*.h5ad  (pre-treatment gastric cells)
         Round_4 CellPhoneDB v5 database zip
Outputs: ceacam_family_expression.csv, ceacam_coexpression.csv,
         cpdb_ceacam_interactions.csv, interaction_target_overlap.csv,
         family_interactions_report.txt, R_family_and_interactions.png
"""

from pathlib import Path
import shutil
import sys
import zipfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import (CELL_TYPE_H5AD, CPDB_DB_ZIP, EPITHELIAL_H5AD)  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
# CellPhoneDB's inputs and raw output land beside outputs/ rather than in it,
# so that outputs/ holds results only and the working files do not travel
# with the submission archive.
WORK = Path(__file__).resolve().parents[1] / "_cpdb_work"

SEED = 0
# CellPhoneDB compares group means, so every group is capped at the same number
# of cells; the cap is set by the smallest population we want to keep separate.
CELLS_PER_GROUP = 600
MIN_CELLS_PER_GROUP = 50
ITERATIONS = 1000
THREADS = 8

STATES = {
    "CEACAM5+ only": "Epithelial_CEACAM5_only",
    "CEACAM6+ only": "Epithelial_CEACAM6_only",
    "Double positive": "Epithelial_CEACAM5_6_double",
    "Double negative": "Epithelial_CEACAM_negative",
}
SENDER_5 = "Epithelial_CEACAM5_only"
SENDER_6 = "Epithelial_CEACAM6_only"
SENDERS = [SENDER_5, SENDER_6, "Epithelial_CEACAM5_6_double"]

SCALE, CM, DPI = 4, 1 / 2.54, 300
TOP_FAMILY = 8

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def raw_adata(path):
    """The raw layer carries the full gene space and the normalised values."""
    ad = sc.read_h5ad(path)
    ad = ad[(ad.obs["Sample site"] == "Stomach")
            & (ad.obs["Treatment phase"] == "Pre")]
    src = ad.raw.to_adata() if ad.raw is not None else ad
    src.obs = ad.obs.copy()
    return src.copy()


def gene_vector(ad, gene):
    i = list(ad.var_names).index(gene)
    x = ad.X[:, i]
    return np.asarray(x.todense()).flatten() if hasattr(x, "todense") else \
        np.asarray(x).flatten()


def cpdb_genes():
    with zipfile.ZipFile(CPDB_DB_ZIP) as z:
        with z.open("gene_table.csv") as fh:
            g = pd.read_csv(fh)
    col = "hgnc_symbol" if "hgnc_symbol" in g.columns else "gene_name"
    return set(g[col].dropna().astype(str))


def main():
    rng = np.random.default_rng(SEED)
    L = ["CEACAM FAMILY EXPRESSION AND INTERACTION TARGETS - Reviewer 1 point R1.5",
         "=" * 94, ""]

    # =============================================== A. family expression
    epi = raw_adata(EPITHELIAL_H5AD)
    fam = sorted(g for g in epi.var_names if str(g).startswith("CEACAM"))
    if "CEACAM5" not in fam or "CEACAM6" not in fam:
        raise SystemExit("CEACAM5/CEACAM6 absent from the epithelial gene space")

    rows = []
    vals = {}
    for g in fam:
        v = gene_vector(epi, g)
        vals[g] = v
        rows.append(dict(gene=g, pct_cells_detected=100.0 * float((v > 0).mean()),
                         mean_expression=float(v.mean()),
                         mean_in_expressing=float(v[v > 0].mean()) if (v > 0).any()
                         else 0.0))
    family = pd.DataFrame(rows).sort_values("pct_cells_detected", ascending=False)
    family["rank"] = np.arange(1, len(family) + 1)
    family.to_csv(OUT / "ceacam_family_expression.csv", index=False)

    L += ["A. THE CEACAM FAMILY IN PRE-TREATMENT GASTRIC TUMOUR EPITHELIUM",
          "-" * 94,
          f"   Cells: {epi.n_obs:,}   family members in the gene space: {len(fam)}",
          "",
          "   gene            % cells detected   mean expression"]
    for _, r in family.iterrows():
        L.append(f"   {r['gene']:<14} {r['pct_cells_detected']:>10.2f}       "
                 f"{r['mean_expression']:>10.4f}")
    top2 = set(family.head(2)["gene"])
    L += ["",
          f"   The two most widely detected members are {', '.join(sorted(top2))}."
          if top2 == {"CEACAM5", "CEACAM6"} else
          f"   The two most widely detected members are {', '.join(sorted(top2))}, "
          "NOT CEACAM5 and CEACAM6.", ""]

    # ==================================================== B. co-expression
    c5, c6 = vals["CEACAM5"], vals["CEACAM6"]
    rho, p_rho = stats.spearmanr(c5, c6)
    samples = epi.obs["sample"].astype(str).values
    per_sample = pd.DataFrame({"sample": samples, "CEACAM5": c5, "CEACAM6": c6})
    sm = per_sample.groupby("sample").mean()
    rho_s, p_s = stats.spearmanr(sm["CEACAM5"], sm["CEACAM6"])

    pos5, pos6 = c5 > 0, c6 > 0
    either = pos5 | pos6
    coexp = pd.DataFrame([
        dict(level="single cell", n=len(c5), statistic="Spearman rho",
             value=float(rho), p=float(p_rho)),
        dict(level="sample mean", n=len(sm), statistic="Spearman rho",
             value=float(rho_s), p=float(p_s)),
        dict(level="single cell", n=int(either.sum()),
             statistic="% of CEACAM-positive cells that are double positive",
             value=100.0 * float((pos5 & pos6).sum() / either.sum()), p=np.nan),
        dict(level="single cell", n=int(pos5.sum()),
             statistic="% of CEACAM5+ cells that are also CEACAM6+",
             value=100.0 * float((pos5 & pos6).sum() / pos5.sum()), p=np.nan),
        dict(level="single cell", n=int(pos6.sum()),
             statistic="% of CEACAM6+ cells that are also CEACAM5+",
             value=100.0 * float((pos5 & pos6).sum() / pos6.sum()), p=np.nan),
    ])
    coexp.to_csv(OUT / "ceacam_coexpression.csv", index=False)

    L += ["B. CO-EXPRESSION OF CEACAM5 AND CEACAM6", "-" * 94]
    for _, r in coexp.iterrows():
        tail = "" if np.isnan(r["p"]) else f"   P = {r['p']:.3g}"
        L.append(f"   {r['statistic']:<52} {r['value']:>8.3f}"
                 f"   (n = {r['n']:,}){tail}")
    L.append("")

    # ============================== C. interaction targets, CellPhoneDB v5
    state = np.select([pos5 & ~pos6, ~pos5 & pos6, pos5 & pos6],
                      ["CEACAM5+ only", "CEACAM6+ only", "Double positive"],
                      default="Double negative")
    epi.obs["cpdb_label"] = pd.Series(state, index=epi.obs_names).map(STATES).values

    keep_genes = cpdb_genes()
    frames = [epi[:, [g for g in epi.var_names if g in keep_genes]].copy()]
    del epi

    for name, path in CELL_TYPE_H5AD.items():
        if name == "Epithelial":
            continue
        ad = raw_adata(path)
        ad.obs["cpdb_label"] = ad.obs["major_cell_type"].astype(str).values
        frames.append(ad[:, [g for g in ad.var_names if g in keep_genes]].copy())
        del ad

    shared = set(frames[0].var_names)
    for f in frames[1:]:
        shared &= set(f.var_names)
    shared = sorted(shared)
    frames = [f[:, shared].copy() for f in frames]

    merged = sc.concat(frames, join="outer", index_unique=None)
    del frames

    # Equal-sized groups, seeded, so no cell type dominates the null.
    idx = []
    counts = merged.obs["cpdb_label"].value_counts()
    dropped = []
    for label, n in counts.items():
        if n < MIN_CELLS_PER_GROUP:
            dropped.append((label, int(n)))
            continue
        pos = np.flatnonzero((merged.obs["cpdb_label"] == label).values)
        if n > CELLS_PER_GROUP:
            pos = rng.choice(pos, CELLS_PER_GROUP, replace=False)
        idx.append(pos)
    sub = merged[np.sort(np.concatenate(idx))].copy()
    del merged

    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    counts_h5ad = WORK / "counts.h5ad"
    meta_tsv = WORK / "meta.tsv"
    sub.obs = sub.obs[["cpdb_label"]].copy()
    sub.write_h5ad(counts_h5ad)
    pd.DataFrame({"Cell": sub.obs_names,
                  "cell_type": sub.obs["cpdb_label"].astype(str).values}
                 ).to_csv(meta_tsv, sep="\t", index=False)

    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

    cpdb_statistical_analysis_method.call(
        cpdb_file_path=str(CPDB_DB_ZIP),
        meta_file_path=str(meta_tsv),
        counts_file_path=str(counts_h5ad),
        counts_data="hgnc_symbol",
        output_path=str(WORK),
        iterations=ITERATIONS,
        threshold=0.1,
        threads=THREADS,
        debug_seed=SEED,
        pvalue=0.05,
        output_suffix="ceacam",
    )

    means = pd.read_csv(next(WORK.glob("statistical_analysis_means_*.txt")), sep="\t")
    pvals = pd.read_csv(next(WORK.glob("statistical_analysis_pvalues_*.txt")), sep="\t")

    meta_cols = [c for c in means.columns if "|" not in c]
    pair_cols = [c for c in means.columns if "|" in c]

    is_ceacam = means["gene_a"].isin(["CEACAM5", "CEACAM6"]) | \
        means["gene_b"].isin(["CEACAM5", "CEACAM6"])
    m = means[is_ceacam].set_index("id_cp_interaction")
    pv = pvals[pvals["id_cp_interaction"].isin(m.index)].set_index("id_cp_interaction")

    long = []
    for iid in m.index:
        for col in pair_cols:
            sender, receiver = col.split("|", 1)
            long.append(dict(
                interaction=m.loc[iid, "interacting_pair"],
                gene_a=m.loc[iid, "gene_a"], gene_b=m.loc[iid, "gene_b"],
                classification=m.loc[iid, "classification"],
                sender=sender, receiver=receiver,
                mean=float(m.loc[iid, col]), p=float(pv.loc[iid, col])))
    inter = pd.DataFrame(long)
    inter["significant"] = inter["p"] < 0.05
    inter.to_csv(OUT / "cpdb_ceacam_interactions.csv", index=False)

    # Which partner does each marker reach, and in which receiver population?
    # gene_a sits in the sender, gene_b in the receiver, so the CEACAM5-only and
    # CEACAM6-only senders are compared on the partner-receiver pairs they hit.
    def targets(sender, marker):
        sel = inter[(inter["sender"] == sender) & inter["significant"]
                    & (inter["gene_a"] == marker)]
        return set(zip(sel["gene_b"], sel["receiver"]))

    t5 = targets(SENDER_5, "CEACAM5")
    t6 = targets(SENDER_6, "CEACAM6")
    shared_t = sorted(t5 & t6)
    overlap = pd.DataFrame([
        dict(sender=SENDER_5, marker="CEACAM5", n_targets=len(t5)),
        dict(sender=SENDER_6, marker="CEACAM6", n_targets=len(t6)),
        dict(sender="shared", marker="CEACAM5 and CEACAM6",
             n_targets=len(shared_t)),
    ])
    overlap["jaccard"] = (len(t5 & t6) / len(t5 | t6)) if (t5 | t6) else np.nan
    overlap.to_csv(OUT / "interaction_target_overlap.csv", index=False)

    L += ["C. INTERACTION TARGETS, CellPhoneDB v5 STATISTICAL ANALYSIS", "-" * 94,
          f"   Cells analysed: {sub.n_obs:,} in {sub.obs['cpdb_label'].nunique()} "
          f"groups, capped at {CELLS_PER_GROUP} cells per group",
          f"   {ITERATIONS} permutations, seed {SEED}"]
    if dropped:
        L.append("   groups below the minimum and therefore dropped: "
                 + ", ".join(f"{d} (n = {n})" for d, n in dropped))
    L += ["",
          "   Annotated partners of each marker in the database:"]
    for marker in ("CEACAM5", "CEACAM6"):
        partners = sorted(set(
            inter.loc[inter["gene_a"] == marker, "gene_b"]) | set(
            inter.loc[inter["gene_b"] == marker, "gene_a"]))
        L.append(f"      {marker:<9} {', '.join(partners)}")
    L += ["",
          f"   Significant partner-receiver pairs (P < 0.05) engaged by the",
          f"   CEACAM5-only sender: {len(t5)}; by the CEACAM6-only sender: "
          f"{len(t6)}; shared: {len(shared_t)}",
          f"   Jaccard overlap: {overlap['jaccard'].iloc[0]:.3f}", ""]
    for partner, receiver in shared_t:
        L.append(f"      {partner:<10} on {receiver}")
    L.append("")

    (OUT / "family_interactions_report.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    _panel(family, inter)


def _panel(family, inter):
    """One figure for the response letter: what is expressed, and what it talks to."""
    pretty = {"Epithelial_CEACAM5_only": "CEACAM5-only epi.",
              "Epithelial_CEACAM6_only": "CEACAM6-only epi.",
              "Epithelial_CEACAM5_6_double": "double-positive epi.",
              "Epithelial_CEACAM_negative": "CEACAM-negative epi.",
              "T & NK cell": "T/NK cells", "Pericyte": "pericytes"}
    short = {SENDER_5: "CEACAM5\nonly", SENDER_6: "CEACAM6\nonly",
             "Epithelial_CEACAM5_6_double": "CEACAM5/6\ndouble"}

    fig, axes = plt.subplots(
        1, 2, figsize=(19.0 * SCALE * CM, 6.5 * SCALE * CM),
        gridspec_kw=dict(width_ratios=[0.75, 1.0]))

    ax = axes[0]
    top = family.head(TOP_FAMILY).iloc[::-1]
    colors = ["#B2182B" if g in ("CEACAM5", "CEACAM6") else "#9e9e9e"
              for g in top["gene"]]
    ax.barh(np.arange(len(top)), top["pct_cells_detected"], color=colors,
            edgecolor="#333333", linewidth=0.4, height=0.7)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["gene"], fontsize=5.5 * SCALE)
    ax.set_xlabel("% of pre-treatment gastric epithelial cells\n"
                  "with the transcript detected", fontsize=6 * SCALE)
    ax.set_title("CEACAM family, ranked", fontsize=6.5 * SCALE)

    # One row per annotated interaction, not per partner: CEACAM5 and CEACAM6
    # can both reach the same partner, and collapsing them would hide which
    # marker carries the pair.
    ax = axes[1]
    sig = inter[inter["sender"].isin(SENDERS) & inter["significant"]].copy()
    sig["row"] = (sig["gene_a"] + " \u2192 " + sig["gene_b"] + " on "
                  + sig["receiver"].map(lambda r: pretty.get(r, r)))
    rows = sorted(sig["row"].unique())
    for j, sender in enumerate(SENDERS):
        sub = sig[sig["sender"] == sender]
        for _, r in sub.iterrows():
            ax.scatter(j, rows.index(r["row"]), s=r["mean"] * 22 * SCALE,
                       c="#B2182B" if r["gene_a"] == "CEACAM5" else "#2166AC",
                       edgecolors="white", linewidths=0.4 * SCALE, zorder=3)
    ax.set_xticks(range(len(SENDERS)))
    ax.set_xticklabels([short[s] for s in SENDERS], fontsize=5.5 * SCALE)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=5 * SCALE)
    ax.set_xlim(-0.6, len(SENDERS) - 0.4)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlabel("Sending epithelial population", fontsize=6 * SCALE)
    ax.set_title("Significant ligand-receptor pairs, CellPhoneDB v5 (P < 0.05)\n"
                 "red = through CEACAM5, blue = through CEACAM6",
                 fontsize=5.5 * SCALE)
    ax.grid(axis="y", color="#e6e6e6", linewidth=0.6, zorder=0)

    for a in axes:
        a.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.84, bottom=0.20, wspace=1.4)
    out = OUT / "R_family_and_interactions.png"
    fig.savefig(out, dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    if "--plot-only" in sys.argv:
        _panel(pd.read_csv(OUT / "ceacam_family_expression.csv"),
               pd.read_csv(OUT / "cpdb_ceacam_interactions.csv"))
    else:
        main()
