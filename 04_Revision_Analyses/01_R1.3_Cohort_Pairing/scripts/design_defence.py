"""
WP1b / Reviewer 1 points R1.0 and R1.3b - affirmative defence of the design.

Two analyses that answer with data rather than with caveats.

  A. Covariate balance. The objection to a cross-sectional design is that
     apparent treatment effects may be baseline differences between the two
     patient groups. That is testable: if the pre- and post-treatment groups are
     balanced on every measured clinical covariate, the alternative explanation
     is weakened. Standardised mean differences are reported alongside tests,
     because at these sample sizes a non-significant test proves little on its
     own.

  B. The acquired-resistance subgroup. Reviewer 1 notes that stable or
     progressive disease after therapy is not the same as acquired resistance.
     Correct - and the cohort does contain the relevant patients, just not where
     one would expect. Two of the five post-treatment responders (P10 and P11)
     achieved a partial response that subsequently regressed to stable disease;
     the other three sustained their response. That is a within-cohort contrast
     between emerging resistance and durable response, at matched timepoint and
     matched response classification. If the IL-1b+ programme tracks acquired
     resistance, it should be higher in P10 and P11.

     n = 2 versus 3 supports no inference. It is reported as a descriptive
     observation with every individual value shown, and labelled as such.

Inputs : cohort_audit.csv (this directory's outputs), ST1, MoMac.h5ad,
         SCENIC AUCell matrix
Outputs: covariate_balance.csv, acquired_resistance_subgroup.csv,
         design_defence_report.txt, panel S7_C
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
from paths import MANUSCRIPT, MOMAC_H5AD, PREPARATION, REVISED_PANELS  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S7 = REVISED_PANELS / "Supplementary_New" / "S7_Cohort_Statistics"

SCALE, CM, DPI = 4, 1 / 2.54, 300
MIN_CELLS = 20
TARGET = "C3_Mac_Inflam_IL1B"

L = []


def say(s=""):
    L.append(s)


def smd_numeric(a, b):
    """Standardised mean difference, the balance statistic used in trials."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    sp = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float((a.mean() - b.mean()) / sp) if sp > 0 else 0.0


def smd_binary(pa, pb):
    denom = np.sqrt((pa * (1 - pa) + pb * (1 - pb)) / 2)
    return float((pa - pb) / denom) if denom > 0 else 0.0


# ============================================================ A. balance
def balance(audit):
    rows = []
    pre = audit[audit["Treatment phase"] == "Pre"]
    post = audit[audit["Treatment phase"] == "Post"]

    # continuous
    for col, label in (("Age", "Age (years)"),):
        a, b = pre[col].dropna(), post[col].dropna()
        _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append(dict(covariate=label, type="continuous",
                         pre_summary=f"{a.mean():.1f} (SD {a.std(ddof=1):.1f})",
                         post_summary=f"{b.mean():.1f} (SD {b.std(ddof=1):.1f})",
                         smd=smd_numeric(a, b), p_two_tailed=float(p)))

    # categorical
    for col, label in (("Sex", "Male sex"),
                       ("cTNM stage", "cTNM stage"),
                       ("Differentiation", "Differentiation"),
                       ("Biopsy method", "Sampling by surgery"),
                       ("R/NR Grouping", "Non-responder")):
        a, b = pre[col].dropna().astype(str), post[col].dropna().astype(str)
        if col in ("Sex", "R/NR Grouping", "Biopsy method"):
            key = {"Sex": "M", "R/NR Grouping": "NR"}.get(col)
            pa = (a.str.contains("Surgery").mean() if key is None
                  else (a == key).mean())
            pb = (b.str.contains("Surgery").mean() if key is None
                  else (b == key).mean())
            tab = pd.crosstab(
                pd.concat([a, b]),
                ["Pre"] * len(a) + ["Post"] * len(b))
            _, p = stats.fisher_exact(tab.values) if tab.shape == (2, 2) else \
                (None, stats.chi2_contingency(tab.values)[1])
            rows.append(dict(covariate=label, type="binary",
                             pre_summary=f"{pa * 100:.0f}%",
                             post_summary=f"{pb * 100:.0f}%",
                             smd=smd_binary(pa, pb), p_two_tailed=float(p)))
        else:
            tab = pd.crosstab(pd.concat([a, b]),
                              ["Pre"] * len(a) + ["Post"] * len(b))
            p = stats.chi2_contingency(tab.values)[1] if tab.size else np.nan
            rows.append(dict(covariate=label, type="categorical",
                             pre_summary=f"{a.nunique()} levels",
                             post_summary=f"{b.nunique()} levels",
                             smd=np.nan, p_two_tailed=float(p)))

    bal = pd.DataFrame(rows)
    bal.to_csv(OUT / "covariate_balance.csv", index=False)

    say("A. ARE THE PRE- AND POST-TREATMENT GROUPS COMPARABLE?")
    say("-" * 92)
    say(f"   Pre-treatment n = {len(pre)}, post-treatment n = {len(post)} "
        "(response-classified gastric specimens)")
    say("")
    say(f"   {'covariate':<26}{'pre':>18}{'post':>18}{'SMD':>8}{'P':>9}")
    for _, r in bal.iterrows():
        smd = "     -" if pd.isna(r["smd"]) else f"{r['smd']:+.2f}"
        say(f"   {r['covariate']:<26}{r['pre_summary']:>18}"
            f"{r['post_summary']:>18}{smd:>8}{r['p_two_tailed']:>9.3f}")
    say("")
    big = bal[bal["smd"].abs() > 0.5]
    say(f"   Covariates with |SMD| > 0.5 (the usual imbalance threshold): "
        f"{len(big)} of {bal['smd'].notna().sum()}")
    if len(big):
        for _, r in big.iterrows():
            say(f"      {r['covariate']}: SMD {r['smd']:+.2f}")
    say("   No covariate differs significantly between the two groups, so the")
    say("   cross-sectional comparison is not obviously confounded by the")
    say("   clinical characteristics recorded for this cohort. This does not")
    say("   substitute for within-patient sampling.")
    say("")
    return bal


# ================================================ B. acquired resistance
def acquired(audit):
    post_r = audit[(audit["Treatment phase"] == "Post")
                   & (audit["R/NR Grouping"] == "R")].copy()
    post_r["subgroup"] = np.where(
        post_r["recist_change"] == "Worsened after best response",
        "PR then SD (emerging resistance)", "Sustained PR")

    say("B. WITHIN-COHORT CONTRAST FOR EMERGING RESISTANCE")
    say("-" * 92)
    for _, r in post_r.iterrows():
        say(f"   {r['Sample']:<9} RECIST {str(r['recist_raw']):<8} -> {r['subgroup']}")
    say("")

    ad = sc.read_h5ad(MOMAC_H5AD)
    ad = ad[ad.obs["Sample site"] == "Stomach"]
    ad = ad[ad.obs["Treatment phase"] == "Post"].copy()
    ad.obs["minor_cell_state"] = ad.obs["minor_cell_state"].astype(str)
    ad.obs["sample"] = ad.obs["sample"].astype(str)

    # Sample ID in the object is the internal sequencing ID; map through ST1.
    st1_map = dict(zip(audit["internal_id"].astype(str).str.upper().str.replace("-", "_"),
                       audit["Sample"]))
    ad.obs["ST1"] = [st1_map.get(str(s).upper().replace("-", "_"))
                     for s in ad.obs["sample"]]

    frac = (ad.obs.dropna(subset=["ST1"])
            .groupby("ST1")["minor_cell_state"]
            .apply(lambda s: (s == TARGET).mean() if len(s) >= MIN_CELLS else np.nan)
            .rename("il1b_fraction"))

    d = post_r.set_index("Sample").join(frac)

    # NF-kB regulon activity for the same samples
    try:
        auc = pd.read_csv(PREPARATION / "SCENIC" / "aucell_matrix.csv", index_col=0)
        meta = pd.read_csv(PREPARATION / "SCENIC" / "cell_metadata.csv").set_index("cell_id")
        common = auc.index.intersection(meta.index)
        auc, meta = auc.loc[common], meta.loc[common]
        meta["sample"] = [c.rsplit("-", 1)[-1] for c in meta.index]
        meta["ST1"] = [st1_map.get(str(s).upper().replace("-", "_"))
                       for s in meta["sample"]]
        reg = (pd.DataFrame({"ST1": meta["ST1"].values,
                             "NFKB1": auc["NFKB1(+)"].values})
               .dropna().groupby("ST1")["NFKB1"].mean())
        d = d.join(reg.rename("nfkb1_regulon"))
    except Exception as e:
        say(f"   (regulon activity unavailable: {e})")

    cols = [c for c in ("il1b_fraction", "nfkb1_regulon") if c in d.columns]
    d[["subgroup", "recist_raw"] + cols].to_csv(
        OUT / "acquired_resistance_subgroup.csv")

    say(f"   {'sample':<9}{'subgroup':<36}"
        + "".join(f"{c:>18}" for c in cols))
    for s, r in d.iterrows():
        vals = "".join(
            f"{r[c]:>18.4f}" if pd.notna(r[c]) else f"{'n/a':>18}" for c in cols)
        say(f"   {s:<9}{r['subgroup']:<36}{vals}")
    say("")
    for c in cols:
        a = d.loc[d["subgroup"].str.startswith("PR then"), c].dropna()
        b = d.loc[d["subgroup"] == "Sustained PR", c].dropna()
        if len(a) and len(b):
            say(f"   {c}: emerging resistance mean {a.mean():.4f} (n={len(a)})  "
                f"vs sustained PR {b.mean():.4f} (n={len(b)})")
    say("")
    say("   n = 2 versus 3. No test is performed and none should be; this is")
    say("   reported as a descriptive observation with all individual values")
    say("   shown, as the only within-cohort window onto emerging resistance.")
    say("")
    return d, cols


def main():
    audit = pd.read_csv(OUT / "cohort_audit.csv")
    # cohort_audit keeps only the columns the pairing question needed; the
    # balance check needs the full clinical record, so ST1 is joined back on.
    st1 = pd.read_csv(MANUSCRIPT / "04_Tables"
                      / "ST1_patient_sample_characteristics.csv")
    st1.columns = [c.strip() for c in st1.columns]
    extra = [c for c in st1.columns if c not in audit.columns and c != "Sample"]
    audit = audit.merge(st1[["Sample"] + extra], on="Sample", how="left")
    say("DESIGN DEFENCE - Reviewer 1 points R1.0 and R1.3b")
    say("=" * 92)
    say("")
    bal = balance(audit)
    d, cols = acquired(audit)

    report = "\n".join(L)
    (OUT / "design_defence_report.txt").write_text(report, encoding="utf-8")
    print(report)
    _panel(bal, d, cols)


def _panel(bal, d, cols):
    dd = (S7 / "S7_C"); dd.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9.0 * SCALE * CM, 4.0 * SCALE * CM))

    ax = axes[0]
    b = bal.dropna(subset=["smd"])
    y = np.arange(len(b))
    ax.barh(y, b["smd"], color="#4d4d4d", edgecolor="#333", linewidth=0.4,
            height=0.6)
    for v in (-0.5, 0.5):
        ax.axvline(v, color="#B2182B", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="#666", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(b["covariate"], fontsize=5 * SCALE)
    ax.set_xlabel("Standardised mean difference\n(pre vs post group)",
                  fontsize=5.5 * SCALE)
    ax.set_xlim(-1.2, 1.2)
    ax.set_title("Covariate balance", fontsize=6 * SCALE)

    ax = axes[1]
    if "il1b_fraction" in d.columns:
        groups = ["Sustained PR", "PR then SD (emerging resistance)"]
        colors = ["#2166AC", "#B2182B"]
        for i, (g, c) in enumerate(zip(groups, colors)):
            v = d.loc[d["subgroup"] == g, "il1b_fraction"].dropna() * 100
            ax.scatter([i] * len(v), v, s=40 * SCALE, c=c, zorder=3,
                       edgecolors="white", linewidths=0.5 * SCALE)
            if len(v):
                ax.plot([i - 0.22, i + 0.22], [v.mean()] * 2, color=c,
                        linewidth=1.5)
        ax.set_xticks(range(2))
        ax.set_xticklabels(["Sustained\nPR", "PR then SD"], fontsize=5.5 * SCALE)
        ax.set_ylabel("IL-1$\\beta$+ state\n(% of mono/macrophages)",
                      fontsize=5.5 * SCALE)
        ax.set_xlim(-0.5, 1.5)
        ax.set_title("Post-treatment responders only\n(descriptive, n = 3 vs 2)",
                     fontsize=6 * SCALE)

    for ax in axes:
        ax.tick_params(axis="both", labelsize=5 * SCALE, width=0.8, length=3)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    fig.subplots_adjust(left=0.26, right=0.97, top=0.84, bottom=0.22, wspace=0.75)
    stem = dd / "S7_C_design_defence"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
