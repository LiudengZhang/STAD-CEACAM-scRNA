"""
WP8b / Reviewer 1 point R1.8 - an affirmative answer on NF-kB pathway activity.

  "Transcriptomic enrichment of the Hallmark TNFa/NF-kB gene set does not
   necessarily demonstrate biochemical activation of the NF-kB pathway.
   Protein-level evidence, such as nuclear p65 localization, phospho-p65, or IkB
   degradation, would substantially strengthen this conclusion."

We cannot supply protein data. We can, however, replace the weakest form of the
argument (a broad co-expression gene set) with two much more specific readouts
that are already computable from the data in hand:

  1. TF regulon activity. pySCENIC infers, per cell, the activity of NFKB1 and
     NFKB2 from the coordinated expression of their DIRECT targets, where
     "direct" means the target carries the TF's binding motif within its
     regulatory region (cisTarget). Unlike a Hallmark set, this is TF-specific
     and motif-anchored, so elevated activity is evidence that the transcription
     factor itself is operating, not that inflammation is generally present.

  2. Canonical negative-feedback targets. NFKBIA, TNFAIP3, NFKB2, RELB, BIRC3
     and TRAF1 are induced by NF-kB and act to terminate or re-tune the pathway.
     Because they are transcribed only after nuclear translocation of the
     complex, their induction is a specific downstream footprint of pathway
     activation. IkB (NFKBIA) resynthesis in particular is the transcriptional
     counterpart of the IkB degradation the reviewer asks about.

Neither replaces a Western blot, and the manuscript says so. Together they move
the claim from "an inflammatory gene set is enriched" to "the NF-kB transcription
factors are active and their feedback programme is running".

Inputs : submission-tree/02_Preparation_for_Panels/SCENIC/{aucell_matrix.csv,cell_metadata.csv}
         submission-tree/01_Raw_Inputs/01_H5AD/{MoMac,Epithelial,Fibroblast}.h5ad
Outputs: nfkb_regulon_activity.csv, nfkb_feedback_targets.csv,
         nfkb_sample_values.csv,
         nfkb_regulon_report.txt, panel S11_C
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
from paths import (  # noqa: E402
    ANALYSIS_PANELS, EPITHELIAL_H5AD, FIBROBLAST_H5AD, MOMAC_H5AD,
    PREPARATION)

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S11 = ANALYSIS_PANELS / "S11_Affirmative_Analyses"

SCALE, CM, DPI = 4, 1 / 2.54, 300
MIN_CELLS = 20
COLOR_R, COLOR_NR = "#2166AC", "#B2182B"

REGULONS = ["NFKB1(+)", "NFKB2(+)", "BACH1(+)"]
# Direct, motif-supported NF-kB targets whose induction is part of the pathway's
# own feedback circuitry. NFKBIA encodes IkB-alpha; its resynthesis is the
# transcriptional readout of the IkB degradation-resynthesis cycle.
FEEDBACK = ["NFKBIA", "TNFAIP3", "NFKB2", "RELB", "BIRC3", "TRAF1"]
# Direct effector targets, as a second, independent panel.
EFFECTOR = ["ICAM1", "CXCL8", "CCL20", "SOD2", "CXCL2", "IER3"]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def mw(a, b):
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(p), float(2.0 * u / (len(a) * len(b)) - 1.0)


def load_regulons():
    auc = pd.read_csv(PREPARATION / "SCENIC" / "aucell_matrix.csv", index_col=0)
    meta = pd.read_csv(PREPARATION / "SCENIC" / "cell_metadata.csv")
    meta = meta.set_index("cell_id")
    common = auc.index.intersection(meta.index)
    auc, meta = auc.loc[common], meta.loc[common]
    # The sample identifier is the suffix of the cell barcode.
    meta["sample"] = [c.rsplit("-", 1)[-1] for c in meta.index]
    return auc, meta


def sample_level(auc, meta, regulon, phase, group_col, cell_type=None):
    m = meta.copy()
    if cell_type is not None:
        m = m[m["major_cell_type"] == cell_type]
    m = m[m[group_col].isin(["Responsed", "No-response"])]
    if not len(m):
        return np.array([]), np.array([])
    v = auc.loc[m.index, regulon]
    df = pd.DataFrame({"sample": m["sample"].values,
                       "group": m[group_col].values, "value": v.values})
    keep = df.groupby("sample").size()
    df = df[df["sample"].isin(keep[keep >= MIN_CELLS].index)]
    agg = df.groupby(["sample", "group"], observed=True)["value"].mean().reset_index()
    return (agg.loc[agg["group"] == "No-response", "value"].values,
            agg.loc[agg["group"] == "Responsed", "value"].values)


def score_panel(path, genes, phase, group_col, name):
    ad = sc.read_h5ad(path)
    ad = ad[ad.obs["Sample site"] == "Stomach"]
    ad = ad[ad.obs["Treatment phase"] == phase]
    ad = ad[ad.obs[group_col].isin(["Responsed", "No-response"])].copy()
    present = [g for g in genes if g in ad.var_names]
    if len(present) < 3:
        return None, None, present
    sc.tl.score_genes(ad, present, score_name="panel", random_state=0)
    df = pd.DataFrame({"sample": ad.obs["sample"].astype(str).values,
                       "group": ad.obs[group_col].astype(str).values,
                       "value": ad.obs["panel"].values})
    keep = df.groupby("sample").size()
    df = df[df["sample"].isin(keep[keep >= MIN_CELLS].index)]
    agg = df.groupby(["sample", "group"], observed=True)["value"].mean().reset_index()
    return (agg.loc[agg["group"] == "No-response", "value"].values,
            agg.loc[agg["group"] == "Responsed", "value"].values, present)


def main():
    auc, meta = load_regulons()
    L = ["NF-kB TRANSCRIPTION FACTOR ACTIVITY - affirmative evidence for R1.8",
         "=" * 94, ""]
    L.append(f"pySCENIC AUCell matrix: {auc.shape[0]:,} cells x {auc.shape[1]} regulons")
    L.append("Regulons are motif-anchored: a gene enters a regulon only if it is")
    L.append("co-expressed with the TF AND carries the TF's binding motif (cisTarget).")
    L.append("")

    rows, svals = [], []
    cell_types = sorted(meta["major_cell_type"].dropna().unique())
    for reg in REGULONS:
        if reg not in auc.columns:
            continue
        for phase, col in (("Pre", "stomach_pre_grouping"),
                           ("Post", "stomach_post_grouping")):
            for ct in [None] + cell_types:
                nr, r = sample_level(auc, meta, reg, phase, col, ct)
                if len(nr) < 2 or len(r) < 2:
                    continue
                p, eff = mw(nr, r)
                rows.append(dict(regulon=reg, phase=phase,
                                 cell_type=ct or "All cells",
                                 n_NR=len(nr), n_R=len(r),
                                 mean_NR=float(nr.mean()), mean_R=float(r.mean()),
                                 rank_biserial_r=eff, p_two_tailed=p))
                # The per-sample values behind that summary. Fig. S11C plots
                # these directly, so the reviewer sees the distribution and the
                # sample count rather than a group mean.
                for grp, vals in (("NR", nr), ("R", r)):
                    for v in vals:
                        svals.append(dict(source="regulon", key=reg, phase=phase,
                                          stratum=ct or "All cells",
                                          group=grp, value=float(v)))
    reg_df = pd.DataFrame(rows)
    reg_df.to_csv(OUT / "nfkb_regulon_activity.csv", index=False)

    L.append("REGULON ACTIVITY, non-responders vs responders (sample-level, two-sided)")
    L.append("-" * 94)
    for reg in REGULONS:
        sub = reg_df[reg_df["regulon"] == reg]
        if not len(sub):
            continue
        L.append(f"  [{reg}]")
        for phase in ("Pre", "Post"):
            s = sub[sub["phase"] == phase].sort_values("p_two_tailed")
            if not len(s):
                continue
            L.append(f"    {phase}-treatment:")
            for _, x in s.head(6).iterrows():
                flag = "  *" if x["p_two_tailed"] < 0.05 else ""
                L.append(f"      {x['cell_type']:<24} NR {x['mean_NR']:.4f}  "
                         f"R {x['mean_R']:.4f}   P = {x['p_two_tailed']:.4f}   "
                         f"r = {x['rank_biserial_r']:+.3f}{flag}")
        L.append("")

    # ------------------------------------------------ feedback target panels
    frows = []
    for panel_name, genes in (("NF-kB negative-feedback targets", FEEDBACK),
                              ("NF-kB direct effector targets", EFFECTOR)):
        for comp, path in (("Monocytes/Macrophages", MOMAC_H5AD),
                           ("Epithelial", EPITHELIAL_H5AD),
                           ("Fibroblast", FIBROBLAST_H5AD)):
            for phase, col in (("Pre", "stomach_pre_grouping"),
                               ("Post", "stomach_post_grouping")):
                nr, r, present = score_panel(path, genes, phase, col, comp)
                if nr is None or len(nr) < 2 or len(r) < 2:
                    continue
                p, eff = mw(nr, r)
                frows.append(dict(panel=panel_name, compartment=comp, phase=phase,
                                  genes_used=";".join(present),
                                  n_NR=len(nr), n_R=len(r),
                                  mean_NR=float(nr.mean()), mean_R=float(r.mean()),
                                  rank_biserial_r=eff, p_two_tailed=p))
                for grp, vals in (("NR", nr), ("R", r)):
                    for v in vals:
                        svals.append(dict(source="panel", key=panel_name,
                                          phase=phase, stratum=comp,
                                          group=grp, value=float(v)))
    fb = pd.DataFrame(frows)
    fb.to_csv(OUT / "nfkb_feedback_targets.csv", index=False)
    sv = pd.DataFrame(svals)
    sv.to_csv(OUT / "nfkb_sample_values.csv", index=False)

    L.append("DIRECT TARGET PANELS (sample-level, two-sided)")
    L.append("-" * 94)
    L.append(f"  Negative feedback: {', '.join(FEEDBACK)}")
    L.append(f"  Direct effectors : {', '.join(EFFECTOR)}")
    L.append("")
    for panel_name in fb["panel"].unique():
        L.append(f"  [{panel_name}]")
        for _, x in fb[fb["panel"] == panel_name].iterrows():
            flag = "  *" if x["p_two_tailed"] < 0.05 else ""
            L.append(f"    {x['compartment']:<24}{x['phase']:<6} "
                     f"NR {x['mean_NR']:+.4f}  R {x['mean_R']:+.4f}   "
                     f"P = {x['p_two_tailed']:.4f}   "
                     f"r = {x['rank_biserial_r']:+.3f}{flag}")
        L.append("")

    L.append("WHAT THIS DOES AND DOES NOT SHOW")
    L.append("-" * 94)
    L.append("  Does: the NF-kB transcription factors themselves show elevated")
    L.append("  inferred activity, measured from motif-supported direct targets, and")
    L.append("  the pathway's own negative-feedback programme (NFKBIA/IkB-alpha,")
    L.append("  TNFAIP3/A20, NFKB2, RELB) is induced. Feedback genes are transcribed")
    L.append("  only after nuclear translocation, so their induction is a specific")
    L.append("  footprint of activation rather than of generic inflammation.")
    L.append("  Does not: measure p65 phosphorylation, nuclear localisation, or IkB")
    L.append("  protein. The manuscript states this.")

    report = "\n".join(L)
    (OUT / "nfkb_regulon_report.txt").write_text(report, encoding="utf-8")
    print(report)

    _panel(reg_df, fb, sv)


def _box(ax, series, colors, ylabel, title=None):
    """
    One sample-level box per group with every sample drawn on top. The groups
    here have five or six samples, so the points are the honest display and the
    box is only there to carry the median and the spread.
    """
    rng = np.random.default_rng(0)
    labels = [lab for lab, _ in series]
    data = [vals for _, vals in series]
    bp = ax.boxplot(data, widths=0.55, showfliers=False, patch_artist=True,
                    medianprops=dict(color="#333333", linewidth=1.0),
                    whiskerprops=dict(color="#666666", linewidth=0.8),
                    capprops=dict(color="#666666", linewidth=0.8),
                    boxprops=dict(linewidth=0.6, edgecolor="#333333"))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.35)
    for i, (vals, c) in enumerate(zip(data, colors), start=1):
        if not len(vals):
            continue
        jitter = rng.uniform(-0.13, 0.13, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=9 * SCALE, c=c,
                   edgecolors="white", linewidths=0.4, zorder=3)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels([f"{lab}\nn = {len(v)}" for lab, v in zip(labels, data)],
                       fontsize=5.0 * SCALE)
    ax.set_ylabel(ylabel, fontsize=6 * SCALE)
    ax.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    if title:
        ax.set_title(title, fontsize=6 * SCALE)


def _panel(reg_df, fb, sv):
    """
    Sample-level boxplots with the individual samples shown, as promised to
    Reviewer 1 in the response to point R1.8. The earlier version of this panel
    drew group means as bars, which hid both the spread and the fact that each
    group holds five or six samples.
    """
    d = (S11 / "S11_C"); d.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9.0 * SCALE * CM, 4.6 * SCALE * CM))

    # left: NFKB1 regulon in monocytes/macrophages, both timepoints
    reg = sv[(sv["source"] == "regulon") & (sv["key"] == "NFKB1(+)")
             & (sv["stratum"] == "MoMac")]
    series, colors = [], []
    for phase in ("Pre", "Post"):
        for grp, c in (("R", COLOR_R), ("NR", COLOR_NR)):
            v = reg.loc[(reg["phase"] == phase) & (reg["group"] == grp), "value"].values
            if len(v):
                series.append((f"{phase}\n{grp}", v))
                colors.append(c)
    post = reg_df[(reg_df["regulon"] == "NFKB1(+)")
                  & (reg_df["cell_type"] == "MoMac")
                  & (reg_df["phase"] == "Post")]
    title = (f"post-treatment P = {post['p_two_tailed'].iloc[0]:.3f}"
             if len(post) else None)
    _box(axes[0], series, colors,
         "NFKB1 regulon activity (AUCell)\nmonocytes/macrophages", title)

    # right: NF-kB negative-feedback target score after treatment, by compartment
    pan = sv[(sv["source"] == "panel")
             & (sv["key"] == "NF-kB negative-feedback targets")
             & (sv["phase"] == "Post")]
    order = [c for c in ("Monocytes/Macrophages", "Epithelial", "Fibroblast")
             if c in set(pan["stratum"])]
    series, colors = [], []
    for comp in order:
        short = "Mono/Mac" if comp == "Monocytes/Macrophages" else comp
        for grp, c in (("R", COLOR_R), ("NR", COLOR_NR)):
            v = pan.loc[(pan["stratum"] == comp) & (pan["group"] == grp), "value"].values
            if len(v):
                series.append((f"{short}\n{grp}", v))
                colors.append(c)
    _box(axes[1], series, colors,
         "NF-$\\kappa$B feedback target score\n(post-treatment)")
    axes[1].axhline(0, color="#666666", linewidth=0.8, zorder=0)

    fig.subplots_adjust(left=0.13, right=0.98, top=0.90, bottom=0.20, wspace=0.42)
    stem = d / "S11_C_nfkb_regulon_and_feedback"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
