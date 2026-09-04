"""
Reviewer 1 point R1.5, the premise rather than the request.

The reviewer's objection opens with a reason: the manuscript treats CEACAM5 and
CEACAM6 as one entity "because their expression is highly correlated". The
number that claim rested on was a per-cell Spearman rho of 0.93. It is not a
per-cell correlation - it is the correlation between sample means, and the
per-cell value is 0.44. Reporting the sample-level figure as if it were
cell-level is what this module corrects, and it does so by answering the
question the reviewer is really asking: are these two genes expressed by the
same cells, or do they merely rise and fall together across patients?

Dropout makes that question harder than it looks. A cell expressing both
transcripts is routinely sequenced as expressing one, or neither, so the
per-cell correlation is attenuated and single-positive cells are manufactured.
Three analyses separate the artefact from the biology, each with its own
control:

  A. Sequencing depth. Cells are split into quintiles of total UMI. If dropout
     drives the picture, then with depth the correlation rises, double-positive
     cells become commoner and double-negative cells rarer - and the apparent
     single-positive population shrinks.
  B. Co-detection against independence, within each depth stratum. Detection
     rates rise with depth, so comparing raw double-positive fractions across
     strata proves nothing. What matters is whether the two genes are detected
     together more often than their own detection rates predict. This test is
     computed inside a stratum and is therefore immune to both dropout and to
     any averaging.
  C. Between-patient covariation. A per-cell correlation can be manufactured
     entirely by patient-level differences. Two controls bound it: the same
     correlation computed inside each sample, and a null in which CEACAM6 is
     permuted against CEACAM5 within sample, which destroys cell-level pairing
     while preserving every between-patient difference.

A metacell sweep is included with a control of its own, because metacells are
the obvious response to dropout and the obvious thing a reader will ask for.
Cells are pooled by kNN within sample at a range of sizes, against pooling the
same number of cells at random. The two curves are reported together: past a
small k they coincide, which says that the rise of rho with metacell size is
arithmetic - averaging shrinks the zeros and pulls every pool toward the sample
mean - rather than recovered biology. The metacell view is therefore a
presentation of the data, not a dropout-corrected measurement of it, and the
co-expression claim rests on B.

Everything is on pre-treatment stomach epithelial cells, the population the
manuscript's CEACAM statements are about, with no response filter: this asks
whether the genes mark one state, not whether that state predicts response.

Inputs : Round_5/01_Raw_Inputs/01_H5AD/Epithelial.h5ad
Outputs: dropout_depth_strata.csv, coexpression_within_sample.csv,
         metacell_sweep.csv, dropout_coexpression_report.txt
         panel S8_H
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD, REVISED_PANELS  # noqa: E402
from shared.sample_ids import sample_id_map, to_study_ids  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
S8 = REVISED_PANELS / "Supplementary_New" / "S8_CEACAM_Metaprogram"

SCALE, CM, DPI = 4, 1 / 2.54, 300
SEED = 42
N_STRATA = 5
KS = [1, 5, 10, 20, 50, 100, 200]
N_PCS = 30
COLOR_5, COLOR_6 = "#2166AC", "#B2182B"
COLOR_REAL, COLOR_SHUF = "#333333", "#B0B0B0"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def load():
    """Pre-treatment stomach epithelium, expression from .raw."""
    ad = sc.read_h5ad(EPITHELIAL_H5AD)
    # Resolved before subsetting, so the crosswalk is checked against every
    # specimen in the object rather than only the pre-treatment stomach ones.
    ids = sample_id_map(ad.obs)
    ad = ad[(ad.obs["Sample site"] == "Stomach")
            & (ad.obs["Treatment phase"] == "Pre")].copy()
    if ad.raw is None:
        raise SystemExit("Epithelial.h5ad has no .raw - refusing to use .X")
    src = ad.raw
    out = {}
    for g in ("CEACAM5", "CEACAM6"):
        i = list(src.var_names).index(g)
        x = src.X[:, i]
        out[g] = x.toarray().flatten() if hasattr(x, "toarray") \
            else np.asarray(x).flatten()
    df = pd.DataFrame(out)
    df["sample"] = ad.obs["sample"].astype(str).values
    df["total_counts"] = ad.obs["total_counts"].astype(float).values
    if "X_pca" not in ad.obsm:
        raise SystemExit("Epithelial.h5ad has no X_pca - cannot build metacells")
    return df, np.asarray(ad.obsm["X_pca"])[:, :N_PCS], ids


def depth_strata(df):
    """A. and B.: dropout signature, and co-detection against independence."""
    q = pd.qcut(df["total_counts"], N_STRATA, labels=False)
    rows = []
    for k in range(N_STRATA):
        d = df[q == k]
        a, b = d["CEACAM5"] > 0, d["CEACAM6"] > 0
        rho, p_rho = stats.spearmanr(d["CEACAM5"], d["CEACAM6"])
        obs, exp = float((a & b).mean()), float(a.mean() * b.mean())
        table = np.array([[int((a & b).sum()), int((a & ~b).sum())],
                          [int((~a & b).sum()), int((~a & ~b).sum())]])
        odds, p_fisher = stats.fisher_exact(table)
        any_one = int((a | b).sum())
        rows.append(dict(
            stratum=f"Q{k + 1}", n_cells=len(d),
            median_total_counts=float(d["total_counts"].median()),
            rho=float(rho), p_rho=float(p_rho),
            pct_ceacam5_pos=100 * float(a.mean()),
            pct_ceacam6_pos=100 * float(b.mean()),
            pct_double_pos=100 * obs,
            pct_single_pos=100 * float((a ^ b).mean()),
            pct_double_neg=100 * float((~a & ~b).mean()),
            pct_double_of_expressing=100 * float((a & b).sum()) / max(any_one, 1),
            expected_double_pos=100 * exp,
            obs_over_expected=obs / exp if exp else np.nan,
            odds_ratio=float(odds), p_fisher=float(p_fisher)))
    return pd.DataFrame(rows)


def within_sample(df, rng, ids):
    """
    C.: the correlation with every between-patient difference removed.

    Grouping stays on the internal label and the study sample ID is written
    into the row instead. Relabelling before the loop would reorder the groups
    and so change which permutation each sample draws from `rng`, moving the
    null column for no reason; the identifiers written out are the same either
    way.
    """
    rows = []
    for s, d in df.groupby("sample", sort=True):
        if len(d) < 20:
            continue
        rho, p = stats.spearmanr(d["CEACAM5"], d["CEACAM6"])
        # Permuting CEACAM6 against CEACAM5 inside the sample destroys the
        # cell-level pairing and keeps everything else.
        perm = rng.permutation(len(d))
        rho_null, _ = stats.spearmanr(d["CEACAM5"].values,
                                      d["CEACAM6"].values[perm])
        rows.append(dict(sample=to_study_ids([s], ids)[0], n_cells=len(d),
                         rho=float(rho), p=float(p),
                         rho_permuted=float(rho_null)))
    return pd.DataFrame(rows).sort_values("rho").reset_index(drop=True)


def _pools_knn(emb, k):
    """A disjoint cover of the sample by kNN pools of exactly k cells."""
    nn = NearestNeighbors(n_neighbors=k).fit(emb)
    _, nbr = nn.kneighbors(emb)
    used, pools = np.zeros(len(emb), bool), []
    for row in nbr:
        if used[row[0]]:
            continue
        take = row[~used[row]]
        if len(take) < k:
            continue
        used[take] = True
        pools.append(take)
    return pools


def metacell_sweep(df, emb, rng):
    """
    Rho against metacell size, with pooling at random as the control.

    The control is matched on pool count as well as pool size: the kNN cover
    leaves cells over, so drawing floor(n/k) random pools instead would give the
    control more pools than the real arm and make the two curves incomparable.
    """
    samples = df["sample"].values
    c5, c6 = df["CEACAM5"].values, df["CEACAM6"].values
    rows = []
    for k in KS:
        m5, m6, s5, s6 = [], [], [], []
        for s in sorted(set(samples)):
            idx = np.where(samples == s)[0]
            if len(idx) < k:
                continue
            if k == 1:
                pools = [np.array([i]) for i in range(len(idx))]
            else:
                pools = _pools_knn(emb[idx], k)
            if not pools:
                continue
            m5 += [c5[idx[p]].mean() for p in pools]
            m6 += [c6[idx[p]].mean() for p in pools]
            perm = rng.permutation(len(idx))
            sh = [perm[i * k:(i + 1) * k] for i in range(len(pools))]
            s5 += [c5[idx[p]].mean() for p in sh]
            s6 += [c6[idx[p]].mean() for p in sh]
        rho, p = stats.spearmanr(m5, m6)
        rho_s, _ = stats.spearmanr(s5, s6)
        rows.append(dict(k=k, n_pools=len(m5), rho=float(rho), p=float(p),
                         rho_random_pools=float(rho_s)))
    agg = df.groupby("sample")[["CEACAM5", "CEACAM6"]].mean()
    rho, p = stats.spearmanr(agg["CEACAM5"], agg["CEACAM6"])
    rows.append(dict(k=np.nan, n_pools=len(agg), rho=float(rho), p=float(p),
                     rho_random_pools=np.nan))
    return pd.DataFrame(rows)


def panel(strata, ws, sweep):
    d = S8 / "S8_H"
    d.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(17.0 * SCALE * CM, 5.0 * SCALE * CM))
    fs = dict(fontsize=6 * SCALE)
    tk = dict(labelsize=5.5 * SCALE, width=0.8, length=3)

    ax = axes[0]
    x = np.arange(len(strata))
    ax.plot(x, strata["rho"], "-o", color=COLOR_REAL, lw=1.2 * SCALE * 0.5,
            ms=4 * SCALE * 0.5, label="Spearman $\\rho$, per cell")
    ax.plot(x, strata["pct_double_of_expressing"] / 100, "-s", color=COLOR_5,
            lw=1.2 * SCALE * 0.5, ms=4 * SCALE * 0.5,
            label="double positive, of cells\nexpressing either gene")
    ax.plot(x, strata["pct_single_pos"] / 100, "-^", color=COLOR_6,
            lw=1.2 * SCALE * 0.5, ms=4 * SCALE * 0.5,
            label="single positive, of all cells")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:,.0f}" for v in strata["median_total_counts"]],
                       rotation=0)
    ax.set_xlabel("median UMI per cell, quintile", **fs)
    ax.set_ylabel("proportion, or $\\rho$", **fs)
    ax.set_title("Deeper cells look more double positive", **fs)
    ax.legend(fontsize=4.6 * SCALE, frameon=False, loc="upper left")
    ax.set_ylim(0, 1)

    # The odds ratio rather than observed/expected: the expected value itself
    # climbs with detection rate, so the ratio shrinks even as the association
    # strengthens. The odds ratio does not have that dependence.
    ax = axes[1]
    ax.plot(x, strata["odds_ratio"], "-o", color=COLOR_REAL,
            lw=1.2 * SCALE * 0.5, ms=4 * SCALE * 0.5)
    ax.axhline(1.0, color="#999999", lw=0.8 * SCALE * 0.5, ls="--")
    ax.annotate("independence", (x[0], 1.0), textcoords="offset points",
                xytext=(2, 4 * SCALE * 0.5), fontsize=4.6 * SCALE,
                color="#777777")
    # Labels sit below the line and alternate side, so they clear the y axis
    # on the left and the last point on the right.
    for i, r in strata.reset_index().iterrows():
        last = i == len(strata) - 1
        ax.annotate(f"{r['obs_over_expected']:.2f}\u00d7 expected",
                    (i, r["odds_ratio"]), textcoords="offset points",
                    xytext=(-4 * SCALE * 0.5 if last else 4 * SCALE * 0.5,
                            -9 * SCALE * 0.5),
                    ha="right" if last else "left", va="top",
                    fontsize=4.6 * SCALE)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:,.0f}" for v in strata["median_total_counts"]])
    ax.set_yscale("log")
    ax.set_yticks([1, 2, 5, 10, 20])
    ax.set_yticklabels(["1", "2", "5", "10", "20"])
    ax.set_ylim(0.8, 30)
    ax.set_xlabel("median UMI per cell, quintile", **fs)
    ax.set_ylabel("odds ratio for co-detection", **fs)
    ax.set_title("Co-detected far above independence, at every depth", **fs)

    ax = axes[2]
    s = sweep.dropna(subset=["k"])
    ax.plot(s["k"], s["rho"], "-o", color=COLOR_REAL, lw=1.2 * SCALE * 0.5,
            ms=4 * SCALE * 0.5, label="kNN metacells")
    ax.plot(s["k"], s["rho_random_pools"], "-o", color=COLOR_SHUF,
            lw=1.2 * SCALE * 0.5, ms=4 * SCALE * 0.5, label="random pools, same size")
    lvl = float(sweep.loc[sweep["k"].isna(), "rho"].iloc[0])
    ax.axhline(lvl, color=COLOR_5, lw=0.8 * SCALE * 0.5, ls=":")
    ax.annotate(f"whole-sample means, $\\rho$ = {lvl:.2f}", (KS[1], lvl),
                textcoords="offset points", xytext=(0, -11 * SCALE * 0.5),
                fontsize=4.6 * SCALE, color=COLOR_5)
    ax.set_xscale("log")
    ax.set_xticks(KS)
    ax.set_xticklabels([str(k) for k in KS])
    ax.set_xlabel("cells pooled per metacell", **fs)
    ax.set_ylabel("Spearman $\\rho$", **fs)
    ax.set_title("Pooling raises $\\rho$, neighbours or not", **fs)
    ax.legend(fontsize=4.6 * SCALE, frameon=False, loc="lower right")
    ax.set_ylim(0, 1)

    for ax in axes:
        ax.tick_params(axis="both", **tk)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.suptitle(
        f"Pre-treatment stomach epithelial cells (n = {int(strata['n_cells'].sum()):,}"
        f", {len(ws)} samples); expression from the log1p CP10K matrix",
        fontsize=6 * SCALE, y=0.99)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.80, bottom=0.15, wspace=0.32)
    stem = d / "S8_H_dropout_and_coexpression"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}.[svg|pdf|png]")


def main():
    rng = np.random.default_rng(SEED)
    df, emb, ids = load()
    L = ["CEACAM5/CEACAM6 CO-EXPRESSION AND DROPOUT - Reviewer 1 point R1.5",
         "=" * 90, "",
         f"Pre-treatment stomach epithelial cells: {len(df):,} "
         f"from {df['sample'].nunique()} samples", f"seed = {SEED}", ""]

    strata = depth_strata(df)
    strata.to_csv(OUT / "dropout_depth_strata.csv", index=False)
    L += ["A. SEQUENCING DEPTH", "-" * 90,
          strata[["stratum", "n_cells", "median_total_counts", "rho",
                  "pct_double_pos", "pct_single_pos", "pct_double_neg",
                  "pct_double_of_expressing"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"), "",
          "B. CO-DETECTION AGAINST INDEPENDENCE, WITHIN STRATUM", "-" * 90,
          strata[["stratum", "pct_double_pos", "expected_double_pos",
                  "obs_over_expected", "odds_ratio", "p_fisher"]]
          .to_string(index=False, float_format=lambda v: f"{v:.4g}"), ""]

    ws = within_sample(df, rng, ids)
    ws.to_csv(OUT / "coexpression_within_sample.csv", index=False)

    # The pooled per-cell correlation, and the part of it that survives when
    # only between-patient structure is left: CEACAM6 is permuted within each
    # sample, so every sample keeps its own mean for both genes while no cell
    # keeps its partner. Whatever the pooled rho retains under that permutation
    # is what patient-level covariation alone can produce.
    rho_pooled, _ = stats.spearmanr(df["CEACAM5"], df["CEACAM6"])
    perm6 = df.groupby("sample")["CEACAM6"].transform(
        lambda v: v.values[rng.permutation(len(v))])
    rho_floor, _ = stats.spearmanr(df["CEACAM5"], perm6)
    L += ["C. WITHIN SAMPLE, PER CELL", "-" * 90,
          f"pooled across samples, per cell            rho {rho_pooled:.4f}",
          f"  the same, CEACAM6 permuted within sample rho {rho_floor:.4f}"
          "   <- what between-patient differences alone give",
          f"median rho {ws['rho'].median():.4f}   "
          f"IQR {ws['rho'].quantile(.25):.4f}-{ws['rho'].quantile(.75):.4f}   "
          f"range {ws['rho'].min():.4f}-{ws['rho'].max():.4f}",
          f"positive in {int((ws['rho'] > 0).sum())} of {len(ws)} samples",
          f"within-sample permutation null, median {ws['rho_permuted'].median():.4f}",
          "", ws.to_string(index=False, float_format=lambda v: f"{v:.4f}"), ""]

    sweep = metacell_sweep(df, emb, rng)
    sweep.to_csv(OUT / "metacell_sweep.csv", index=False)
    L += ["D. METACELL SIZE, AGAINST POOLING AT RANDOM", "-" * 90,
          sweep.to_string(index=False, float_format=lambda v: f"{v:.4f}"), "",
          "Past a small k the two curves coincide, so the rise of rho with pool",
          "size is the arithmetic of averaging rather than recovered signal. The",
          "metacell view is a presentation of the data; the co-expression claim",
          "rests on B, which is computed inside a depth stratum and is unaffected",
          "by pooling.", ""]

    (OUT / "dropout_coexpression_report.txt").write_text("\n".join(L),
                                                         encoding="utf-8")
    print("\n".join(L))
    panel(strata, ws, sweep)


if __name__ == "__main__":
    main()
