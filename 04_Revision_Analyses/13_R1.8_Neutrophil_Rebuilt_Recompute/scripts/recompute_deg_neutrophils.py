"""
Neutrophil differential expression and Hallmark GSEA, NR vs R at each
timepoint, on paths.NEUTROPHILS_SOUND_H5AD. Same method as
12_R1.8_DEG_Recompute/work/recompute_deg.py (model, gene filter, contrast,
seed, Hallmark file and GSEA settings).

Writes to this module's outputs/: deg/<cell>_<phase>_{mast,ttest}.csv,
gsea/<cell>_<phase>_<method>_hallmark.csv and
recompute_summary_<cell>_<phase>.csv.

Run: python recompute_deg_neutrophils.py [--cell CELL|ALL] [--phase pre|post]
(scripts/run.sh sets up LD_LIBRARY_PATH for MAST and runs both phases.)
"""

from pathlib import Path
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import (CELL_TYPE_H5AD, HALLMARK_GMT,  # noqa: E402
                   NEUTROPHILS_SOUND_H5AD, PROJECT_ROOT)

warnings.filterwarnings("ignore")

OUT = Path(__file__).resolve().parents[1] / "outputs"
(OUT / "deg").mkdir(parents=True, exist_ok=True)
(OUT / "gsea").mkdir(parents=True, exist_ok=True)

SEED = 42
MIN_PCT = 0.10          # a gene must be detected in this fraction of the cells
MIN_CELLS_PER_GROUP = 20
HALLMARK = str(HALLMARK_GMT)
NFKB_TERM = "TNF-alpha Signaling via NF-kB"

# Mast and plasma cells come from separate integrated objects.
ROUND_4_MAJOR = (PROJECT_ROOT.parent / "upstream-pipeline" / "04_Final_Panels"
                 / "00_Set_Ups" / "00_Data" / "01_Major_Cell_Types")

# The thirteen populations and their input files (one population per file).
CELL_SOURCES = {
    "B_cells": CELL_TYPE_H5AD["B_cells"],
    "DC_cells": CELL_TYPE_H5AD["DC_cells"],
    "Endothelial_cells": CELL_TYPE_H5AD["Endothelial"],
    "Epithelial": CELL_TYPE_H5AD["Epithelial"],
    "Fibroblast": CELL_TYPE_H5AD["Fibroblast"],
    "Mast_cells": ROUND_4_MAJOR / "mast_cell_integrated.h5ad",
    "MoMac": CELL_TYPE_H5AD["MoMac"],
    # Built from the per-sample counts (matched on Ensembl ID); see
    # 06_Clean_Data/rebuild_neutrophils_sound.py.
    "Neutrophils": NEUTROPHILS_SOUND_H5AD,
    "NK_cells": CELL_TYPE_H5AD["NK_cells"],
    "Pericyte": CELL_TYPE_H5AD["Pericyte"],
    "Plasma_cells": ROUND_4_MAJOR / "plasma_cell_integrated.h5ad",
    "TCD4_cells": CELL_TYPE_H5AD["TCD4"],
    "TCD8_cells": CELL_TYPE_H5AD["TCD8"],
}

# Non-responders vs responders (reference): positive logFC = higher in NR.
PHASES = {
    "pre": dict(column="stomach_pre_grouping", ref="Responsed", test="No-response"),
    "post": dict(column="stomach_post_grouping", ref="Responsed", test="No-response"),
}


LADDER_CELLS = 40
LADDER_TOL = 0.01       # log1p CP10K matrices sit near 3e-4 (float32 precision)


class BadMatrix(RuntimeError):
    """The input for one cell type is not log1p CP10K. Raised instead of
    SystemExit so main() can record it per cell type."""


def assert_log1p_cp10k(X, source):
    """Refuse a matrix that is not log1p CP10K, by the integer-ladder test:
    within one cell, expm1(value) / expm1(smallest value) must be 1, 2, 3, ...
    """
    Xc = X.tocsr() if sparse.issparse(X) else np.asarray(X)
    rng = np.random.default_rng(SEED)
    rows = rng.choice(Xc.shape[0], min(LADDER_CELLS, Xc.shape[0]), replace=False)
    worst, tested, integral = 0.0, 0, 0
    for i in rows:
        v = (Xc.data[Xc.indptr[i]:Xc.indptr[i + 1]] if sparse.issparse(Xc)
             else Xc[i][Xc[i] > 0]).astype(np.float64)
        if v.size < 20:
            continue
        tested += 1
        if np.abs(v - np.round(v)).max() < 1e-6:
            integral += 1
            continue
        c = np.expm1(v) / np.expm1(v.min())      # one count is the smallest step
        worst = max(worst, float(np.abs(c - np.round(c)).max()))
    if tested == 0:
        raise BadMatrix(f"{source}: too few detected genes to test the matrix")
    if integral > tested / 2:
        raise BadMatrix(
            f"{source}: .raw holds integer counts, not log1p CP10K. "
            "Normalising it here would hide which file was wrong; fix the input.")
    if worst >= LADDER_TOL:
        raise BadMatrix(
            f"{source}: .raw is NOT log1p CP10K - implied counts miss the "
            f"integers by up to {worst:.3f} (correct matrices sit at 3e-4). This "
            "is the signature of a second normalisation on already normalised "
            "values. "
            "Refusing to test differential expression on it.")


def load_subset(path, phase):
    """The cells of one population at one timepoint, expression from .raw.

    The file holds one population, so the only filter is the timepoint's
    response grouping.
    """
    cfg = PHASES[phase]
    adata = sc.read_h5ad(path)
    if adata.raw is None:
        raise BadMatrix(f"{path.name}: no .raw - refusing to fall back to .X")
    m = adata.obs[cfg["column"]].astype(str).isin([cfg["ref"], cfg["test"]])
    if m.sum() == 0:
        return None
    sub = adata[m]
    raw = sub.raw.to_adata()[sub.obs_names]
    raw.obs = sub.obs.copy()
    assert_log1p_cp10k(raw.X, f"{path.name} .raw")
    return raw


def filter_genes(X):
    """Detected in at least MIN_PCT of the cells."""
    keep = (np.asarray((X > 0).mean(axis=0)).ravel() >= MIN_PCT)
    return keep


def run_mast(expr, groups, samples, genes, ref, test):
    """MAST hurdle model via rpy2. Any NaN or infinity aborts."""
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter

    bad = ~np.isfinite(expr)
    if bad.any():
        raise SystemExit(
            f"{int(bad.sum())} non-finite values in the expression matrix. "
            "The input is not log1p CP10K; refusing to substitute zeros.")

    # Pass strings as dtype=object: rpy2 turns a numpy '<U' array into an R
    # array, which MAST's FromMatrix rejects.
    def chars(a):
        return np.asarray([str(x) for x in a], dtype=object)

    with localconverter(ro.default_converter + pandas2ri.converter
                        + numpy2ri.converter):
        ro.r.assign("expr_matrix", expr)
        ro.r.assign("cell_metadata", chars(groups))
        ro.r.assign("cell_names", chars([f"c{i}" for i in range(len(groups))]))
        ro.r.assign("gene_names", chars(genes))
        ro.r.assign("reference_group", str(ref))
        ro.r.assign("comparison_group", str(test))
        ro.r.assign("sample_id", chars(samples))

    ro.r('''
    suppressPackageStartupMessages({library(MAST); library(SingleCellExperiment)})
    expr_mat <- t(expr_matrix)
    rownames(expr_mat) <- as.character(gene_names)
    colnames(expr_mat) <- as.character(cell_names)
    cell_meta <- data.frame(condition = as.character(cell_metadata),
                            sample_id = as.factor(as.character(sample_id)),
                            wellKey = as.character(cell_names),
                            stringsAsFactors = FALSE)
    rownames(cell_meta) <- as.character(cell_names)
    cell_meta$cngeneson <- as.numeric(scale(colSums(expr_mat > 0)))
    gene_meta <- data.frame(primerid = as.character(gene_names),
                            stringsAsFactors = FALSE)
    rownames(gene_meta) <- as.character(gene_names)
    sca <- FromMatrix(exprsArray = expr_mat, cData = cell_meta, fData = gene_meta)
    sca$condition <- relevel(factor(sca$condition), ref = as.character(reference_group))
    zlmCond <- zlm(~ condition + sample_id + cngeneson, sca)
    contrast_name <- paste0("condition", as.character(comparison_group))
    s <- summary(zlmCond, doLRT = contrast_name)$datatable
    fc <- merge(s[contrast==contrast_name & component=="H", .(primerid, `Pr(>Chisq)`)],
                s[contrast==contrast_name & component=="logFC", .(primerid, coef)],
                by="primerid")
    results <- data.frame(gene = fc$primerid,
                          logfoldchanges = fc$coef,
                          pvals = fc$`Pr(>Chisq)`,
                          pvals_adj = p.adjust(fc$`Pr(>Chisq)`, method="fdr"),
                          stringsAsFactors = FALSE)
    results <- results[!is.na(results$pvals) & !is.na(results$logfoldchanges), ]
    ''')
    with localconverter(ro.default_converter + pandas2ri.converter
                        + numpy2ri.converter):
        out = ro.conversion.rpy2py(ro.r["results"])
    return out if isinstance(out, pd.DataFrame) else pd.DataFrame(out)


def run_ttest(raw, keep, cfg):
    """Welch t-test on exactly the cells and genes MAST just used."""
    sub = raw[:, keep].copy()
    sub.obs["_g"] = sub.obs[cfg["column"]].astype(str)
    sc.tl.rank_genes_groups(sub, groupby="_g", groups=[cfg["test"]],
                            reference=cfg["ref"], method="t-test_overestim_var",
                            use_raw=False)
    d = sc.get.rank_genes_groups_df(sub, group=cfg["test"])
    return d.rename(columns={"names": "gene", "pvals_adj": "pvals_adj"})


def gsea(df, tag):
    """Prerank on logFC x -log10(P). Positive NES means enriched in NR."""
    import gseapy as gp
    d = df.dropna(subset=["logfoldchanges", "pvals"]).copy()
    d["metric"] = d["logfoldchanges"] * -np.log10(d["pvals"].clip(lower=1e-300))
    rnk = (d.set_index("gene")["metric"].groupby(level=0).first()
           .sort_values(ascending=False))
    if len(rnk) < 15:
        return None
    res = gp.prerank(rnk=rnk, gene_sets=HALLMARK, permutation_num=1000,
                     min_size=15, max_size=500, seed=SEED, no_plot=True,
                     outdir=None, threads=2)
    out = res.res2d.copy()
    out["Term"] = out["Term"].astype(str).str.replace(r"^.*__", "", regex=True)
    return out


def one(cell, phase):
    cfg = PHASES[phase]
    t0 = time.time()
    raw = load_subset(CELL_SOURCES[cell], phase)
    if raw is None or raw.n_obs == 0:
        return dict(cell=cell, phase=phase, status="no cells")
    g = raw.obs[cfg["column"]].astype(str)
    n_ref, n_test = int((g == cfg["ref"]).sum()), int((g == cfg["test"]).sum())
    if min(n_ref, n_test) < MIN_CELLS_PER_GROUP:
        return dict(cell=cell, phase=phase, status=f"too few cells {n_ref}/{n_test}")

    X = raw.X.toarray() if hasattr(raw.X, "toarray") else np.asarray(raw.X)
    keep = filter_genes(X)
    genes = np.asarray(raw.var_names)[keep]
    expr = np.ascontiguousarray(X[:, keep], dtype=np.float64)
    samples = raw.obs["sample"].astype(str).values

    # MAST is the sensitivity analysis: if R fails to load, record the failure
    # and still run the t-test.
    mast, mast_status = None, "ok"
    try:
        mast = run_mast(expr, g.values, samples, genes, cfg["ref"], cfg["test"])
        mast.to_csv(OUT / "deg" / f"{cell}_{phase}_mast.csv", index=False)
    except Exception as exc:
        mast_status = f"ERROR {type(exc).__name__}: {exc}"[:200]
        print(f"[{cell} {phase}] MAST unavailable: {mast_status}", flush=True)

    tt = run_ttest(raw, keep, cfg)
    tt.to_csv(OUT / "deg" / f"{cell}_{phase}_ttest.csv", index=False)

    row = dict(cell=cell, phase=phase, status="ok", mast_status=mast_status,
               source=CELL_SOURCES[cell].name, n_ref=n_ref, n_test=n_test,
               n_genes_tested=int(keep.sum()), n_genes_total=int(len(keep)),
               pct_genes_kept=100 * float(keep.mean()),
               minutes=round((time.time() - t0) / 60, 1))
    for name, df in (("mast", mast), ("ttest", tt)):
        if df is None:
            continue
        e = gsea(df, f"{cell}_{phase}_{name}")
        if e is None:
            continue
        e.to_csv(OUT / "gsea" / f"{cell}_{phase}_{name}_hallmark.csv", index=False)
        hit = e[e["Term"].str.contains("TNF-alpha", case=False, na=False)]
        if len(hit):
            row[f"{name}_nfkb_nes"] = float(hit.iloc[0]["NES"])
            row[f"{name}_nfkb_fdr"] = float(hit.iloc[0]["FDR q-val"])
    return row


def main():
    ap = argparse.ArgumentParser()
    # Neutrophils by default; the other twelve are computed in module 12.
    ap.add_argument("--cell", default="Neutrophils",
                    choices=sorted(CELL_SOURCES) + ["ALL"])
    ap.add_argument("--phase", default=None, choices=list(PHASES))
    args = ap.parse_args()

    cells = list(CELL_SOURCES) if args.cell == "ALL" else [args.cell]
    phases = [args.phase] if args.phase else list(PHASES)
    missing = [c for c in cells if not CELL_SOURCES[c].exists()]
    if missing:
        raise SystemExit("no input file for: "
                         + ", ".join(f"{c} ({CELL_SOURCES[c]})" for c in missing))
    rows = []
    for c in cells:
        for p in phases:
            print(f"[{c} {p}] start from {CELL_SOURCES[c].name}", flush=True)
            try:
                r = one(c, p)
            except Exception as exc:
                r = dict(cell=c, phase=p, status=f"ERROR {type(exc).__name__}: {exc}")
            rows.append(r)
            print(f"[{c} {p}] {r.get('status')} "
                  f"{r.get('minutes','')}min", flush=True)

    df = pd.DataFrame(rows)
    # One summary file per (cell, phase).
    for r in rows:
        pd.DataFrame([r]).to_csv(
            OUT / f"recompute_summary_{r['cell']}_{r['phase']}.csv", index=False)
    print(df.to_string(index=False))

    # Exit non-zero if any job failed; "no cells" and "too few cells" are
    # expected skips.
    expected_skips = ("no cells", "too few cells")
    broken = [r for r in rows
              if r.get("status") != "ok"
              and not str(r.get("status", "")).startswith(expected_skips)]
    no_mast = [r for r in rows if r.get("mast_status", "ok") != "ok"]
    if no_mast:
        print(f"\nMAST did not run for {len(no_mast)} of {len(rows)} jobs; the "
              f"t-test results the paper reports are unaffected, and the "
              f"reason is in the mast_status column.")
    if broken:
        print(f"\n{len(broken)} of {len(rows)} jobs failed:")
        for r in broken:
            print(f"  {r['cell']:<20} {r['phase']:<5} {r['status']}")
        sys.exit(1)


if __name__ == "__main__":
    main()
