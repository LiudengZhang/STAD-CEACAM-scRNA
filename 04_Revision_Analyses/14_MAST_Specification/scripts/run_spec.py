"""
Task 2: run one differential-expression specification for one cell type and
one timepoint, then the project's own prerank GSEA on the result.

Copied from 12_R1.8_DEG_Recompute/scripts/recompute_deg.py. The gene filter,
the contrast, the seed, the Hallmark file, the rank metric and the GSEA
settings are the original's, unchanged. What varies is the model, and only
the model:

  dep  zlm(~ condition + sample_id + cngeneson, sca)
         the deposited call, kept here so its GSEA is produced by exactly the
         same code path as the alternatives and the comparison is not
         confounded by anything else.

  A0   zlm(~ condition, sca)
         neither sample_id nor cngeneson. Not a defensible analysis - the
         cellular detection rate is a real confounder and MAST's own
         documentation says to adjust for it - but it is the only way to see
         how much of the difference between MAST and the Welch t-test is the
         cngeneson term and how much is MAST.

  A    zlm(~ condition + cngeneson, sca)
         sample_id dropped. Pseudoreplication untreated - every cell counts as
         an independent observation - but the condition coefficient is
         estimable.

  B    zlm(~ condition + cngeneson + (1 | sample_id), sca,
           method = "glmer", ebayes = FALSE)
         a random intercept for sample. MAST's own documented mixed-model
         call. ebayes must be FALSE: MAST's empirical-Bayes variance shrinkage
         is defined for the fixed-effect fits only.

  T    Welch t-test, scanpy rank_genes_groups, on exactly the same cells and
         genes. Not a specification of MAST; the comparator.

Sample identifiers are replaced by a positional label (S01, S02, ...) before
anything is assigned into R, so no specimen identifier reaches R, a log or a
file. The labelling is per contrast and carries no meaning across contrasts.

Run:
    LD_LIBRARY_PATH=$CONDA_PREFIX/lib conda run -n Liudeng_Python_310 \
        python run_spec.py --cell MoMac --phase post --spec B [--threads 8]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from sources import (CELL_SOURCES, PHASES, HALLMARK_GMT, SEED,
                     MIN_CELLS_PER_GROUP, load_subset, filter_genes)

OUT = Path(__file__).resolve().parents[1] / "outputs"
(OUT / "deg").mkdir(parents=True, exist_ok=True)
(OUT / "gsea").mkdir(parents=True, exist_ok=True)
(OUT / "summaries").mkdir(parents=True, exist_ok=True)

FORMULA = {
    "dep": "~ condition + sample_id + cngeneson",
    "A0":  "~ condition",
    "A":   "~ condition + cngeneson",
    "B":   "~ condition + cngeneson + (1 | sample_id)",
}


def run_mast(expr, groups, samples, genes, ref, test, spec, threads):
    """MAST, as MAST documents itself. No fallback: a failure is a failure."""
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter

    bad = ~np.isfinite(expr)
    if bad.any():
        raise SystemExit(f"{int(bad.sum())} non-finite values in the matrix.")

    def chars(a):
        return np.asarray([str(x) for x in a], dtype=object)

    uniq = sorted(set(samples))
    anon = {s: f"S{i+1:02d}" for i, s in enumerate(uniq)}
    samples_anon = np.array([anon[s] for s in samples])

    with localconverter(ro.default_converter + pandas2ri.converter
                        + numpy2ri.converter):
        ro.r.assign("expr_matrix", expr)
        ro.r.assign("cell_metadata", chars(groups))
        ro.r.assign("cell_names", chars([f"c{i}" for i in range(len(groups))]))
        ro.r.assign("gene_names", chars(genes))
        ro.r.assign("reference_group", str(ref))
        ro.r.assign("comparison_group", str(test))
        ro.r.assign("sample_id", chars(samples_anon))
        ro.r.assign("spec", str(spec))
        ro.r.assign("n_threads", int(threads))

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
    contrast_name <- paste0("condition", as.character(comparison_group))

    # MAST's own parallelism, and the only one it has: zlm() dispatches to
    # parallel::mclapply when getOption("mc.cores") > 1 and falls back to a
    # plain lapply otherwise (see body(MAST::zlm)). A registered foreach /
    # doParallel backend is NOT used by MAST and leaves every fit
    # single-threaded while looking like it is running on N cores.
    options(mc.cores = n_threads)
    set.seed(42)
    zlmCond <- switch(spec,
      "dep" = zlm(~ condition + sample_id + cngeneson, sca),
      "A0"  = zlm(~ condition, sca),
      "A"   = zlm(~ condition + cngeneson, sca),
      "B"   = zlm(~ condition + cngeneson + (1 | sample_id), sca,
                  method = "glmer", ebayes = FALSE),
      stop("unknown spec"))
    conv_C <- sum(zlmCond@converged[, "C"], na.rm = TRUE)
    conv_D <- sum(zlmCond@converged[, "D"], na.rm = TRUE)
    n_fit  <- nrow(zlmCond@converged)
    na_coefC <- sum(is.na(zlmCond@coefC[, contrast_name]))
    na_coefD <- sum(is.na(zlmCond@coefD[, contrast_name]))
    cat(sprintf("[R] spec=%s genes=%d convergedC=%d convergedD=%d naC=%d naD=%d\\n",
                spec, n_fit, conv_C, conv_D, na_coefC, na_coefD))

    # summary()'s own `parallel` argument defaults to FALSE and is handed
    # straight to lrTest(), so the likelihood-ratio refit - which costs about
    # as much as the fit itself - runs on one core no matter what mc.cores
    # says. This is MAST's documented argument for that, and it changes no
    # number: the reduced-model refit is the same per-gene computation either
    # way. Checked, not assumed - see EQUIVALENCE.md.
    s <- summary(zlmCond, doLRT = contrast_name, parallel = n_threads > 1)$datatable
    fc <- merge(s[contrast==contrast_name & component=="H", .(primerid, `Pr(>Chisq)`)],
                s[contrast==contrast_name & component=="logFC", .(primerid, coef)],
                by="primerid")
    results <- data.frame(gene = fc$primerid,
                          logfoldchanges = fc$coef,
                          pvals = fc$`Pr(>Chisq)`,
                          pvals_adj = p.adjust(fc$`Pr(>Chisq)`, method="fdr"),
                          stringsAsFactors = FALSE)
    results <- results[!is.na(results$pvals) & !is.na(results$logfoldchanges), ]
    fit_report <- data.frame(n_genes_fit = n_fit, converged_C = conv_C,
                             converged_D = conv_D, na_coefC = na_coefC,
                             na_coefD = na_coefD)
''')
    with localconverter(ro.default_converter + pandas2ri.converter
                        + numpy2ri.converter):
        out = ro.conversion.rpy2py(ro.r["results"])
        rep = ro.conversion.rpy2py(ro.r["fit_report"])
    rep = rep if isinstance(rep, pd.DataFrame) else pd.DataFrame(rep)
    return (out if isinstance(out, pd.DataFrame) else pd.DataFrame(out)), rep


def run_ttest(raw, keep, cfg):
    sub = raw[:, keep].copy()
    sub.obs["_g"] = sub.obs[cfg["column"]].astype(str)
    sc.tl.rank_genes_groups(sub, groupby="_g", groups=[cfg["test"]],
                            reference=cfg["ref"], method="t-test_overestim_var",
                            use_raw=False)
    d = sc.get.rank_genes_groups_df(sub, group=cfg["test"])
    return d.rename(columns={"names": "gene"})


def gsea(df):
    """The project's own prerank: rank_metric = logFC x -log10(P), Hallmark,
    1000 permutations. recompute_deg.py:343, unchanged."""
    import gseapy as gp
    d = df.dropna(subset=["logfoldchanges", "pvals"]).copy()
    d["metric"] = d["logfoldchanges"] * -np.log10(d["pvals"].clip(lower=1e-300))
    rnk = (d.set_index("gene")["metric"].groupby(level=0).first()
           .sort_values(ascending=False))
    if len(rnk) < 15:
        return None
    res = gp.prerank(rnk=rnk, gene_sets=str(HALLMARK_GMT), permutation_num=1000,
                     min_size=15, max_size=500, seed=SEED, no_plot=True,
                     outdir=None, threads=2)
    out = res.res2d.copy()
    out["Term"] = out["Term"].astype(str).str.replace(r"^.*__", "", regex=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--phase", required=True, choices=list(PHASES))
    ap.add_argument("--spec", required=True, choices=["dep", "A0", "A", "B", "T"])
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--force", action="store_true",
                    help="redo a job whose summary already exists")
    a = ap.parse_args()
    cfg = PHASES[a.phase]
    t0 = time.time()

    # Idempotent: a job whose summary already exists is not redone. This is
    # what makes it safe to run a second driver over the tail of a queue the
    # first driver has not reached yet - whichever gets there first wins, and
    # the other exits without touching the outputs. The run is seeded, so a
    # collision would produce the same numbers anyway.
    done = OUT / "summaries" / f"{a.cell}_{a.phase}_{a.spec}.csv"
    if done.exists() and not a.force:
        print(f"[{a.cell} {a.phase} {a.spec}] already done, skipping", flush=True)
        return

    raw = load_subset(a.cell, a.phase)
    if raw is None or raw.n_obs == 0:
        raise SystemExit(f"{a.cell} {a.phase}: no cells")
    g = raw.obs[cfg["column"]].astype(str)
    n_ref, n_test = int((g == cfg["ref"]).sum()), int((g == cfg["test"]).sum())
    if min(n_ref, n_test) < MIN_CELLS_PER_GROUP:
        raise SystemExit(f"{a.cell} {a.phase}: too few cells {n_ref}/{n_test}")

    X = raw.X.toarray() if hasattr(raw.X, "toarray") else np.asarray(raw.X)
    keep = filter_genes(X)
    genes = np.asarray(raw.var_names)[keep]
    print(f"[{a.cell} {a.phase} {a.spec}] cells {raw.n_obs} "
          f"({n_ref} R / {n_test} NR) genes {int(keep.sum())}", flush=True)

    row = dict(cell=a.cell, phase=a.phase, spec=a.spec,
               formula=FORMULA.get(a.spec, "Welch t-test"),
               source=CELL_SOURCES[a.cell].name, n_ref_cells=n_ref,
               n_test_cells=n_test, n_genes_tested=int(keep.sum()),
               n_samples=int(raw.obs["sample"].astype(str).nunique()))

    if a.spec == "T":
        df = run_ttest(raw, keep, cfg)
    else:
        expr = np.ascontiguousarray(X[:, keep], dtype=np.float64)
        samples = raw.obs["sample"].astype(str).values
        df, rep = run_mast(expr, g.values, samples, genes,
                           cfg["ref"], cfg["test"], a.spec, a.threads)
        for c in rep.columns:
            row[c] = int(rep.iloc[0][c])
    df.to_csv(OUT / "deg" / f"{a.cell}_{a.phase}_{a.spec}.csv", index=False)
    row["n_genes_returned"] = int(len(df))

    e = gsea(df)
    if e is not None:
        e.to_csv(OUT / "gsea" / f"{a.cell}_{a.phase}_{a.spec}_hallmark.csv",
                 index=False)
        hit = e[e["Term"].str.contains("TNF-alpha", case=False, na=False)]
        if len(hit):
            h = hit.iloc[0]
            row.update(nfkb_nes=float(h["NES"]), nfkb_nom_p=float(h["NOM p-val"]),
                       nfkb_fdr_q=float(h["FDR q-val"]))
            e2 = e.copy()
            e2["NES"] = pd.to_numeric(e2["NES"], errors="coerce")
            e2 = e2.sort_values("NES", ascending=False).reset_index(drop=True)
            row["nfkb_rank"] = int(e2.index[e2["Term"] == h["Term"]][0]) + 1
            row["n_sets"] = int(len(e2))
    row["minutes"] = round((time.time() - t0) / 60, 2)
    row["status"] = "ok"
    pd.DataFrame([row]).to_csv(
        OUT / "summaries" / f"{a.cell}_{a.phase}_{a.spec}.csv", index=False)
    print(f"[{a.cell} {a.phase} {a.spec}] ok "
          f"NES={row.get('nfkb_nes')} {row['minutes']}min", flush=True)


if __name__ == "__main__":
    main()
