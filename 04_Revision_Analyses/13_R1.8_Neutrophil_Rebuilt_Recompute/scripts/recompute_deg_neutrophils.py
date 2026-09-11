"""
12_R1.8_DEG_Recompute/scripts/recompute_deg.py, run for neutrophils only, on
the rebuilt input.

The original is copied rather than edited, and exactly two things differ:

  1. `CELL_SOURCES["Neutrophils"]` points at
     06_Clean_Data/02_Rebuilt/Neutrophils_sound.h5ad instead of the shipped
     submission-tree file. Nothing else about the source changes: the rebuilt object
     holds the same 61,167 cells and the same 20,060 genes, in the same order,
     with the same obs. Only the values in `.raw` differ, because they are
     normalised once instead of twice.
  2. Outputs land in this module's own `outputs/`, because OUT is derived from
     __file__. 12_R1.8_DEG_Recompute/outputs/ - which verify_numbers.py reads -
     is not touched.

The model, the gene filter, the contrast, the seed, the Hallmark file and the
GSEA settings are the original's, unchanged. The point is to move one input and
nothing else.

The docstring of the original follows.

Reviewer 1 point R1.8, and the repair of the differential expression it rests on.

R1.8 says the claim that "every population exhibited elevated NF-kB pathway
activity" is "unusually broad and may reflect a common inflammatory state,
sample composition, or an analytical artifact". The audit in
00_Data_Audit/FINDINGS.md section 7 found that the third possibility was the
case, twice over:

  - Seven cell types were tested on the `.X` destroyed on 2025-07-30 (section 1).
  - The other six were destroyed at run time, because the DEG library carried
    the same misfiring guard - `if 'log1p' not in adata.uns` - and log1p'd a
    z-scored matrix (`deg_pathway_analysis_library.py:162`).
  - Neither was visible, because line 521 replaced every NaN with a zero before
    MAST saw it.

So no cell type was ever tested on log1p CP10K data. This module runs the same
test again on the matrix that was always intact.

Three things are deliberately different from the original:

  1. Expression comes from `.raw` - log1p CP10K - and the script refuses to
     start if `.raw` is absent. There is no `.X` branch.
  2. A NaN or an infinity is a fatal error, not something to fill with zero.
     The original's zero-fill is what made the damage invisible for a year.
  3. Every cell type is read from its own per-cell-type object, and the matrix
     is tested before it is used.

Point 3 is not a detail. Taking every cell type from one pooled object,
`full_dataset.h5ad`, on the belief that its `.raw` is clean, is the obvious way
to write this and it is wrong: FINDINGS.md section 12.11 measured that matrix
and it is not clean. It holds the
log-normalised values put through normalisation and log1p a second time. The
damage is monotone within a cell and still sums to 1e4, so the docstring claim
"log1p CP10K" survived a year unchallenged - the only guard here checked for
non-finite values, which a doubly normalised matrix passes.

Section 12.12 audited every matrix in every input. Two are doubly normalised:
the `.raw` of `full_dataset.h5ad` and of `Neutrophils.h5ad`. The other eleven
`.raw` matrices are sound to 3e-4, which is float32. So the escape from a
damaged input is the per-cell-type files, not the pooled one, and every
population reproduces exactly from them - the responder and non-responder
counts of all thirteen match `full_dataset` cell for cell, which is what makes
this a change of source and not a change of cohort.

`assert_log1p_cp10k` now runs on every matrix before anything reads it, using
the same integer-ladder test as `00_Data_Audit/audit_raw_normalisation.py`.
A matrix that fails it stops that cell type with a message naming the file;
nothing is filled, substituted or skipped quietly.

Neutrophils have no sound source. Their `.raw` is damaged in both doors -
`Neutrophils.h5ad` and the neutrophil rows of `full_dataset.h5ad` - and the
sound neutrophil objects that exist elsewhere hold a different, earlier cell
set (51,616 cells against 61,167, and different R/NR counts), so they are not
a substitute. The cell type is kept in the run and fails loudly rather than
being dropped, because a missing row is easier to miss than a failing one.

The MAST model is unchanged, because the Methods describe it and the point of
this run is to repair the input, not to re-specify the test:

    zlm(~ condition + sample_id + cngeneson)

with the hurdle LRT on the condition contrast, FDR by Benjamini-Hochberg.

A Welch t-test on the same cells and the same gene set runs alongside it. That
is not redundancy: it is the answer to R1.8. Two methods that share nothing but
the input either agree, in which case the enrichment is not an artefact of
either, or they do not, in which case we would rather know. The concordance is
tabulated in deg_method_concordance.csv.

Inputs : submission-tree/01_Raw_Inputs/01_H5AD/<cell type>.h5ad     (.raw only), and
         upstream-pipeline/04_Final_Panels/00_Set_Ups/00_Data/01_Major_Cell_Types/
         {mast,plasma}_cell_integrated.h5ad for the two populations that have
         no file in 01_H5AD
Outputs: per cell type and timepoint,
           deg/<cell>_<phase>_mast.csv     MAST hurdle results
           deg/<cell>_<phase>_ttest.csv    Welch t-test on the same matrix
           gsea/<cell>_<phase>_<method>_hallmark.csv
         deg_method_concordance.csv, nfkb_by_celltype.csv, recompute_report.txt

The tables in outputs/ are NOT the output of this script as it now stands. They
were produced from full_dataset.h5ad, they are what verify_numbers.py checks the
manuscript against, and they were left in place deliberately. This script was run
on the sound inputs and the result is archived under
07_Archive/2026-08-31_deg_recompute_on_sound_per_cell_type_inputs/: it moves
fourteen checked numbers, including MoMac's pre-treatment NES from +1.06 to
-1.00, so adopting it is an author's decision and not a re-run. Running this
script overwrites outputs/ with the sound-input version and verify_numbers.py
will then fail until the manuscript is brought into line with it.

Run: python recompute_deg.py [--cell CELL] [--phase pre|post] [--jobs N]

MAST needs the environment's own libstdc++ ahead of the system one, or rpy2
cannot open libR.so and every MAST job is recorded as unavailable while the
t-test half still runs:

    LD_LIBRARY_PATH=$CONDA_PREFIX/lib conda run -n stad_ceacam python ...
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
# Hallmark gene sets come from the pinned file, not from Enrichr at run
# time. Passing the library name made the analysis depend on a remote
# service and recorded nothing about which release was used; the pinned
# copy is content-identical to the live one, so this changes no result.
# See 00_Reference/README.txt.
HALLMARK = str(HALLMARK_GMT)
NFKB_TERM = "TNF-alpha Signaling via NF-kB"

# Mast and plasma cells are the two populations with no file in 01_H5AD. These
# two upstream-pipeline objects hold them, their .raw carries the same 57,058 genes as
# the other eleven, and both pass audit_raw_normalisation.py. The sibling
# directory 01.1_Major_Cell_Types_Raw_Counts is NOT used: those files have no
# .raw at all, and the raw_counts layer beside them is broken - two nonzero
# genes and a library of three in a typical cell.
ROUND_4_MAJOR = (PROJECT_ROOT.parent / "upstream-pipeline" / "04_Final_Panels"
                 / "00_Set_Ups" / "00_Data" / "01_Major_Cell_Types")

# The thirteen populations the manuscript reports, and the object each is read
# from. Every file holds exactly one population, so there is no obs label to
# filter on: taking every cell in the file reproduces the responder and
# non-responder counts full_dataset.h5ad gives for that population, for all
# thirteen and at both timepoints. The files' own `major_cell_type` is not that
# partition - NK_cells, TCD4 and TCD8 all label their cells "T & NK cell", and
# MoMac.h5ad carries 487 cells labelled "DC" that full_dataset counts as
# Monocytes/Macrophages - so filtering on it would silently drop cells.
CELL_SOURCES = {
    "B_cells": CELL_TYPE_H5AD["B_cells"],
    "DC_cells": CELL_TYPE_H5AD["DC_cells"],
    "Endothelial_cells": CELL_TYPE_H5AD["Endothelial"],
    "Epithelial": CELL_TYPE_H5AD["Epithelial"],
    "Fibroblast": CELL_TYPE_H5AD["Fibroblast"],
    "Mast_cells": ROUND_4_MAJOR / "mast_cell_integrated.h5ad",
    "MoMac": CELL_TYPE_H5AD["MoMac"],
    # The rebuilt object: the shipped cell set and gene set, counts recovered
    # from the 70 per-sample aligner-era files on Ensembl ID, normalised once.
    # Built by 06_Clean_Data/rebuild_neutrophils_sound.py, cut to the deposit
    # obs by 06_Clean_Data/build_deposit_neutrophils_sound.py. Named through
    # paths.py rather than as a literal: a working-tree path
    # (06_Clean_Data/02_Rebuilt/) resolves to nothing in the deposit, so this
    # script would not run there at all.
    "Neutrophils": NEUTROPHILS_SOUND_H5AD,
    "NK_cells": CELL_TYPE_H5AD["NK_cells"],
    "Pericyte": CELL_TYPE_H5AD["Pericyte"],
    "Plasma_cells": ROUND_4_MAJOR / "plasma_cell_integrated.h5ad",
    "TCD4_cells": CELL_TYPE_H5AD["TCD4"],
    "TCD8_cells": CELL_TYPE_H5AD["TCD8"],
}

# Non-responder versus responder, at each timepoint. group1 is the reference,
# so a positive logFC means higher in non-responders - the direction the
# manuscript reports, and the one verified against six cytokines in
# FINDINGS.md section 7.
PHASES = {
    "pre": dict(column="stomach_pre_grouping", ref="Responsed", test="No-response"),
    "post": dict(column="stomach_post_grouping", ref="Responsed", test="No-response"),
}


LADDER_CELLS = 40
LADDER_TOL = 0.01       # the eleven sound matrices sit at 3e-4, which is float32


class BadMatrix(RuntimeError):
    """The input for one cell type is not what it claims to be.

    Not SystemExit: that is a BaseException, so it would walk past the per-cell
    handler in main() and take the other twelve populations down with it. This
    stops one cell type, records the reason in its status, and still makes the
    module exit non-zero.
    """


def assert_log1p_cp10k(X, source):
    """Refuse a matrix that is not log1p CP10K, by the integer-ladder test.

    This is the test of 00_Data_Audit/audit_raw_normalisation.py, brought
    inside the analysis. log1p CP10K means

        value = log1p(count * 1e4 / library)

    so within one cell `expm1(value) / expm1(smallest value)` has to come back
    1, 2, 3, 4 ... A matrix normalised a second time keeps the ranks and still
    sums to 1e4, which is why every check this script used to carry - monotone,
    finite, non-negative - passes it. Only the spacing gives it away.

    The scale factor is taken from the matrix rather than from
    obs["total_counts"]: Neutrophils.h5ad was gene-filtered after its library
    was recorded, so obs would fail it for the wrong reason.
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
            f"integers by up to {worst:.3f} (sound matrices sit at 3e-4). This "
            "is the signature of a second normalisation on already normalised "
            "values; see 00_Data_Audit/FINDINGS.md sections 12.11 and 12.12. "
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
    """Detected in at least MIN_PCT of the cells. Same rule as the original."""
    keep = (np.asarray((X > 0).mean(axis=0)).ravel() >= MIN_PCT)
    return keep


def run_mast(expr, groups, samples, genes, ref, test):
    """
    The original model, on a matrix that is what it claims to be.

    Any NaN or infinity aborts. The original replaced them with zero, which is
    how a matrix that was 31% NaN produced a table with no NaN in it.
    """
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter

    bad = ~np.isfinite(expr)
    if bad.any():
        raise SystemExit(
            f"{int(bad.sum())} non-finite values in the expression matrix. "
            "The input is not log1p CP10K; refusing to substitute zeros.")

    # Every string vector goes over as dtype=object. A numpy '<U' array is
    # converted by rpy2 into an R *array*, not a character vector, and MAST's
    # FromMatrix then rejects the object with "Order of `exprsArray` and
    # `cData` doesn't match" - the names are right, the type is not.
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

    # MAST is the sensitivity analysis; the Welch t-test is what the paper
    # reports. Run MAST first and an R that will not load - rpy2 unable to open
    # libR.so, which is the ordinary case without the LD_LIBRARY_PATH above -
    # kills the job before the primary analysis runs, while the module still
    # exits 0. MAST is isolated here, and its failure recorded rather than fatal.
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
    # Neutrophils by default, both phases: exactly what scripts/run.sh runs,
    # and exactly what outputs/ holds. This module exists for that one
    # population; the other twelve are module 12's, and their sound-input
    # tables are the archived sound recompute rather than anything this script
    # produces. With no arguments it used to start all thirteen, which needs
    # two upstream-pipeline objects that are not deposited - so the figure driver, which
    # launches every scripts/*.py with no arguments, killed it on a missing
    # file before it reached the population it is named after.
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
    # One summary per (cell, phase). run.sh invoked this once per phase, so
    # that is the name the shipped tables carry; writing it per row keeps a
    # two-phase run reproducing the same file names instead of inventing
    # recompute_summary_Neutrophils_None.csv.
    for r in rows:
        pd.DataFrame([r]).to_csv(
            OUT / f"recompute_summary_{r['cell']}_{r['phase']}.csv", index=False)
    print(df.to_string(index=False))

    # Exit non-zero when a job did not produce its primary result. Without this
    # the module reported success after all 26 jobs failed, and the driver
    # recorded it as reproduced. "no cells" and "too few cells" are documented
    # design decisions, not failures.
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
