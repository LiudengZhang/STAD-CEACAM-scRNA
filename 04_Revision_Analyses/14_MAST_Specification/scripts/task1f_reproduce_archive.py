"""
The closing test: the deposited MAST result is not stable under renaming the
samples, and that is demonstrable against the deposit's own archived output.

For one pre contrast, the deposited model is fitted twice on identical data.
The two runs differ in one thing: the order in which the sample factor's
levels fall.

  py_order : positional labels assigned in Python's sorted() order, which is
             what this module's harness used
  r_order  : positional labels assigned in R's factor() level order, which is
             the order the deposited script's real identifiers fall into

Each is compared, gene by gene, against
07_Archive/2026-08-31_deg_recompute_on_sound_per_cell_type_inputs/.../deg/.
If the archive is reproduced under one ordering and not the other, then the
number the deposit reports is a function of how its sample names happen to
sort - which is the definition of a coefficient that is not identified.

No specimen identifier is read into a variable that is written or printed.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from sources import PHASES, load_subset, filter_genes

MOD = Path(__file__).resolve().parents[1]
OUT = MOD / "outputs"
ARCH = (MOD.parents[1] / "07_Archive"
        / "2026-08-31_deg_recompute_on_sound_per_cell_type_inputs"
        / "04_Revision_Analyses" / "12_R1.8_DEG_Recompute" / "outputs" / "deg")
CELL = sys.argv[1] if len(sys.argv) > 1 else "Endothelial_cells"
PHASE = sys.argv[2] if len(sys.argv) > 2 else "pre"
cfg = PHASES[PHASE]

raw = load_subset(CELL, PHASE)
X = raw.X.toarray() if hasattr(raw.X, "toarray") else np.asarray(raw.X)
keep = filter_genes(X)
genes = np.asarray(raw.var_names)[keep]
expr = np.ascontiguousarray(X[:, keep], dtype=np.float64)
g = raw.obs[cfg["column"]].astype(str).values
samples = raw.obs["sample"].astype(str).values
print(f"{CELL} {PHASE}: {expr.shape[0]} cells, {expr.shape[1]} genes, "
      f"{len(set(samples))} samples", flush=True)

import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter


def chars(a):
    return np.asarray([str(x) for x in a], dtype=object)


uniq = list(pd.unique(samples))
py_sorted = sorted(uniq)
with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
    ro.r.assign("s_tmp", chars(uniq))
r_sorted = list(ro.r('levels(factor(as.character(s_tmp)))'))

labellings = {}
for tag, order in (("py_order", py_sorted), ("r_order", r_sorted)):
    pos = {s: i for i, s in enumerate(order)}
    labellings[tag] = np.array([f"S{pos[s]+1:02d}" for s in samples])

arch = pd.read_csv(ARCH / f"{CELL}_{PHASE}_mast.csv").set_index("gene")
rows = []
for tag, lab in labellings.items():
    with localconverter(ro.default_converter + pandas2ri.converter
                        + numpy2ri.converter):
        ro.r.assign("expr_matrix", expr)
        ro.r.assign("cell_metadata", chars(g))
        ro.r.assign("cell_names", chars([f"c{i}" for i in range(len(g))]))
        ro.r.assign("gene_names", chars(genes))
        ro.r.assign("reference_group", str(cfg["ref"]))
        ro.r.assign("comparison_group", str(cfg["test"]))
        ro.r.assign("sample_id", chars(lab))
    ro.r('''
    suppressPackageStartupMessages({library(MAST); library(SingleCellExperiment)})
    expr_mat <- t(expr_matrix)
    rownames(expr_mat) <- as.character(gene_names)
    colnames(expr_mat) <- as.character(cell_names)
    cm <- data.frame(condition = as.character(cell_metadata),
                     sample_id = as.factor(as.character(sample_id)),
                     wellKey = as.character(cell_names), stringsAsFactors = FALSE)
    rownames(cm) <- as.character(cell_names)
    cm$cngeneson <- as.numeric(scale(colSums(expr_mat > 0)))
    gm <- data.frame(primerid = as.character(gene_names), stringsAsFactors = FALSE)
    rownames(gm) <- as.character(gene_names)
    sca <- FromMatrix(exprsArray = expr_mat, cData = cm, fData = gm)
    sca$condition <- relevel(factor(sca$condition), ref = as.character(reference_group))
    ctr <- paste0("condition", as.character(comparison_group))
    mm <- model.matrix(~ condition + sample_id + cngeneson, cm)
    z <- zlm(~ condition + sample_id + cngeneson, sca)
    dropped <- paste(setdiff(colnames(mm), colnames(z@coefC)), collapse=", ")
    s <- summary(z, doLRT = ctr)$datatable
    fc <- merge(s[contrast==ctr & component=="H", .(primerid, `Pr(>Chisq)`)],
                s[contrast==ctr & component=="logFC", .(primerid, coef)], by="primerid")
    results <- data.frame(gene = fc$primerid, logfoldchanges = fc$coef,
                          pvals = fc$`Pr(>Chisq)`, stringsAsFactors = FALSE)
    results <- results[!is.na(results$pvals) & !is.na(results$logfoldchanges), ]
    ''')
    with localconverter(ro.default_converter + pandas2ri.converter
                        + numpy2ri.converter):
        res = pd.DataFrame(ro.conversion.rpy2py(ro.r["results"])).set_index("gene")
        dropped = str(ro.r["dropped"][0])
    c = arch.index.intersection(res.index)
    dlf = (arch.loc[c, "logfoldchanges"] - res.loc[c, "logfoldchanges"]).abs()
    dlp = (np.log10(arch.loc[c, "pvals"].clip(lower=1e-300))
           - np.log10(res.loc[c, "pvals"].clip(lower=1e-300))).abs()
    rows.append(dict(cell=CELL, phase=PHASE, labelling=tag,
                     dropped_column=dropped, n_genes=len(c),
                     max_abs_dlogFC_vs_archive=float(dlf.max()),
                     median_abs_dlogFC_vs_archive=float(dlf.median()),
                     max_abs_dlog10P_vs_archive=float(dlp.max()),
                     reproduces_archive=bool(dlf.max() < 1e-9)))
    print(f"  [{tag}] MAST dropped {dropped!r}; vs archive "
          f"max|dlogFC|={dlf.max():.3g} max|dlog10P|={dlp.max():.3g} "
          f"reproduces={dlf.max() < 1e-9}", flush=True)

df = pd.DataFrame(rows)
p = OUT / "task1_archive_reproduction_by_labelling.csv"
df.to_csv(p, mode="a", header=not p.exists(), index=False)
print()
print(df.to_string(index=False))
print("\ndone", flush=True)
