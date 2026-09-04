"""
Task 1, the mechanism: MAST does not report the aliasing. It silently removes
one column.

model.matrix(~ condition + sample_id + cngeneson) has 13 columns and rank 12.
This prints the columns MAST's fit actually kept, for the base labelling and
for three random relabellings of the same samples. If the kept set differs
between labellings, the `condition` coefficient is a different contrast each
time - which is what the 95.5% sign flips in task1_relabelling_spread.csv are.

Small contrast, few genes: this is a question about the design, not the data.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from sources import PHASES, load_subset, filter_genes

OUT = Path(__file__).resolve().parents[1] / "outputs"
CELL = sys.argv[1] if len(sys.argv) > 1 else "Mast_cells"
PHASE = sys.argv[2] if len(sys.argv) > 2 else "pre"
cfg = PHASES[PHASE]

raw = load_subset(CELL, PHASE)
X = raw.X.toarray() if hasattr(raw.X, "toarray") else np.asarray(raw.X)
keep = filter_genes(X)
genes = np.asarray(raw.var_names)[keep][:10]
expr = np.ascontiguousarray(X[:, keep][:, :10], dtype=np.float64)
cng = np.asarray((X[:, keep] > 0).sum(axis=1), dtype=float)
g = raw.obs[cfg["column"]].astype(str).values
samples = raw.obs["sample"].astype(str).values
uniq = sorted(set(samples))
idx = {s: i for i, s in enumerate(uniq)}
base = np.array([f"S{idx[s]+1:02d}" for s in samples])
rng = np.random.default_rng(11)
perms = {}
for k in range(3):
    o = rng.permutation(len(uniq))
    perms[f"perm{k+1}"] = np.array([f"P{o[idx[s]]+1:02d}" for s in samples])
# condition of each base label, so the report can say whether the dropped
# sample is an R or an NR one. No specimen identifier is used or written.
cond_by_label = {}
for lab, c in zip(base, g):
    cond_by_label[lab] = "R" if c == cfg["ref"] else "NR"
print(f"{CELL} {PHASE}: {expr.shape[0]} cells, {len(uniq)} samples, "
      f"R/NR by positional label: "
      f"{ {k: cond_by_label[k] for k in sorted(cond_by_label)} }", flush=True)

import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter


def chars(a):
    return np.asarray([str(x) for x in a], dtype=object)


with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
    ro.r.assign("expr_matrix", expr)
    ro.r.assign("cngeneson_src", cng)
    ro.r.assign("cell_metadata", chars(g))
    ro.r.assign("cell_names", chars([f"c{i}" for i in range(len(g))]))
    ro.r.assign("gene_names", chars(genes))
    ro.r.assign("reference_group", str(cfg["ref"]))
    ro.r.assign("sample_id", chars(base))
    for k, v in perms.items():
        ro.r.assign(f"sid_{k}", chars(v))
    ro.r.assign("perm_names", chars(list(perms)))

ro.r('''
suppressPackageStartupMessages({library(MAST); library(SingleCellExperiment)})
expr_mat <- t(expr_matrix)
rownames(expr_mat) <- as.character(gene_names)
colnames(expr_mat) <- as.character(cell_names)
cm <- data.frame(condition = as.character(cell_metadata),
                 sample_id = as.factor(as.character(sample_id)),
                 wellKey = as.character(cell_names), stringsAsFactors = FALSE)
rownames(cm) <- as.character(cell_names)
cm$cngeneson <- as.numeric(scale(cngeneson_src))
cm$condition <- relevel(factor(cm$condition), ref = as.character(reference_group))
gm <- data.frame(primerid = as.character(gene_names), stringsAsFactors = FALSE)
rownames(gm) <- as.character(gene_names)
sca <- FromMatrix(exprsArray = expr_mat, cData = cm, fData = gm)
sca$condition <- relevel(factor(sca$condition), ref = as.character(reference_group))

mm <- model.matrix(~ condition + sample_id + cngeneson, cm)
cat("\\ndesign built by model.matrix : ", ncol(mm), "columns, rank", qr(mm)$rank, "\\n")
cat("its columns                  : ", paste(colnames(mm), collapse=", "), "\\n\\n")

rows <- list()
z <- zlm(~ condition + sample_id + cngeneson, sca)
kept <- colnames(z@coefC)
cat("MAST kept                    : ", ncol(z@coefC), "columns\\n")
cat("dropped by MAST, base labels : ",
    paste(setdiff(colnames(mm), kept), collapse=", "), "\\n")
rows[["base"]] <- paste(setdiff(colnames(mm), kept), collapse=", ")

for (nm in perm_names) {
  s2 <- sca
  colData(s2)$sample_id <- as.factor(as.character(get(paste0("sid_", nm))))
  cm2 <- cm; cm2$sample_id <- as.factor(as.character(get(paste0("sid_", nm))))
  mm2 <- model.matrix(~ condition + sample_id + cngeneson, cm2)
  zz <- zlm(~ condition + sample_id + cngeneson, s2)
  d <- setdiff(colnames(mm2), colnames(zz@coefC))
  cat("dropped by MAST,", nm, "      : ", paste(d, collapse=", "), "\\n")
  rows[[nm]] <- paste(d, collapse=", ")
}
dropped <- data.frame(labelling = names(rows), dropped_column = unlist(rows),
                      stringsAsFactors = FALSE)
''')
with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
    dr = pd.DataFrame(ro.conversion.rpy2py(ro.r["dropped"]))
dr.insert(0, "phase", PHASE); dr.insert(0, "cell", CELL)
dr["dropped_sample_is"] = dr["dropped_column"].map(
    lambda s: cond_by_label.get(s.replace("sample_id", ""), "") if s.startswith("sample_id") else "")
dr.to_csv(OUT / "task1_column_dropped_by_MAST.csv", index=False)
print()
print(dr.to_string(index=False))
print("\ndone", flush=True)
