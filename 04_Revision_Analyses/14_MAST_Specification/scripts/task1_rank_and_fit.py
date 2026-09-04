"""
Task 1, parts 2 and 3.

2. Rank of the deposited design matrix `~ condition + sample_id + cngeneson`
   against its column count, and which columns are aliased.
3. What MAST actually returns for the `condition` coefficient under that
   design: coefficients, standard errors, convergence.

Plus one test the brief does not ask for and which settles the question
without argument: refit with the *same* design and the same data, changing
only which sample_id level R treats as the reference. In an identifiable
model the `condition` coefficient cannot depend on that. If it moves, the
coefficient is not an estimate of anything.

Representative contrast: MoMac post, the contrast the brief singles out
(+2.228 by t-test, +0.892 by MAST).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from sources import CELL_SOURCES, PHASES, load_subset, filter_genes

OUT = Path(__file__).resolve().parents[1] / "outputs"
CELL, PHASE = (sys.argv[1] if len(sys.argv) > 1 else "MoMac",
               sys.argv[2] if len(sys.argv) > 2 else "post")
cfg = PHASES[PHASE]

print(f"### representative contrast: {CELL} {PHASE}", flush=True)
raw = load_subset(CELL, PHASE)
X = raw.X.toarray() if hasattr(raw.X, "toarray") else np.asarray(raw.X)
keep = filter_genes(X)
genes = np.asarray(raw.var_names)[keep]
expr = np.ascontiguousarray(X[:, keep], dtype=np.float64)
g = raw.obs[cfg["column"]].astype(str).values
samples = raw.obs["sample"].astype(str).values
print(f"cells {expr.shape[0]}  genes tested {expr.shape[1]}  "
      f"samples {len(set(samples))}", flush=True)

import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter


def chars(a):
    return np.asarray([str(x) for x in a], dtype=object)


# Sample identifiers never leave this process: they are replaced by a
# positional label before anything is assigned into R or printed.
uniq = sorted(set(samples))
anon = {s: f"S{i+1:02d}" for i, s in enumerate(uniq)}
samples_anon = np.array([anon[s] for s in samples])
# The reverse-order relabelling used for the invariance test. Same partition
# of cells, different alphabetical order, so factor() picks a different
# reference level.
rev = {s: f"S{len(uniq)-i:02d}z" for i, s in enumerate(uniq)}
samples_rev = np.array([rev[s] for s in samples])

NG = 40   # genes for the fit inspection; the rank question needs none

with localconverter(ro.default_converter + pandas2ri.converter
                    + numpy2ri.converter):
    ro.r.assign("expr_matrix", expr[:, :NG])
    ro.r.assign("cngeneson_src", np.asarray((expr > 0).sum(axis=1), dtype=float))
    ro.r.assign("cell_metadata", chars(g))
    ro.r.assign("cell_names", chars([f"c{i}" for i in range(len(g))]))
    ro.r.assign("gene_names", chars(genes[:NG]))
    ro.r.assign("reference_group", str(cfg["ref"]))
    ro.r.assign("comparison_group", str(cfg["test"]))
    ro.r.assign("sample_id", chars(samples_anon))
    ro.r.assign("sample_id_rev", chars(samples_rev))

ro.r('''
suppressPackageStartupMessages({library(MAST); library(SingleCellExperiment)})
expr_mat <- t(expr_matrix)
rownames(expr_mat) <- as.character(gene_names)
colnames(expr_mat) <- as.character(cell_names)
cell_meta <- data.frame(condition = as.character(cell_metadata),
                        sample_id = as.factor(as.character(sample_id)),
                        sample_id_rev = as.factor(as.character(sample_id_rev)),
                        wellKey = as.character(cell_names),
                        stringsAsFactors = FALSE)
rownames(cell_meta) <- as.character(cell_names)
# cngeneson from the FULL gene set, as recompute_deg.py computes it
cell_meta$cngeneson <- as.numeric(scale(cngeneson_src))
cell_meta$condition <- relevel(factor(cell_meta$condition),
                               ref = as.character(reference_group))

cat("\\n=== 2. RANK OF THE DEPOSITED DESIGN ===\\n")
mm <- model.matrix(~ condition + sample_id + cngeneson, cell_meta)
cat("columns          :", ncol(mm), "\\n")
cat("qr()$rank        :", qr(mm)$rank, "\\n")
cat("rank deficiency  :", ncol(mm) - qr(mm)$rank, "\\n")
cat("column names     :", paste(colnames(mm), collapse=", "), "\\n")

cat("\\n--- alias() on the equivalent lm ---\\n")
y <- as.numeric(expr_mat[1, ])
al <- alias(lm(y ~ condition + sample_id + cngeneson, data = cell_meta))
print(al$Complete)

cat("\\n--- for comparison, the two candidate correct designs ---\\n")
mmA <- model.matrix(~ condition + cngeneson, cell_meta)
cat("A ~condition+cngeneson        cols", ncol(mmA), " rank", qr(mmA)$rank, "\\n")
mmS <- model.matrix(~ sample_id + cngeneson, cell_meta)
cat("  ~sample_id+cngeneson        cols", ncol(mmS), " rank", qr(mmS)$rank, "\\n")

cat("\\n=== 3. WHAT MAST RETURNS UNDER THE ALIASED DESIGN ===\\n")
gene_meta <- data.frame(primerid = as.character(gene_names), stringsAsFactors = FALSE)
rownames(gene_meta) <- as.character(gene_names)
sca <- FromMatrix(exprsArray = expr_mat, cData = cell_meta, fData = gene_meta)
sca$condition <- relevel(factor(sca$condition), ref = as.character(reference_group))
contrast_name <- paste0("condition", as.character(comparison_group))

cat("\\n-- deposited call: zlm(~condition+sample_id+cngeneson, sca)  [default method=bayesglm]\\n")
zAlias <- zlm(~ condition + sample_id + cngeneson, sca)
cat("converged C (discrete/continuous slots):",
    sum(zAlias@converged[, "C"]), "of", nrow(zAlias@converged), "\\n")
cat("converged D:", sum(zAlias@converged[, "D"]), "of", nrow(zAlias@converged), "\\n")
cC <- zAlias@coefC[, contrast_name]
cD <- zAlias@coefD[, contrast_name]
cat("NA in coefC[condition]:", sum(is.na(cC)), " NA in coefD[condition]:", sum(is.na(cD)), "\\n")
cat("coefC[condition] range: ", paste(round(range(cC, na.rm=TRUE), 4), collapse=" .. "), "\\n")
cat("coefD[condition] range: ", paste(round(range(cD, na.rm=TRUE), 4), collapse=" .. "), "\\n")
vC <- sqrt(zAlias@vcovC[contrast_name, contrast_name, ])
cat("se(coefC[condition]) range:", paste(signif(range(vC, na.rm=TRUE), 4), collapse=" .. "),
    " zero SEs:", sum(vC == 0, na.rm=TRUE), "\\n")

cat("\\n-- same design, method='glm' (no prior; R drops aliased columns)\\n")
zGlm <- zlm(~ condition + sample_id + cngeneson, sca, method = "glm", ebayes = FALSE)
gC <- zGlm@coefC[, contrast_name]
gD <- zGlm@coefD[, contrast_name]
cat("NA in coefC[condition]:", sum(is.na(gC)), "of", length(gC), "\\n")
cat("NA in coefD[condition]:", sum(is.na(gD)), "of", length(gD), "\\n")

cat("\\n-- INVARIANCE TEST: identical data and design, different sample_id reference level\\n")
sca2 <- sca
colData(sca2)$sample_id <- colData(sca2)$sample_id_rev
zAlias2 <- zlm(~ condition + sample_id + cngeneson, sca2)
c2 <- zAlias2@coefC[, contrast_name]
cat("max |coefC[condition] - coefC_relabelled[condition]| :",
    signif(max(abs(cC - c2), na.rm = TRUE), 5), "\\n")
cat("median |difference| :", signif(median(abs(cC - c2), na.rm = TRUE), 5), "\\n")
cat("Spearman rho between the two :",
    signif(cor(cC, c2, method = "spearman", use = "complete.obs"), 4), "\\n")

cat("\\n-- the same invariance test on spec A ~condition+cngeneson (full rank)\\n")
zA <- zlm(~ condition + cngeneson, sca)
aC <- zA@coefC[, contrast_name]
cat("coefC[condition] range:", paste(round(range(aC, na.rm=TRUE), 4), collapse=" .. "), "\\n")
cat("median |coefC_aliased - coefC_specA| :",
    signif(median(abs(cC - aC), na.rm = TRUE), 5), "\\n")

cat("\\n-- LRT the deposited script extracts, on these", length(gene_names), "genes\\n")
s <- summary(zAlias, doLRT = contrast_name)$datatable
p <- s[contrast == contrast_name & component == "H", ][["Pr(>Chisq)"]]
cat("hurdle P: min", signif(min(p, na.rm=TRUE), 4), " median", signif(median(p, na.rm=TRUE), 4),
    " n NA", sum(is.na(p)), "\\n")
lf <- s[contrast == contrast_name & component == "logFC", ][["coef"]]
cat("logFC   : range", paste(round(range(lf, na.rm=TRUE), 4), collapse=" .. "),
    " n NA", sum(is.na(lf)), "\\n")
''')
print("\ndone", flush=True)
