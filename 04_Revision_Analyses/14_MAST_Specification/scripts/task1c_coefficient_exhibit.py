"""
Task 1, part 3: what the deposited fit actually returns for `condition`, and
what that number is a contrast of.

The design is rank deficient by one (task1_rank_and_fit.py). R resolves that
by pivoting: it scans the design's columns left to right and sets to NA the
first column that is linearly dependent on the ones before it. Which column
that is depends on the ORDER of the columns, and the order of the sample
dummies is the alphabetical order of the sample factor's levels - an arbitrary
labelling choice that carries no information about the experiment.

So the test is: relabel the sample factor's levels - same cells, same samples,
same partition, same design space, nothing about the data changed - and see
whether the `condition` coefficient moves. In an identifiable model it cannot.

A straight REVERSAL of the labels is not a valid probe: reversing maps the
last-dependent-column to the mirror position and can land on the same pair of
samples. Random permutations are used instead, and the first one is a
reversal only by accident.

Writes:
  outputs/task1_condition_coefficients.csv   per gene, per fit
  outputs/task1_glm_na_pattern.csv           which coefficient R pivots away
  outputs/task1_relabelling_spread.csv       the coefficient under K labellings
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from sources import PHASES, load_subset, filter_genes

OUT = Path(__file__).resolve().parents[1] / "outputs"
CELL = sys.argv[1] if len(sys.argv) > 1 else "MoMac"
PHASE = sys.argv[2] if len(sys.argv) > 2 else "post"
NG = int(sys.argv[3]) if len(sys.argv) > 3 else 200
NPERM = 6
cfg = PHASES[PHASE]

raw = load_subset(CELL, PHASE)
X = raw.X.toarray() if hasattr(raw.X, "toarray") else np.asarray(raw.X)
keep = filter_genes(X)
genes = np.asarray(raw.var_names)[keep][:NG]
expr = np.ascontiguousarray(X[:, keep][:, :NG], dtype=np.float64)
cng = np.asarray((X[:, keep] > 0).sum(axis=1), dtype=float)
g = raw.obs[cfg["column"]].astype(str).values
samples = raw.obs["sample"].astype(str).values
uniq = sorted(set(samples))
print(f"{CELL} {PHASE}: {expr.shape[0]} cells, {expr.shape[1]} genes examined, "
      f"{len(uniq)} samples", flush=True)

# Positional labels only; no specimen identifier is assigned into R or written.
idx = {s: i for i, s in enumerate(uniq)}
base = np.array([f"S{idx[s]+1:02d}" for s in samples])

rng = np.random.default_rng(7)
perm_labels = {}
for k in range(NPERM):
    order = rng.permutation(len(uniq))
    perm_labels[f"perm{k+1}"] = np.array([f"P{order[idx[s]]+1:02d}" for s in samples])
# condition of each positional sample, for reading the pivot off
cond_of = {f"S{idx[s]+1:02d}": str(c) for s, c in zip(samples, g)}

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
    ro.r.assign("comparison_group", str(cfg["test"]))
    ro.r.assign("sample_id", chars(base))
    for k, v in perm_labels.items():
        ro.r.assign(f"sid_{k}", chars(v))
    ro.r.assign("perm_names", chars(list(perm_labels)))

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
gm <- data.frame(primerid = as.character(gene_names), stringsAsFactors = FALSE)
rownames(gm) <- as.character(gene_names)
sca <- FromMatrix(exprsArray = expr_mat, cData = cm, fData = gm)
sca$condition <- relevel(factor(sca$condition), ref = as.character(reference_group))
ctr <- paste0("condition", as.character(comparison_group))

grab <- function(z, tag) {
  d <- data.frame(gene = rownames(z@coefC),
                  coefC = as.numeric(z@coefC[, ctr]),
                  coefD = as.numeric(z@coefD[, ctr]),
                  seC = as.numeric(sqrt(z@vcovC[ctr, ctr, ])),
                  convC = as.logical(z@converged[, "C"]),
                  convD = as.logical(z@converged[, "D"]),
                  stringsAsFactors = FALSE)
  names(d)[-1] <- paste0(names(d)[-1], "_", tag)
  d
}

z1 <- zlm(~ condition + sample_id + cngeneson, sca)
z3 <- zlm(~ condition + sample_id + cngeneson, sca, method = "glm", ebayes = FALSE)
z4 <- zlm(~ condition + cngeneson, sca)
z5 <- zlm(~ condition + cngeneson, sca, method = "glm", ebayes = FALSE)

# ---- which coefficients does the fit set to NA (i.e. pivot away)? ----
na_pattern <- rbind(
 data.frame(method = "bayesglm", coefficient = colnames(z1@coefC),
            n_NA_continuous = as.integer(colSums(is.na(z1@coefC))),
            n_NA_discrete   = as.integer(colSums(is.na(z1@coefD))),
            n_genes = nrow(z1@coefC), stringsAsFactors = FALSE),
 data.frame(method = "glm", coefficient = colnames(z3@coefC),
            n_NA_continuous = as.integer(colSums(is.na(z3@coefC))),
            n_NA_discrete   = as.integer(colSums(is.na(z3@coefD))),
            n_genes = nrow(z3@coefC), stringsAsFactors = FALSE))
cat("\\n--- which coefficients come back NA under the deposited formula ---\\n")
print(na_pattern[na_pattern$n_NA_continuous > 0 | na_pattern$n_NA_discrete > 0, ])
cat("(a coefficient with n_NA equal to n_genes is the column the fit pivoted away)\\n")

out <- Reduce(function(a, b) merge(a, b, by = "gene"),
              list(grab(z1, "dep_bayesglm"), grab(z3, "dep_glm"),
                   grab(z4, "A_bayesglm"), grab(z5, "A_glm")))

# ---- relabelling spread ----
spread <- grab(z1, "base")[, c("gene", "coefC_base")]
pivoted <- data.frame(labelling = "base",
                      pivoted_away = paste(colnames(z1@coefC)[colSums(is.na(z1@coefC)) == nrow(z1@coefC)],
                                           collapse = ";"), stringsAsFactors = FALSE)
for (nm in perm_names) {
  s2 <- sca
  colData(s2)$sample_id <- as.factor(as.character(get(paste0("sid_", nm))))
  zz <- zlm(~ condition + sample_id + cngeneson, s2)
  d <- data.frame(gene = rownames(zz@coefC), v = as.numeric(zz@coefC[, ctr]),
                  stringsAsFactors = FALSE)
  names(d)[2] <- paste0("coefC_", nm)
  spread <- merge(spread, d, by = "gene")
  pivoted <- rbind(pivoted, data.frame(labelling = nm,
      pivoted_away = paste(colnames(zz@coefC)[colSums(is.na(zz@coefC)) == nrow(zz@coefC)],
                           collapse = ";"), stringsAsFactors = FALSE))
}
cat("\\n--- column the fit pivoted away, per labelling ---\\n")
print(pivoted)
''')
with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
    tab = pd.DataFrame(ro.conversion.rpy2py(ro.r["out"]))
    nap = pd.DataFrame(ro.conversion.rpy2py(ro.r["na_pattern"]))
    spr = pd.DataFrame(ro.conversion.rpy2py(ro.r["spread"]))
    piv = pd.DataFrame(ro.conversion.rpy2py(ro.r["pivoted"]))

tab.insert(0, "phase", PHASE); tab.insert(0, "cell", CELL)
tab.to_csv(OUT / "task1_condition_coefficients.csv", index=False)
nap.to_csv(OUT / "task1_glm_na_pattern.csv", index=False)

cols = [c for c in spr.columns if c.startswith("coefC_")]
spr["range_across_labellings"] = spr[cols].max(axis=1) - spr[cols].min(axis=1)
spr["sign_flips"] = spr[cols].apply(
    lambda r: int(len(set(np.sign(r.values))) > 1), axis=1)
spr.insert(0, "phase", PHASE); spr.insert(0, "cell", CELL)
spr.to_csv(OUT / "task1_relabelling_spread.csv", index=False)
piv.to_csv(OUT / "task1_pivoted_column_per_labelling.csv", index=False)

print("\n--- the condition coefficient across", len(cols), "arbitrary labellings "
      "of the same samples ---")
print(f"genes examined                    : {len(spr)}")
print(f"median range across labellings    : {spr['range_across_labellings'].median():.4f}")
print(f"90th percentile of that range     : {spr['range_across_labellings'].quantile(.9):.4f}")
print(f"max range                         : {spr['range_across_labellings'].max():.4f}")
print(f"genes whose coefficient changes SIGN with the labelling: "
      f"{int(spr['sign_flips'].sum())} of {len(spr)} "
      f"({100*spr['sign_flips'].mean():.1f}%)")

from scipy import stats
recs = []
for a, b in [("dep_bayesglm", "dep_glm"), ("dep_bayesglm", "A_bayesglm"),
             ("A_bayesglm", "A_glm"), ("dep_glm", "A_bayesglm")]:
    va, vb = tab[f"coefC_{a}"], tab[f"coefC_{b}"]
    m = va.notna() & vb.notna()
    recs.append(dict(fit_a=a, fit_b=b, n=int(m.sum()),
                     pearson=round(float(stats.pearsonr(va[m], vb[m]).statistic), 4),
                     spearman=round(float(stats.spearmanr(va[m], vb[m]).statistic), 4),
                     median_abs_diff=round(float((va[m] - vb[m]).abs().median()), 4),
                     max_abs_diff=round(float((va[m] - vb[m]).abs().max()), 4),
                     pct_same_sign=round(float(100 * (np.sign(va[m]) == np.sign(vb[m])).mean()), 1)))
r = pd.DataFrame(recs)
r.to_csv(OUT / "task1_coefficient_agreement.csv", index=False)
print()
print(r.to_string(index=False))
print("\ndone", flush=True)
