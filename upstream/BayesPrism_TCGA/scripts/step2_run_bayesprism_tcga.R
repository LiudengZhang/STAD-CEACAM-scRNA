#!/usr/bin/env Rscript
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# Step 2: Run BayesPrism deconvolution on TCGA-STAD bulk data
#
# Input:  sc_counts.tsv, sc_cell_types.tsv, tcga_bulk_counts.tsv
# Output: tcga_bayesprism_fractions.tsv
#         tcga_bayesprism_fibroblast_expression.tsv
#         tcga_bayesprism_momac_expression.tsv
#         tcga_bayesprism_result.rds

suppressPackageStartupMessages({
    library(BayesPrism)
    library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
    WORK_DIR <- dirname(script_path)
} else {
    WORK_DIR <- getwd()
}
setwd(WORK_DIR)

sc_dt <- fread("sc_counts.tsv", sep="\t", header=TRUE)
sc_rownames <- sc_dt[[1]]
sc_dt[, 1 := NULL]
sc_mat <- as.matrix(sc_dt)
rownames(sc_mat) <- sc_rownames
rm(sc_dt); gc()

sc_labels <- fread("sc_cell_types.tsv", sep="\t", header=TRUE)

bk_dt <- fread("tcga_bulk_counts.tsv", sep="\t", header=TRUE)
bk_rownames <- bk_dt[[1]]
bk_dt[, 1 := NULL]
bk_mat <- as.matrix(bk_dt)
rownames(bk_mat) <- bk_rownames
rm(bk_dt); gc()

common_genes <- intersect(colnames(sc_mat), colnames(bk_mat))
sc_mat <- sc_mat[, common_genes]
bk_mat <- bk_mat[, common_genes]

# Cell type labels (use cell_type for both to avoid sparse state issue)
cell_type_labels <- sc_labels$cell_type
cell_state_labels <- sc_labels$cell_type

myPrism <- new.prism(
    reference = sc_mat,
    mixture = bk_mat,
    input.type = "GEP",
    cell.type.labels = cell_type_labels,
    cell.state.labels = cell_state_labels,
    key = NULL,
    outlier.cut = 0.01,
    outlier.fraction = 0.1
)

rm(sc_mat, bk_mat); gc()

bp_result <- run.prism(prism = myPrism, n.cores = 4)

theta <- get.fraction(bp = bp_result,
                      which.theta = "final",
                      state.or.type = "type")
write.table(theta, "tcga_bayesprism_fractions.tsv", sep="\t", quote=FALSE)

for (ct in colnames(theta)) {
    cat(sprintf("  %s: %.4f (%.1f%%)\n", ct, mean(theta[, ct]), mean(theta[, ct]) * 100))
}

fib_expr <- get.exp(bp = bp_result,
                    state.or.type = "type",
                    cell.name = "Fibroblast")
write.table(fib_expr, "tcga_bayesprism_fibroblast_expression.tsv", sep="\t", quote=FALSE)

momac_expr <- get.exp(bp = bp_result,
                      state.or.type = "type",
                      cell.name = "Monocytes/Macrophages")
write.table(momac_expr, "tcga_bayesprism_momac_expression.tsv", sep="\t", quote=FALSE)

key_genes <- c("IL1B", "IL6", "OSM", "TNF", "NFKB1", "CEACAM5", "CEACAM6")
for (g in key_genes) {
    if (g %in% colnames(fib_expr)) {
        cat(sprintf("  %s: mean=%.4f, median=%.4f\n", g, mean(fib_expr[, g]), median(fib_expr[, g])))
    }
}
for (g in key_genes) {
    if (g %in% colnames(momac_expr)) {
        cat(sprintf("  %s: mean=%.4f, median=%.4f\n", g, mean(momac_expr[, g]), median(momac_expr[, g])))
    }
}

saveRDS(bp_result, "tcga_bayesprism_result.rds")
