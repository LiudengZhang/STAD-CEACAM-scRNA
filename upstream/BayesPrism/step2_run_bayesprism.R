#!/usr/bin/env Rscript
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# Step 2: Run BayesPrism deconvolution on TIGER bulk data
#
# Input:  sc_counts.tsv, sc_cell_types.tsv, bulk_counts.tsv
# Output: bayesprism_fractions.tsv, bayesprism_epithelial_expression.tsv

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
bk_dt <- fread("bulk_counts.tsv", sep="\t", header=TRUE)
bk_rownames <- bk_dt[[1]]
bk_dt[, 1 := NULL]
bk_mat <- as.matrix(bk_dt)
rownames(bk_mat) <- bk_rownames
rm(bk_dt); gc()

# Cell type and state labels
# Use cell_type for both to avoid sparse state issue (some states have <10 cells)
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

# Free memory
rm(sc_mat, bk_mat); gc()

bp_result <- run.prism(prism = myPrism, n.cores = 4)

theta <- get.fraction(bp = bp_result,
                      which.theta = "final",
                      state.or.type = "type")
write.table(theta, "bayesprism_fractions.tsv", sep="\t", quote=FALSE)

epi_expr <- get.exp(bp = bp_result,
                    state.or.type = "type",
                    cell.name = "Epithelial cells")

write.table(epi_expr, "bayesprism_epithelial_expression.tsv", sep="\t", quote=FALSE)

ceacam_genes <- grep("^CEACAM[0-9]", colnames(epi_expr), value=TRUE)

if (length(ceacam_genes) > 0) {
    ceacam_df <- epi_expr[, ceacam_genes, drop=FALSE]
    write.table(ceacam_df, "bayesprism_epithelial_ceacam.tsv", sep="\t", quote=FALSE)
    cat("Saved: bayesprism_epithelial_ceacam.tsv\n")

    cat("\nCEACAM expression summary (epithelial-specific):\n")
    for (g in ceacam_genes) {
        cat(sprintf("  %s: mean=%.4f, median=%.4f\n", g, mean(epi_expr[, g]), median(epi_expr[, g])))
    }
}

saveRDS(bp_result, "bayesprism_result.rds")
