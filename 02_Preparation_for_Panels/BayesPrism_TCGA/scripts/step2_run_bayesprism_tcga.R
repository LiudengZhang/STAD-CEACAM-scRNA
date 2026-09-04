#!/usr/bin/env Rscript
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

cat("============================================================\n")
cat("Step 2: BayesPrism Deconvolution — TCGA-STAD\n")
cat("============================================================\n")

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
    WORK_DIR <- dirname(script_path)
} else {
    WORK_DIR <- getwd()
}
setwd(WORK_DIR)
cat(sprintf("Working directory: %s\n", WORK_DIR))

# --- Load scRNA reference ---
cat("\nLoading scRNA-seq reference...\n")
sc_dt <- fread("sc_counts.tsv", sep="\t", header=TRUE)
sc_rownames <- sc_dt[[1]]
sc_dt[, 1 := NULL]
sc_mat <- as.matrix(sc_dt)
rownames(sc_mat) <- sc_rownames
rm(sc_dt); gc()
cat(sprintf("  scRNA: %d cells x %d genes\n", nrow(sc_mat), ncol(sc_mat)))

sc_labels <- fread("sc_cell_types.tsv", sep="\t", header=TRUE)
cat(sprintf("  Cell types: %s\n", paste(unique(sc_labels$cell_type), collapse=", ")))

# --- Load TCGA bulk ---
cat("\nLoading TCGA-STAD bulk data...\n")
bk_dt <- fread("tcga_bulk_counts.tsv", sep="\t", header=TRUE)
bk_rownames <- bk_dt[[1]]
bk_dt[, 1 := NULL]
bk_mat <- as.matrix(bk_dt)
rownames(bk_mat) <- bk_rownames
rm(bk_dt); gc()
cat(sprintf("  Bulk: %d samples x %d genes\n", nrow(bk_mat), ncol(bk_mat)))

# --- Ensure gene overlap ---
common_genes <- intersect(colnames(sc_mat), colnames(bk_mat))
cat(sprintf("  Common genes: %d\n", length(common_genes)))
sc_mat <- sc_mat[, common_genes]
bk_mat <- bk_mat[, common_genes]

# Cell type labels (use cell_type for both to avoid sparse state issue)
cell_type_labels <- sc_labels$cell_type
cell_state_labels <- sc_labels$cell_type

# --- Create prism object ---
cat("\nCreating BayesPrism object...\n")
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

# --- Run BayesPrism ---
cat("\nRunning BayesPrism (410 samples, this will take a while)...\n")
cat(sprintf("Start time: %s\n", Sys.time()))
bp_result <- run.prism(prism = myPrism, n.cores = 4)
cat(sprintf("End time: %s\n", Sys.time()))

# --- Extract results ---
cat("\nExtracting results...\n")

# 1. Cell type fractions
theta <- get.fraction(bp = bp_result,
                      which.theta = "final",
                      state.or.type = "type")
cat(sprintf("Fractions: %d samples x %d types\n", nrow(theta), ncol(theta)))
write.table(theta, "tcga_bayesprism_fractions.tsv", sep="\t", quote=FALSE)
cat("Saved: tcga_bayesprism_fractions.tsv\n")

# Print mean fractions
cat("\nMean cell type fractions:\n")
for (ct in colnames(theta)) {
    cat(sprintf("  %s: %.4f (%.1f%%)\n", ct, mean(theta[, ct]), mean(theta[, ct]) * 100))
}

# 2. Fibroblast-specific expression
cat("\nExtracting Fibroblast-specific expression...\n")
fib_expr <- get.exp(bp = bp_result,
                    state.or.type = "type",
                    cell.name = "Fibroblast")
cat(sprintf("Fibroblast expression: %d samples x %d genes\n", nrow(fib_expr), ncol(fib_expr)))
write.table(fib_expr, "tcga_bayesprism_fibroblast_expression.tsv", sep="\t", quote=FALSE)
cat("Saved: tcga_bayesprism_fibroblast_expression.tsv\n")

# 3. Monocytes/Macrophages-specific expression
cat("\nExtracting Monocytes/Macrophages-specific expression...\n")
momac_expr <- get.exp(bp = bp_result,
                      state.or.type = "type",
                      cell.name = "Monocytes/Macrophages")
cat(sprintf("MoMac expression: %d samples x %d genes\n", nrow(momac_expr), ncol(momac_expr)))
write.table(momac_expr, "tcga_bayesprism_momac_expression.tsv", sep="\t", quote=FALSE)
cat("Saved: tcga_bayesprism_momac_expression.tsv\n")

# Key gene checks
key_genes <- c("IL1B", "IL6", "OSM", "TNF", "NFKB1", "CEACAM5", "CEACAM6")
cat("\nKey gene expression (Fibroblast-specific):\n")
for (g in key_genes) {
    if (g %in% colnames(fib_expr)) {
        cat(sprintf("  %s: mean=%.4f, median=%.4f\n", g, mean(fib_expr[, g]), median(fib_expr[, g])))
    }
}
cat("\nKey gene expression (MoMac-specific):\n")
for (g in key_genes) {
    if (g %in% colnames(momac_expr)) {
        cat(sprintf("  %s: mean=%.4f, median=%.4f\n", g, mean(momac_expr[, g]), median(momac_expr[, g])))
    }
}

# Save full result
saveRDS(bp_result, "tcga_bayesprism_result.rds")
cat("\nSaved: tcga_bayesprism_result.rds\n")
cat("\nStep 2 complete!\n")
