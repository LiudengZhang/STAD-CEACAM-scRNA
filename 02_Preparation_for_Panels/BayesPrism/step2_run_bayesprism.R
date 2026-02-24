#!/usr/bin/env Rscript
# Step 2: Run BayesPrism deconvolution on TIGER bulk data
#
# Input:  sc_counts.tsv, sc_cell_types.tsv, bulk_counts.tsv
# Output: bayesprism_fractions.tsv, bayesprism_epithelial_expression.tsv

suppressPackageStartupMessages({
    library(BayesPrism)
    library(data.table)
})

cat("============================================================\n")
cat("Step 2: BayesPrism Deconvolution\n")
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

# --- Load data with fread (much faster than read.table) ---
cat("\nLoading scRNA-seq reference (fread)...\n")
sc_dt <- fread("sc_counts.tsv", sep="\t", header=TRUE)
sc_rownames <- sc_dt[[1]]
sc_dt[, 1 := NULL]
sc_mat <- as.matrix(sc_dt)
rownames(sc_mat) <- sc_rownames
rm(sc_dt); gc()
cat(sprintf("  scRNA: %d cells x %d genes\n", nrow(sc_mat), ncol(sc_mat)))

sc_labels <- fread("sc_cell_types.tsv", sep="\t", header=TRUE)
cat(sprintf("  Cell types: %s\n", paste(unique(sc_labels$cell_type), collapse=", ")))

cat("\nLoading bulk data (fread)...\n")
bk_dt <- fread("bulk_counts.tsv", sep="\t", header=TRUE)
bk_rownames <- bk_dt[[1]]
bk_dt[, 1 := NULL]
bk_mat <- as.matrix(bk_dt)
rownames(bk_mat) <- bk_rownames
rm(bk_dt); gc()
cat(sprintf("  Bulk: %d samples x %d genes\n", nrow(bk_mat), ncol(bk_mat)))

# Cell type and state labels
# Use cell_type for both to avoid sparse state issue (some states have <10 cells)
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

# Free memory
rm(sc_mat, bk_mat); gc()

# --- Run BayesPrism ---
cat("\nRunning BayesPrism (this may take a while)...\n")
bp_result <- run.prism(prism = myPrism, n.cores = 4)

# --- Extract results ---
cat("\nExtracting results...\n")

# 1. Cell type fractions
theta <- get.fraction(bp = bp_result,
                      which.theta = "final",
                      state.or.type = "type")
cat(sprintf("Fractions: %d samples x %d types\n", nrow(theta), ncol(theta)))
write.table(theta, "bayesprism_fractions.tsv", sep="\t", quote=FALSE)
cat("Saved: bayesprism_fractions.tsv\n")

# 2. Epithelial-specific expression
cat("\nExtracting epithelial-specific expression...\n")
epi_expr <- get.exp(bp = bp_result,
                    state.or.type = "type",
                    cell.name = "Epithelial cells")
cat(sprintf("Epithelial expression: %d samples x %d genes\n", nrow(epi_expr), ncol(epi_expr)))

# Save full epithelial expression
write.table(epi_expr, "bayesprism_epithelial_expression.tsv", sep="\t", quote=FALSE)
cat("Saved: bayesprism_epithelial_expression.tsv\n")

# Check CEACAM genes
ceacam_genes <- grep("^CEACAM[0-9]", colnames(epi_expr), value=TRUE)
cat(sprintf("\nCEACAM genes found: %s\n", paste(ceacam_genes, collapse=", ")))

if (length(ceacam_genes) > 0) {
    ceacam_df <- epi_expr[, ceacam_genes, drop=FALSE]
    write.table(ceacam_df, "bayesprism_epithelial_ceacam.tsv", sep="\t", quote=FALSE)
    cat("Saved: bayesprism_epithelial_ceacam.tsv\n")

    cat("\nCEACAM expression summary (epithelial-specific):\n")
    for (g in ceacam_genes) {
        cat(sprintf("  %s: mean=%.4f, median=%.4f\n", g, mean(epi_expr[, g]), median(epi_expr[, g])))
    }
}

# Save result object
saveRDS(bp_result, "bayesprism_result.rds")
cat("\nSaved: bayesprism_result.rds\n")
cat("\nStep 2 complete!\n")
