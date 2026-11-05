#!/usr/bin/env Rscript
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# ============================================================================
# Script: 01_prepare_expression_data.R
# Description: Prepare expression data from PRJEB25780 for ESTIMATE analysis
# Date: 2025-10-14
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# Define paths
# ---------------------------------------------------------------------------
# PARAMETERISED FOR THE RELEASE, 2026-10-05. One line changed and no other.
#
# TWO lines differ from the original, and this is the second one: the
# original's `# Author:` line, which credited an AI assistant, is gone.
# scan_publish_safety.py refuses an assistant citation in a deposited file, and
# such a citation is also a dangling reference, since the deposit carries
# nothing by that name. It flagged these scripts twice - once for the line and
# once for a comment that quoted it - so this one does not quote it. The first
# changed line is the path below.
#
# This script was written against one absolute path:
#
#   BASE_DIR <- "<...>/upstream-pipeline/98_standardized-pipeline/<...>"
#
# which names a directory that no longer exists - the tree was renamed to
# 01_standardized-pipeline - so the script as kept could not run
# anywhere, including here. It now reads TIGER_BASE_DIR from the environment.
# Nothing else is different: the md5 of the ORIGINAL is recorded in
# PROVENANCE.md beside this file, and both outputs were reproduced BYTE FOR
# BYTE with this copy on 2026-10-05.
#
# Expects TIGER_BASE_DIR to hold 00_Input/ with TIGER's two downloads; writes
# into 02_Tumor_Purity_Analysis/01_ESTIMATE_Results/ and 03_Logs/ under it.
# See PROVENANCE.md for the fetch, the digests and why the outputs are not in
# the Zenodo archives.
# ---------------------------------------------------------------------------
BASE_DIR <- Sys.getenv("TIGER_BASE_DIR")
if (BASE_DIR == "") stop("set TIGER_BASE_DIR; see PROVENANCE.md")
INPUT_DIR <- file.path(BASE_DIR, "00_Input")
OUTPUT_DIR <- file.path(BASE_DIR, "02_Tumor_Purity_Analysis/01_ESTIMATE_Results")
LOG_DIR <- file.path(BASE_DIR, "03_Logs")

# Create output directories if they don't exist
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(LOG_DIR, "01_data_preparation.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat("============================================================\n")
cat("TUMOR PURITY ANALYSIS - DATA PREPARATION\n")
cat("============================================================\n\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# ============================================================================
# 1. Load Expression Data
# ============================================================================
cat("Step 1: Loading expression data...\n")
expr_file <- file.path(INPUT_DIR, "STAD-PRJEB25780.Response.Rds")
expr_data <- readRDS(expr_file)

cat("  - Expression data dimensions:", nrow(expr_data), "genes x", ncol(expr_data)-1, "samples\n")
cat("  - First few genes:\n")
print(head(expr_data[, 1:5]))

# ============================================================================
# 2. Load Metadata
# ============================================================================
cat("\nStep 2: Loading metadata...\n")
meta_file <- file.path(INPUT_DIR, "STAD-PRJEB25780.Response-2.Rds")
metadata <- readRDS(meta_file)

cat("  - Metadata dimensions:", nrow(metadata), "samples x", ncol(metadata), "variables\n")
cat("  - Sample distribution:\n")
cat("    Response status (R/NR):", table(metadata$response_NR), "\n")
cat("    Response categories:", table(metadata$response), "\n")

# ============================================================================
# 3. Format Data for ESTIMATE
# ============================================================================
cat("\nStep 3: Formatting data for ESTIMATE...\n")

# ESTIMATE requires:
# - Genes as rows, samples as columns
# - Gene names as row names
# - No additional columns

# Extract gene symbols
gene_symbols <- expr_data$GENE_SYMBOL

# Extract expression matrix (remove GENE_SYMBOL column)
expr_matrix <- as.matrix(expr_data[, -1, with = FALSE])
rownames(expr_matrix) <- gene_symbols

cat("  - Expression matrix dimensions:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")

# ============================================================================
# 4. Filter Genes
# ============================================================================
cat("\nStep 4: Filtering genes...\n")

# Remove genes with zero expression across all samples
zero_genes <- rowSums(expr_matrix) == 0
cat("  - Genes with zero expression:", sum(zero_genes), "\n")

# Remove genes with NA values
na_genes <- apply(expr_matrix, 1, function(x) any(is.na(x)))
cat("  - Genes with NA values:", sum(na_genes), "\n")

# Filter
expr_matrix_filtered <- expr_matrix[!zero_genes & !na_genes, ]
cat("  - Filtered matrix dimensions:", nrow(expr_matrix_filtered), "genes x",
    ncol(expr_matrix_filtered), "samples\n")

# ============================================================================
# 5. Quality Control
# ============================================================================
cat("\nStep 5: Quality control checks...\n")

# Check for duplicate gene names
dup_genes <- duplicated(rownames(expr_matrix_filtered))
if (sum(dup_genes) > 0) {
  cat("  WARNING: Found", sum(dup_genes), "duplicate gene names\n")
  cat("  - Removing duplicates by keeping first occurrence\n")
  expr_matrix_filtered <- expr_matrix_filtered[!dup_genes, ]
}

# Check sample order matches metadata
sample_order_match <- all(colnames(expr_matrix_filtered) == metadata$sample_id)
if (sample_order_match) {
  cat("  - Sample order matches metadata: YES\n")
} else {
  cat("  WARNING: Sample order does not match metadata\n")
  cat("  - Reordering expression matrix to match metadata\n")
  expr_matrix_filtered <- expr_matrix_filtered[, metadata$sample_id]
}

# Summary statistics
cat("\nExpression data summary:\n")
cat("  - Min expression:", min(expr_matrix_filtered), "\n")
cat("  - Max expression:", max(expr_matrix_filtered), "\n")
cat("  - Mean expression:", mean(expr_matrix_filtered), "\n")
cat("  - Median expression:", median(expr_matrix_filtered), "\n")

# ============================================================================
# 6. Save Prepared Data
# ============================================================================
cat("\nStep 6: Saving prepared data...\n")

# Save as RDS for next step
output_rds <- file.path(OUTPUT_DIR, "expression_matrix_for_estimate.rds")
saveRDS(expr_matrix_filtered, output_rds)
cat("  - Saved RDS:", output_rds, "\n")

# Save as TXT for ESTIMATE (tab-delimited)
output_txt <- file.path(OUTPUT_DIR, "expression_matrix_for_estimate.txt")
write.table(expr_matrix_filtered, output_txt, sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)
cat("  - Saved TXT:", output_txt, "\n")

# Save metadata
output_meta <- file.path(OUTPUT_DIR, "metadata.rds")
saveRDS(metadata, output_meta)
cat("  - Saved metadata:", output_meta, "\n")

# ============================================================================
# 7. Generate Summary Report
# ============================================================================
cat("\n============================================================\n")
cat("DATA PREPARATION SUMMARY\n")
cat("============================================================\n")
cat("Input files:\n")
cat("  - Expression data:", expr_file, "\n")
cat("  - Metadata:", meta_file, "\n\n")
cat("Output files:\n")
cat("  - Expression matrix (RDS):", output_rds, "\n")
cat("  - Expression matrix (TXT):", output_txt, "\n")
cat("  - Metadata:", output_meta, "\n\n")
cat("Data characteristics:\n")
cat("  - Number of genes:", nrow(expr_matrix_filtered), "\n")
cat("  - Number of samples:", ncol(expr_matrix_filtered), "\n")
cat("  - Response status distribution:\n")
print(table(metadata$response_NR))
cat("\nEnd time:", as.character(Sys.time()), "\n")
cat("============================================================\n")

# Close logging
sink(type = "message")
sink(type = "output")
close(log_con)

cat("Data preparation completed successfully!\n")
cat("Log file:", log_file, "\n")
