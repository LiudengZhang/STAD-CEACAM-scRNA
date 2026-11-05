#!/usr/bin/env Rscript
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# ============================================================================
# Script: 02_run_estimate.R
# Description: Run ESTIMATE algorithm to calculate tumor purity scores
# Date: 2025-10-14
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
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
ESTIMATE_DIR <- file.path(BASE_DIR, "02_Tumor_Purity_Analysis/01_ESTIMATE_Results")
LOG_DIR <- file.path(BASE_DIR, "03_Logs")

# Set up logging
log_file <- file.path(LOG_DIR, "02_estimate_analysis.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat("============================================================\n")
cat("TUMOR PURITY ANALYSIS - ESTIMATE\n")
cat("============================================================\n\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# ============================================================================
# 1. Check and Install ESTIMATE Package
# ============================================================================
cat("Step 1: Checking ESTIMATE package...\n")

# ESTIMATE is available from sourceforge
# Download link: https://sourceforge.net/projects/estimateproject/
estimate_installed <- "estimate" %in% rownames(installed.packages())

if (!estimate_installed) {
  cat("  - ESTIMATE package not found\n")
  cat("  - Attempting to install ESTIMATE...\n")

  # Try to download and install ESTIMATE
  estimate_url <- "https://sourceforge.net/projects/estimateproject/files/estimate_1.0.13.tar.gz/download"
  temp_file <- tempfile(fileext = ".tar.gz")

  tryCatch({
    download.file(estimate_url, temp_file, mode = "wb")
    install.packages(temp_file, repos = NULL, type = "source")
    cat("  - ESTIMATE installed successfully\n")
  }, error = function(e) {
    cat("  ERROR: Could not install ESTIMATE automatically\n")
    cat("  Please install manually from: https://sourceforge.net/projects/estimateproject/\n")
    cat("  Error message:", conditionMessage(e), "\n")
    stop("ESTIMATE package installation failed")
  })
} else {
  cat("  - ESTIMATE package found\n")
}

library(estimate)

# ============================================================================
# 2. Load Prepared Data
# ============================================================================
cat("\nStep 2: Loading prepared data...\n")

# Load expression matrix
expr_file <- file.path(ESTIMATE_DIR, "expression_matrix_for_estimate.rds")
expr_matrix <- readRDS(expr_file)
cat("  - Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")

# Load metadata
meta_file <- file.path(ESTIMATE_DIR, "metadata.rds")
metadata <- readRDS(meta_file)
cat("  - Metadata:", nrow(metadata), "samples\n")

# ============================================================================
# 3. Prepare Input File for ESTIMATE
# ============================================================================
cat("\nStep 3: Preparing input file for ESTIMATE...\n")

# ESTIMATE requires a GCT format file
# Create a temporary GCT file
gct_file <- file.path(ESTIMATE_DIR, "expression_for_estimate.gct")

# Write GCT header
cat("#1.2\n", file = gct_file)
cat(nrow(expr_matrix), "\t", ncol(expr_matrix), "\n", file = gct_file, append = TRUE)

# Create data frame with NAME and Description columns
gct_data <- data.frame(
  NAME = rownames(expr_matrix),
  Description = rownames(expr_matrix),
  expr_matrix,
  check.names = FALSE
)

# Write data
write.table(gct_data, gct_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE, append = TRUE)

cat("  - GCT file created:", gct_file, "\n")

# ============================================================================
# 4. Run ESTIMATE Algorithm
# ============================================================================
cat("\nStep 4: Running ESTIMATE algorithm...\n")
cat("  This may take several minutes...\n\n")

# Define output file
estimate_output <- file.path(ESTIMATE_DIR, "estimate_scores.gct")

# Run ESTIMATE
# Platform: affymetrix or illumina or agilent (we use illumina for RNA-seq)
tryCatch({
  estimateScore(
    input.ds = gct_file,
    output.ds = estimate_output,
    platform = "illumina"
  )
  cat("  - ESTIMATE completed successfully\n")
}, error = function(e) {
  cat("  ERROR in ESTIMATE calculation:\n")
  cat("  ", conditionMessage(e), "\n")
  stop("ESTIMATE calculation failed")
})

# ============================================================================
# 5. Parse ESTIMATE Results
# ============================================================================
cat("\nStep 5: Parsing ESTIMATE results...\n")

# Read ESTIMATE output (GCT format)
estimate_lines <- readLines(estimate_output)
# Skip first two lines (GCT header)
estimate_data <- read.table(text = estimate_lines[3:length(estimate_lines)],
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("  - ESTIMATE output dimensions:", nrow(estimate_data), "x", ncol(estimate_data), "\n")

# Extract scores
# ESTIMATE output has rows for: StromalScore, ImmuneScore, ESTIMATEScore
score_names <- estimate_data$NAME

# Create data frame with scores for each sample
scores_df <- data.frame(
  sample_id = colnames(estimate_data)[3:ncol(estimate_data)],
  stringsAsFactors = FALSE
)

# Extract each score type
for (score_type in score_names) {
  score_values <- as.numeric(estimate_data[estimate_data$NAME == score_type, 3:ncol(estimate_data)])
  scores_df[[score_type]] <- score_values
}

# Calculate tumor purity using the formula from ESTIMATE paper
# TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATE score)
scores_df$TumorPurity <- cos(0.6049872018 + 0.0001467884 * scores_df$ESTIMATEScore)

cat("\nScore summary:\n")
cat("  - Stromal Score range:", range(scores_df$StromalScore), "\n")
cat("  - Immune Score range:", range(scores_df$ImmuneScore), "\n")
cat("  - ESTIMATE Score range:", range(scores_df$ESTIMATEScore), "\n")
cat("  - Tumor Purity range:", range(scores_df$TumorPurity), "\n")

# ============================================================================
# 6. Merge with Metadata
# ============================================================================
cat("\nStep 6: Merging scores with metadata...\n")

# Merge scores with metadata
metadata_with_scores <- merge(metadata, scores_df, by = "sample_id", all.x = TRUE)

cat("  - Merged data dimensions:", nrow(metadata_with_scores), "x", ncol(metadata_with_scores), "\n")

# Check for missing values
missing_scores <- sum(is.na(metadata_with_scores$TumorPurity))
if (missing_scores > 0) {
  cat("  WARNING:", missing_scores, "samples have missing purity scores\n")
} else {
  cat("  - All samples have purity scores\n")
}

# ============================================================================
# 7. Save Results
# ============================================================================
cat("\nStep 7: Saving results...\n")

# Save scores only
scores_file <- file.path(ESTIMATE_DIR, "tumor_purity_scores.tsv")
write.table(scores_df, scores_file, sep = "\t", quote = FALSE,
            row.names = FALSE)
cat("  - Scores table:", scores_file, "\n")

# Save RDS
scores_rds <- file.path(ESTIMATE_DIR, "tumor_purity_scores.rds")
saveRDS(scores_df, scores_rds)
cat("  - Scores RDS:", scores_rds, "\n")

# Save merged metadata with scores
metadata_scores_file <- file.path(ESTIMATE_DIR, "metadata_with_purity_scores.tsv")
write.table(metadata_with_scores, metadata_scores_file, sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("  - Metadata with scores:", metadata_scores_file, "\n")

# Save RDS
metadata_scores_rds <- file.path(ESTIMATE_DIR, "metadata_with_purity_scores.rds")
saveRDS(metadata_with_scores, metadata_scores_rds)
cat("  - Metadata RDS:", metadata_scores_rds, "\n")

# ============================================================================
# 8. Generate Summary Statistics
# ============================================================================
cat("\n============================================================\n")
cat("ESTIMATE ANALYSIS SUMMARY\n")
cat("============================================================\n\n")

cat("Overall statistics:\n")
cat("  Number of samples:", nrow(scores_df), "\n\n")

cat("Stromal Score:\n")
print(summary(scores_df$StromalScore))
cat("\n")

cat("Immune Score:\n")
print(summary(scores_df$ImmuneScore))
cat("\n")

cat("ESTIMATE Score:\n")
print(summary(scores_df$ESTIMATEScore))
cat("\n")

cat("Tumor Purity:\n")
print(summary(scores_df$TumorPurity))
cat("\n")

cat("Scores by response status:\n")
by_response <- metadata_with_scores %>%
  group_by(response_NR) %>%
  summarise(
    n = n(),
    mean_purity = mean(TumorPurity, na.rm = TRUE),
    sd_purity = sd(TumorPurity, na.rm = TRUE),
    mean_stromal = mean(StromalScore, na.rm = TRUE),
    mean_immune = mean(ImmuneScore, na.rm = TRUE),
    .groups = 'drop'
  )
print(by_response)

cat("\nEnd time:", as.character(Sys.time()), "\n")
cat("============================================================\n")

# Close logging
sink(type = "message")
sink(type = "output")
close(log_con)

cat("ESTIMATE analysis completed successfully!\n")
cat("Log file:", log_file, "\n")
