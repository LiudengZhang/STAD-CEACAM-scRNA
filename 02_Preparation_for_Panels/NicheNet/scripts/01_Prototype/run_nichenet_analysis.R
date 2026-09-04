#!/usr/bin/env Rscript
# Paths below refer to the upstream Round_4 processing pipeline, which is
# not part of this release. This script is included as a record of how the
# input was produced; it is not called by _run_all_panels.sh.
#
# NicheNet Analysis R Script
# Ligand-Target Prediction for Cell-Cell Communication
#
# This script performs NicheNet analysis to predict target genes
# regulated by ligands from sender cells in receiver cells.
#
# Usage:
#   Rscript run_nichenet_analysis.R <data_dir> <output_dir> <sender> <receivers>
#
# Arguments:
#   data_dir  : Directory containing expression data (CSV files)
#   output_dir: Directory for output files
#   sender    : Sender cell type name
#   receivers : Comma-separated list of receiver cell types
#

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript run_nichenet_analysis.R <data_dir> <output_dir> <sender> <receivers>\n")
  quit(status = 1)
}

data_dir <- args[1]
output_dir <- args[2]
sender_type <- args[3]
receiver_types <- unlist(strsplit(args[4], ","))

cat("================================================================================\n")
cat("NicheNet Ligand-Target Prediction Analysis\n")
cat("================================================================================\n\n")
cat("Data directory:", data_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Sender:", sender_type, "\n")
cat("Receivers:", paste(receiver_types, collapse=", "), "\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(nichenetr)
  library(tidyverse)
  library(Seurat)
  library(yaml)
  library(Matrix)
})

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load NicheNet ligand-target prior model and ligand-receptor network from local files
cat("Loading NicheNet databases from local files...\n")
db_dir <- "/path/to/Project_4_05232025/Round_4/02_Playground/02_NicheNet/00_Databases"
ligand_target_matrix <- readRDS(file.path(db_dir, "ligand_target_matrix_nsga2r_final.rds"))
lr_network <- readRDS(file.path(db_dir, "lr_network_human_21122021.rds"))
weighted_networks <- readRDS(file.path(db_dir, "weighted_networks_nsga2r_final.rds"))

cat("  Ligand-target matrix:", nrow(ligand_target_matrix), "ligands x",
    ncol(ligand_target_matrix), "targets\n")
cat("  Ligand-receptor network:", nrow(lr_network), "interactions\n\n")

# Load MTX paths metadata
cat("Loading MTX paths metadata...\n")
mtx_paths_file <- file.path(data_dir, "mtx_paths.yaml")

if (!file.exists(mtx_paths_file)) {
  cat("ERROR: MTX paths file not found:", mtx_paths_file, "\n")
  quit(status = 1)
}

mtx_paths <- read_yaml(mtx_paths_file)

# Load expression data for sender from MTX
cat("Loading expression data...\n")

if (!sender_type %in% names(mtx_paths)) {
  cat("ERROR: Sender cell type not found in MTX paths:", sender_type, "\n")
  quit(status = 1)
}

sender_mtx_dir <- mtx_paths[[sender_type]]
cat("  Sender MTX directory:", sender_mtx_dir, "\n")

# Read MTX files using Seurat
sender_mat <- ReadMtx(
  mtx = file.path(sender_mtx_dir, "matrix.mtx"),
  features = file.path(sender_mtx_dir, "features.tsv"),
  cells = file.path(sender_mtx_dir, "barcodes.tsv"),
  feature.column = 1
)

# Convert to dense matrix and transpose (cells x genes format)
sender_expr <- as.data.frame(t(as.matrix(sender_mat)))
cat("  Sender (", sender_type, "):", nrow(sender_expr), "cells x", ncol(sender_expr), "genes\n")

# Identify expressed ligands in sender cells
# First filter to ligands in the LR network
lr_ligands <- unique(lr_network$from)
ligands_in_data <- intersect(lr_ligands, colnames(sender_expr))

if (length(ligands_in_data) > 0) {
  # Calculate expression fraction for these ligands
  expr_fraction <- colMeans(sender_expr[, ligands_in_data, drop=FALSE] > 0)
  expressed_ligands <- names(expr_fraction[expr_fraction > 0.1])
  # Also ensure they're in the ligand-target matrix
  expressed_ligands <- intersect(expressed_ligands, rownames(ligand_target_matrix))
} else {
  expressed_ligands <- character(0)
}

cat("  LR network ligands in data:", length(ligands_in_data), "\n")
cat("  Expressed ligands (>10% cells):", length(expressed_ligands), "\n\n")

# Analyze each receiver cell type
for (receiver_type in receiver_types) {
  cat("Analyzing receiver:", receiver_type, "\n")

  if (!receiver_type %in% names(mtx_paths)) {
    cat("  WARNING: Receiver cell type not found in MTX paths, skipping\n\n")
    next
  }

  receiver_mtx_dir <- mtx_paths[[receiver_type]]
  cat("  Receiver MTX directory:", receiver_mtx_dir, "\n")

  # Read MTX files using Seurat
  receiver_mat <- ReadMtx(
    mtx = file.path(receiver_mtx_dir, "matrix.mtx"),
    features = file.path(receiver_mtx_dir, "features.tsv"),
    cells = file.path(receiver_mtx_dir, "barcodes.tsv"),
    feature.column = 1
  )

  # Convert to dense matrix and transpose (cells x genes format)
  receiver_expr <- as.data.frame(t(as.matrix(receiver_mat)))
  cat("  Loaded:", nrow(receiver_expr), "cells x", ncol(receiver_expr), "genes\n")

  # Identify variable genes as potential targets
  gene_vars <- apply(receiver_expr, 2, var)
  target_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(200, length(gene_vars))]
  target_genes <- intersect(target_genes, colnames(ligand_target_matrix))

  cat("  Target genes:", length(target_genes), "\n")

  # Calculate ligand activities
  if (length(target_genes) > 10 && length(expressed_ligands) > 5) {
    ligand_activities <- predict_ligand_activities(
      geneset = target_genes,
      background_expressed_genes = colnames(receiver_expr),
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = expressed_ligands
    )

    # Save results
    result_file <- file.path(output_dir,
                            paste0(sender_type, "_to_", receiver_type, "_ligand_activities.csv"))
    write.csv(ligand_activities, result_file, row.names = FALSE)

    cat("  Top ligands:\n")
    print(head(ligand_activities %>% arrange(desc(pearson)), 10))
    cat("  Results saved to:", result_file, "\n\n")
  } else {
    cat("  WARNING: Insufficient genes for analysis, skipping\n\n")
  }
}

cat("================================================================================\n")
cat("NicheNet analysis completed\n")
cat("================================================================================\n")
