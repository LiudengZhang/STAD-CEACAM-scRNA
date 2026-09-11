#!/usr/bin/env Rscript
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
#
# NicheNet Analysis: C3_Mac_Inflam_IL1B as Sender
#
# Predicts ligand-target interactions from C3 inflammatory macrophages
# to 13 receiver cell types using final Round 4.2 annotations.
#
# Filter: Post-treatment non-responder (stomach_post_grouping == No-response)
# Uses local NicheNet databases (not Zenodo URLs).
#

cat("================================================================================\n")
cat("NicheNet Analysis: C3_Mac_Inflam_IL1B -> All Receivers\n")
cat("Filter: Post-treatment Non-Responder\n")
cat("================================================================================\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(nichenetr)
  library(tidyverse)
  library(Seurat)
  library(yaml)
  library(Matrix)
})

# Get script directory
args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script_dir <- dirname(script_path)
app_dir <- dirname(script_dir)

# Paths
config_file <- file.path(app_dir, "01_Config", "config.yaml")
data_dir <- file.path(app_dir, "prepared_data")
output_dir <- file.path(app_dir, "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load config
config <- read_yaml(config_file)
cat("Config loaded from:", config_file, "\n")

# Sender
sender_type <- config$sender
cat("Sender:", sender_type, "\n\n")

# Load NicheNet databases (local files)
cat("Loading NicheNet databases (local)...\n")
ligand_target_matrix <- readRDS(config$databases$ligand_target_matrix)
lr_network <- readRDS(config$databases$lr_network)
cat("  Ligand-target matrix:", nrow(ligand_target_matrix), "ligands x",
    ncol(ligand_target_matrix), "targets\n")
cat("  Ligand-receptor network:", nrow(lr_network), "interactions\n\n")

# Load MTX paths
mtx_paths_file <- file.path(data_dir, "mtx_paths.yaml")
mtx_paths <- read_yaml(mtx_paths_file)

# Function to read MTX data
read_mtx_data <- function(mtx_dir) {
  mat <- ReadMtx(
    mtx = file.path(mtx_dir, "matrix.mtx"),
    features = file.path(mtx_dir, "features.tsv"),
    cells = file.path(mtx_dir, "barcodes.tsv"),
    feature.column = 1
  )
  return(as.data.frame(t(as.matrix(mat))))
}

# Load sender expression
cat("Loading sender expression data...\n")
sender_dir <- mtx_paths[[sender_type]]
sender_expr <- read_mtx_data(sender_dir)
cat("  Sender (", sender_type, "):", nrow(sender_expr), "cells x", ncol(sender_expr), "genes\n")

# Identify expressed ligands in sender (must be in NicheNet database)
min_expr_frac <- config$nichenet$min_expression_fraction
expressed_genes_sender <- colnames(sender_expr)[colMeans(sender_expr > 0) > min_expr_frac]
expressed_ligands <- intersect(expressed_genes_sender, rownames(ligand_target_matrix))

# Also ensure ligands are in the L-R network
lr_ligands <- unique(lr_network$from)
expressed_ligands <- intersect(expressed_ligands, lr_ligands)

cat("  Expressed genes (>", min_expr_frac*100, "% cells):", length(expressed_genes_sender), "\n")
cat("  Valid NicheNet ligands:", length(expressed_ligands), "\n\n")

# Get receivers from config
receivers <- config$receivers

# Analyze each receiver
all_results <- list()

for (receiver_type in receivers) {
  cat("--------------------------------------------------------------------------------\n")
  cat("Analyzing receiver:", receiver_type, "\n")

  if (!receiver_type %in% names(mtx_paths)) {
    cat("  WARNING: Receiver not found in MTX paths, skipping\n\n")
    next
  }

  receiver_dir <- mtx_paths[[receiver_type]]
  receiver_expr <- read_mtx_data(receiver_dir)
  cat("  Loaded:", nrow(receiver_expr), "cells x", ncol(receiver_expr), "genes\n")

  # Identify expressed receptors
  expressed_receptors <- colnames(receiver_expr)[colMeans(receiver_expr > 0) > min_expr_frac]

  # Get expressed L-R pairs
  lr_pairs <- lr_network %>%
    filter(from %in% expressed_ligands & to %in% expressed_receptors)
  cat("  Expressed L-R pairs:", nrow(lr_pairs), "\n")

  # Identify variable genes as potential targets (must be in NicheNet targets)
  n_target_genes <- config$nichenet$n_target_genes
  nichenet_targets <- colnames(ligand_target_matrix)
  receiver_nichenet_genes <- intersect(colnames(receiver_expr), nichenet_targets)

  # Get most variable among NicheNet targets
  if (length(receiver_nichenet_genes) > 0) {
    gene_vars <- apply(receiver_expr[, receiver_nichenet_genes, drop=FALSE], 2, var)
    target_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_target_genes, length(gene_vars))]
  } else {
    target_genes <- character(0)
  }

  cat("  NicheNet target genes in receiver:", length(receiver_nichenet_genes), "\n")
  cat("  Selected target genes:", length(target_genes), "\n")

  # Calculate ligand activities
  if (length(target_genes) >= 10 && length(expressed_ligands) >= 5) {
    # Background = ALL expressed genes in receiver (NicheNet internally intersects
    # with ligand_target_matrix columns). Must be superset of target_genes.
    background_genes <- colnames(receiver_expr)[colMeans(receiver_expr > 0) > min_expr_frac]

    tryCatch({
      ligand_activities <- predict_ligand_activities(
        geneset = target_genes,
        background_expressed_genes = background_genes,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = expressed_ligands
      )

      # Sort by pearson
      ligand_activities <- ligand_activities %>% arrange(desc(pearson))

      # Get top ligands
      top_n <- config$nichenet$top_n_ligands
      top_ligands <- head(ligand_activities$test_ligand, top_n)

      cat("  Top 5 ligands:\n")
      print(head(ligand_activities, 5))

      # Save results per receiver
      safe_receiver <- gsub("[^A-Za-z0-9]", "_", receiver_type)

      # Create per-receiver output dir
      receiver_out_dir <- file.path(output_dir, safe_receiver)
      dir.create(receiver_out_dir, recursive = TRUE, showWarnings = FALSE)

      # Ligand activities
      write.csv(
        ligand_activities,
        file.path(receiver_out_dir, "ligand_activities.csv"),
        row.names = FALSE
      )

      # L-R interactions
      write.csv(
        lr_pairs,
        file.path(receiver_out_dir, "lr_interactions.csv"),
        row.names = FALSE
      )

      # Ligand-target links for top ligands
      if (length(top_ligands) > 0) {
        active_ligand_target_links <- ligand_target_matrix[top_ligands, target_genes, drop = FALSE]
        active_ligand_target_df <- as.data.frame(as.matrix(active_ligand_target_links)) %>%
          rownames_to_column("ligand") %>%
          pivot_longer(-ligand, names_to = "target", values_to = "weight") %>%
          filter(weight > 0) %>%
          arrange(desc(weight))

        write.csv(
          active_ligand_target_df,
          file.path(receiver_out_dir, "ligand_targets.csv"),
          row.names = FALSE
        )
      }

      all_results[[receiver_type]] <- ligand_activities
      cat("  Results saved to:", receiver_out_dir, "\n\n")
    }, error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      cat("  Skipping this receiver.\n\n")
    })
  } else {
    cat("  WARNING: Insufficient genes/ligands for analysis, skipping\n\n")
  }
}

# Aggregate results: raw CSV (all ligand-receiver pairs)
cat("================================================================================\n")
cat("Aggregating results...\n")
cat("================================================================================\n")

raw_df <- do.call(rbind, lapply(names(all_results), function(receiver) {
  df <- all_results[[receiver]]
  df$receiver <- receiver
  df
}))

write.csv(raw_df,
          file.path(output_dir, "c3mac_ligand_activities_raw.csv"),
          row.names = FALSE)
cat("Raw activities saved:", nrow(raw_df), "rows\n")

# Aggregate: mean AUROC across receivers per ligand
agg_df <- raw_df %>%
  group_by(test_ligand) %>%
  summarise(
    auroc_mean = mean(auroc, na.rm = TRUE),
    auroc_std = sd(auroc, na.rm = TRUE),
    n_receivers = n(),
    pearson_mean = mean(pearson, na.rm = TRUE),
    pearson_std = sd(pearson, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(ligand = test_ligand) %>%
  arrange(desc(auroc_mean))

write.csv(agg_df,
          file.path(output_dir, "c3mac_ligand_activities_aggregated.csv"),
          row.names = FALSE)
cat("Aggregated activities saved:", nrow(agg_df), "unique ligands\n")

# Summary
cat("\n================================================================================\n")
cat("Summary\n")
cat("================================================================================\n")

summary_df <- data.frame(
  receiver = names(all_results),
  n_ligands = sapply(all_results, nrow),
  top_ligand = sapply(all_results, function(x) x$test_ligand[1]),
  top_pearson = sapply(all_results, function(x) x$pearson[1])
)

write.csv(summary_df, file.path(output_dir, "summary.csv"), row.names = FALSE)
print(summary_df)

cat("\nResults saved to:", output_dir, "\n")
cat("================================================================================\n")
cat("NicheNet C3_Mac analysis completed\n")
cat("================================================================================\n")
