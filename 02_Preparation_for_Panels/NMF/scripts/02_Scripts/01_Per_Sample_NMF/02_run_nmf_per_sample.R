#!/usr/bin/env Rscript
# Paths below refer to the upstream Round_4 processing pipeline, which is
# not part of this release. This script is included as a record of how the
# input was produced; it is not called by _run_all_panels.sh.
# 02_run_nmf_per_sample.R - Run NMF per sample (K=4-9)

library(NMF)
library(Matrix)
library(yaml)

# Hardcoded base directory (robust path resolution)
base_dir <- "/path/to/Project_4_05232025/Round_4/999_Panels/Figure_02.1/02_A/3CA_MetaProgram"
config_path <- file.path(base_dir, "00_Config", "config.yaml")
input_dir <- file.path(base_dir, "03_Output", "01_Per_Sample_NMF")

config <- yaml::read_yaml(config_path)
k_values <- config$nmf$k_values
n_iter <- config$nmf$n_iter
seed <- config$nmf$random_seed
top_genes <- config$nmf$top_genes_per_program
nmf_method <- config$nmf$method

cat("PHASE 1.2: RUN PER-SAMPLE NMF\n")

sample_summary <- read.csv(file.path(input_dir, "sample_summary.csv"))
samples <- sample_summary$sample_id

for (sample in samples) {
    cat("\nProcessing:", sample, "\n")
    counts_file <- file.path(input_dir, sample, paste0(sample, "_counts.csv"))
    if (!file.exists(counts_file)) next

    counts <- read.csv(counts_file, row.names=1, check.names=FALSE)
    counts_t <- t(as.matrix(counts))
    counts_t[counts_t < 0] <- 0

    # Remove genes with zero expression across all cells
    row_sums <- rowSums(counts_t)
    counts_t <- counts_t[row_sums > 0, , drop=FALSE]

    col_sums <- colSums(counts_t)
    col_sums[col_sums == 0] <- 1
    counts_norm <- sweep(counts_t, 2, col_sums, "/") * 1e6
    counts_log <- log2(counts_norm + 1)

    # Remove any remaining rows with NA or all zeros after log transform
    row_means <- rowMeans(counts_log)
    counts_log <- counts_log[!is.na(row_means) & row_means > 0, , drop=FALSE]

    cat("  Matrix dimensions:", nrow(counts_log), "genes x", ncol(counts_log), "cells\n")

    for (k in k_values) {
        cat("  K =", k, "\n")
        set.seed(seed)
        nmf_result <- tryCatch({
            nmf(counts_log, rank=k, method=nmf_method, nrun=n_iter, seed=seed, .options='v')
        }, error = function(e) {
            nmf(counts_log, rank=k, method=nmf_method, nrun=10, seed=seed)
        })

        W <- basis(nmf_result)
        H <- coef(nmf_result)

        program_genes <- list()
        for (i in 1:k) {
            gene_scores <- W[, i]
            names(gene_scores) <- rownames(W)
            top_idx <- order(gene_scores, decreasing=TRUE)[1:min(top_genes, length(gene_scores))]
            program_genes[[paste0("P", i)]] <- names(gene_scores)[top_idx]
        }

        W_df <- as.data.frame(W)
        colnames(W_df) <- paste0("Program_", 1:k)
        write.csv(W_df, file.path(input_dir, sample, paste0(sample, "_nmf_k", k, "_W.csv")))

        H_df <- as.data.frame(t(H))
        colnames(H_df) <- paste0("Program_", 1:k)
        write.csv(H_df, file.path(input_dir, sample, paste0(sample, "_nmf_k", k, "_H.csv")))

        genes_df <- do.call(cbind, lapply(program_genes, function(x) c(x, rep(NA, top_genes - length(x)))))
        write.csv(genes_df, file.path(input_dir, sample, paste0(sample, "_nmf_k", k, "_genes.csv")), row.names=FALSE)
    }
}

cat("\nPHASE 1.2 COMPLETE\n")
