#!/usr/bin/env Rscript
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
# Do the re-run and the deposited run agree on what the panels actually read?
#
# The basis matrix agrees to 6.4e-12, which is BLAS noise rather than a
# different answer - the original run's thread count was never recorded, and
# pinning OMP_NUM_THREADS to 1 makes two of our own runs bit-identical. But W
# is not what Figure 2G and 2H consume. They read {sample}_nmf_k{k}_genes.csv,
# the top fifty genes of each program, and a difference of 1e-12 can only
# change that list if two genes tie to within 1e-12 - which is worth checking
# rather than assuming.
#
# This re-runs several samples at the settled parameters (seed 42, brunet,
# nrun 10, established by test_nrun.R against nrun 100) and compares the gene
# lists character for character.
#
# Usage: OMP_NUM_THREADS=1 Rscript verify_nmf.R [n_samples]

suppressPackageStartupMessages({library(NMF)})

DEP <- file.path("/path/to/home/project",
                 "submission-tree/02_Preparation_for_Panels/NMF/per_sample_nmf")
K_VALUES <- 4:9
TOP_GENES <- 50

args <- commandArgs(trailingOnly = TRUE)
n_samples <- if (length(args) >= 1) as.integer(args[1]) else 3L

samples <- list.dirs(DEP, recursive = FALSE, full.names = FALSE)
samples <- samples[order(sapply(samples, function(s)
    file.size(file.path(DEP, s, paste0(s, "_counts.csv")))))]
samples <- head(samples[!is.na(samples)], n_samples)

prepare <- function(sample) {
    counts <- read.csv(file.path(DEP, sample, paste0(sample, "_counts.csv")),
                       row.names = 1, check.names = FALSE)
    m <- t(as.matrix(counts))
    m[m < 0] <- 0
    m <- m[rowSums(m) > 0, , drop = FALSE]
    cs <- colSums(m); cs[cs == 0] <- 1
    m <- log2(sweep(m, 2, cs, "/") * 1e6 + 1)
    rm_ <- rowMeans(m)
    m[!is.na(rm_) & rm_ > 0, , drop = FALSE]
}

agree <- 0; differ <- 0; skipped <- 0
for (sample in samples) {
    m <- prepare(sample)
    cat(sprintf("\n%s  %d genes x %d cells\n", sample, nrow(m), ncol(m)))
    for (k in K_VALUES) {
        want_file <- file.path(DEP, sample, sprintf("%s_nmf_k%d_genes.csv",
                                                    sample, k))
        if (!file.exists(want_file)) { skipped <- skipped + 1; next }

        set.seed(42)
        W <- basis(nmf(m, rank = k, method = "brunet", nrun = 10, seed = 42))
        got <- sapply(seq_len(k), function(i) {
            s <- W[, i]; names(s) <- rownames(W)
            names(s)[order(s, decreasing = TRUE)[seq_len(min(TOP_GENES,
                                                             length(s)))]]
        })
        want <- as.matrix(read.csv(want_file, check.names = FALSE))
        # Program order is arbitrary between runs; compare the sets of lists.
        # unname() matters: apply() over the deposited matrix carries its
        # column names P1..Pk while apply() over the computed one does not,
        # and identical() compares attributes as well as content. Without it
        # every list is reported as differing at a gene overlap of 1.000,
        # which is the comparison failing rather than the data disagreeing.
        got_s <- sort(unname(apply(got, 2, paste, collapse = "|")))
        want_s <- sort(unname(apply(want, 2, function(x)
            paste(x[!is.na(x)], collapse = "|"))))
        if (identical(got_s, want_s)) {
            agree <- agree + 1
            cat(sprintf("  K=%d  identical\n", k))
        } else {
            differ <- differ + 1
            shared <- mean(mapply(function(a, b)
                length(intersect(strsplit(a, "\\|")[[1]],
                                 strsplit(b, "\\|")[[1]])) / TOP_GENES,
                got_s, want_s))
            cat(sprintf("  K=%d  DIFFERS   mean gene overlap %.3f\n", k, shared))
        }
    }
}

cat(sprintf("\n%d gene lists identical   %d differ   %d had no deposited file\n",
            agree, differ, skipped))
quit(status = if (differ > 0) 1 else 0)
