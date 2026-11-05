#!/usr/bin/env Rscript
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
# Which nrun produced the deposited NMF programs?
#
# 02_run_nmf_per_sample.R hard-codes a path to a config that is not the one the
# standardised pipeline carries, and the two disagree on exactly one parameter:
#
#     01.1_meta_program_epithelial/00_Config/config.yaml   n_iter: 10
#     _Archived/999_Panels/.../3CA_MetaProgram/...yaml     n_iter: 100
#
# Both set random_seed: 42 and method: brunet, so the seed was never the open
# question - nrun was. nrun changes the answer: NMF runs the factorisation that
# many times and keeps the best fit. The script also carries a tryCatch that
# silently falls back to nrun=10 on any error, so a failure would look like a
# success at the wrong setting.
#
# This settles it by running one sample both ways at seed 42 and comparing the
# basis matrix against the deposited one, rather than by reading the code.

suppressPackageStartupMessages({library(NMF); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
sample <- if (length(args) >= 1) args[1] else "inLN_1221"
k <- if (length(args) >= 2) as.integer(args[2]) else 4L

dep <- file.path("/path/to/home/project",
                 "submission-tree/02_Preparation_for_Panels/NMF/per_sample_nmf", sample)

# The preprocessing below is copied from 02_run_nmf_per_sample.R unchanged. If
# it drifts, the comparison stops testing nrun and starts testing the copy.
counts <- read.csv(file.path(dep, paste0(sample, "_counts.csv")),
                   row.names = 1, check.names = FALSE)
m <- t(as.matrix(counts))
m[m < 0] <- 0
m <- m[rowSums(m) > 0, , drop = FALSE]
cs <- colSums(m); cs[cs == 0] <- 1
m <- log2(sweep(m, 2, cs, "/") * 1e6 + 1)
rm_ <- rowMeans(m)
m <- m[!is.na(rm_) & rm_ > 0, , drop = FALSE]
cat(sprintf("%s: %d genes x %d cells, K=%d\n", sample, nrow(m), ncol(m), k))

deposited <- as.matrix(read.csv(
    file.path(dep, sprintf("%s_nmf_k%d_W.csv", sample, k)), row.names = 1))
cat(sprintf("deposited W: %d x %d\n", nrow(deposited), ncol(deposited)))

for (nrun in c(10, 100)) {
    t0 <- Sys.time()
    set.seed(42)
    res <- nmf(m, rank = k, method = "brunet", nrun = nrun, seed = 42)
    W <- basis(res)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (!identical(dim(W), dim(deposited)) ||
        !identical(rownames(W), rownames(deposited))) {
        cat(sprintf("  nrun=%-3d  shape or gene order differs; not comparable\n",
                    nrun))
        next
    }
    # Column order is arbitrary between runs, so compare the best matching of
    # programs rather than position for position.
    cc <- cor(W, deposited)
    best <- apply(cc, 2, max)
    cat(sprintf("  nrun=%-3d  max |W - deposited| = %.3e   per-program r = %s   %.0fs\n",
                nrun, max(abs(W - deposited)),
                paste(sprintf("%.4f", best), collapse = " "), secs))
}
