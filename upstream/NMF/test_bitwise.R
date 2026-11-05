#!/usr/bin/env Rscript
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
# Is the NMF step bit-identical run to run, or only identical to 6e-12?
#
# test_nrun.R settled which nrun made the deposited programs: at seed 42,
# method brunet, nrun 10, the basis matrix reproduces to max |dW| = 6.4e-12
# with every program correlating at 1.0000, while nrun 100 gives a different
# answer entirely. 6e-12 is floating-point noise, not disagreement - but noise
# is what a multi-threaded BLAS produces by reordering a sum, and the round's
# standard is byte equality.
#
# So this runs the same factorisation twice in one process with the BLAS pinned
# to a single thread. If the two agree exactly, the pipeline is bit-reproducible
# once OMP_NUM_THREADS is fixed, and run.sh fixes it.

suppressPackageStartupMessages({library(NMF)})

args <- commandArgs(trailingOnly = TRUE)
sample <- if (length(args) >= 1) args[1] else "inLN_1221"
k <- if (length(args) >= 2) as.integer(args[2]) else 4L

dep <- file.path("/path/to/home/project",
                 "submission-tree/02_Preparation_for_Panels/NMF/per_sample_nmf", sample)

counts <- read.csv(file.path(dep, paste0(sample, "_counts.csv")),
                   row.names = 1, check.names = FALSE)
m <- t(as.matrix(counts))
m[m < 0] <- 0
m <- m[rowSums(m) > 0, , drop = FALSE]
cs <- colSums(m); cs[cs == 0] <- 1
m <- log2(sweep(m, 2, cs, "/") * 1e6 + 1)
rm_ <- rowMeans(m)
m <- m[!is.na(rm_) & rm_ > 0, , drop = FALSE]

run_once <- function() {
    set.seed(42)
    basis(nmf(m, rank = k, method = "brunet", nrun = 10, seed = 42))
}

W1 <- run_once()
W2 <- run_once()

cat(sprintf("threads: OMP_NUM_THREADS=%s\n", Sys.getenv("OMP_NUM_THREADS", "unset")))
cat(sprintf("%s K=%d: max |W1 - W2| = %.3e   identical = %s\n",
            sample, k, max(abs(W1 - W2)), identical(W1, W2)))

deposited <- as.matrix(read.csv(
    file.path(dep, sprintf("%s_nmf_k%d_W.csv", sample, k)), row.names = 1))
cat(sprintf("vs deposited: max |dW| = %.3e\n", max(abs(W1 - deposited))))
