# DESeq2 on the same summed counts and the same genes limma-voom kept, every
# contrast in one R session.
#
# Non-responder against responder, responder as the reference level, so a
# positive log2FoldChange means higher in non-responders.
#
# Standard Wald test, DESeq() defaults, no shrinkage applied to the ranking
# metric: the metric recompute_deg.py ranks on is logFC x -log10(P), and
# shrinking the fold change would change the metric rather than the test.
# results() is called with cooksCutoff and independentFiltering at their
# defaults; genes those set to NA are dropped before ranking, as they are for
# every other method here.
#
# The node this ran on was shared with a job holding 26 cores, and DESeq2 was
# taking about five minutes a contrast at one core. It is therefore shardable:
# <shard> of <nshards>, each process taking every nshards-th contrast, and a
# contrast whose output already exists is skipped so an interrupted run resumes.
# Sharding changes no result - the contrasts are independent fits.
#
# Usage: Rscript run_deseq2.R <pseudobulk_dir> <de_out_dir> [shard] [nshards]

suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
pb <- args[1]; de <- args[2]
shard <- if (length(args) >= 4) as.integer(args[3]) else 1
nshards <- if (length(args) >= 4) as.integer(args[4]) else 1
cat(sprintf("DESeq2 %s  shard %d of %d\n", packageVersion("DESeq2"), shard, nshards))

all_f <- sort(Sys.glob(file.path(pb, "*_counts.tsv.gz")))
mine <- all_f[seq_along(all_f) %% nshards == (shard - 1) %% nshards]
for (counts_f in mine) {
  base <- sub("_counts\\.tsv\\.gz$", "", basename(counts_f))
  if (file.exists(file.path(de, paste0(base, "_deseq2.csv")))) {
    cat(sprintf("DESeq2 have    %s\n", base)); flush(stdout()); next
  }
  meta_f <- file.path(pb, paste0(base, "_samples.tsv"))
  genes_f <- file.path(de, paste0(base, "_genes_kept.txt"))
  prefix <- file.path(de, base)
  if (!file.exists(genes_f)) {
    cat(sprintf("DESeq2 SKIPPED %s - no gene list, limma-voom did not run\n", base))
    next
  }

  counts <- as.matrix(read.delim(gzfile(counts_f), row.names = 1, check.names = FALSE))
  meta <- read.delim(meta_f, stringsAsFactors = FALSE)
  stopifnot(identical(colnames(counts), meta$sample_code))

  counts <- counts[readLines(genes_f), , drop = FALSE]
  mode(counts) <- "integer"

  coldata <- data.frame(group = factor(meta$group, levels = c("R", "NR")),
                        row.names = meta$sample_code)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "NR", "R"))

  out <- data.frame(gene = rownames(res),
                    logfoldchanges = res$log2FoldChange,
                    base_mean = res$baseMean,
                    stat = res$stat,
                    pvals = res$pvalue,
                    pvals_adj = res$padj,
                    stringsAsFactors = FALSE)
  write.csv(out, paste0(prefix, "_deseq2.csv"), row.names = FALSE)

  cat(sprintf("DESeq2      %-26s %6d genes, %2d samples (%d R / %d NR), %4d p NA, %5d FDR<0.05\n",
              base, nrow(out), ncol(counts),
              sum(coldata$group == "R"), sum(coldata$group == "NR"),
              sum(is.na(out$pvals)), sum(out$pvals_adj < 0.05, na.rm = TRUE)))
  flush(stdout())
}
cat("DESeq2 done\n")
