# limma-voom on the summed counts, every contrast in one R session.
#
# Non-responder against responder, responder as the reference, so a positive
# logFC means higher in non-responders - the direction recompute_deg.py reports
# and the direction the manuscript reads.
#
# The gene filter is edgeR::filterByExpr on the pseudobulk library sizes, which
# is the filter limma and edgeR document for this design. The kept gene list is
# written out because run_deseq2.R reads it: giving both methods the same gene
# universe is what makes their GSEA ranks comparable to each other, and DESeq2
# is run in a different environment that has no edgeR.
#
# A contrast whose table already exists is skipped, so a killed run resumes
# where it stopped. The fit is deterministic, so a re-run reproduces it.
#
# Usage: Rscript run_limma_voom.R <pseudobulk_dir> <de_out_dir>

suppressPackageStartupMessages({library(limma); library(edgeR)})

args <- commandArgs(trailingOnly = TRUE)
pb <- args[1]; de <- args[2]
dir.create(de, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("limma %s  edgeR %s\n", packageVersion("limma"), packageVersion("edgeR")))

for (counts_f in sort(Sys.glob(file.path(pb, "*_counts.tsv.gz")))) {
  base <- sub("_counts\\.tsv\\.gz$", "", basename(counts_f))
  meta_f <- file.path(pb, paste0(base, "_samples.tsv"))
  prefix <- file.path(de, base)
  if (file.exists(paste0(prefix, "_limma.csv")) &&
      file.exists(paste0(prefix, "_genes_kept.txt"))) {
    cat(sprintf("limma-voom have %s\n", base)); flush(stdout()); next
  }

  counts <- as.matrix(read.delim(gzfile(counts_f), row.names = 1, check.names = FALSE))
  meta <- read.delim(meta_f, stringsAsFactors = FALSE)
  stopifnot(identical(colnames(counts), meta$sample_code))

  group <- factor(meta$group, levels = c("R", "NR"))
  y <- DGEList(counts = counts, group = group)

  keep <- filterByExpr(y, group = group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)                       # TMM

  design <- model.matrix(~ group)               # coef 2 is NR - R
  v <- voom(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  out <- data.frame(gene = rownames(tt),
                    logfoldchanges = tt$logFC,
                    ave_expr = tt$AveExpr,
                    t = tt$t,
                    pvals = tt$P.Value,
                    pvals_adj = tt$adj.P.Val,
                    stringsAsFactors = FALSE)
  write.csv(out, paste0(prefix, "_limma.csv"), row.names = FALSE)
  writeLines(rownames(y), paste0(prefix, "_genes_kept.txt"))

  cat(sprintf("limma-voom  %-26s %6d of %6d genes kept, %2d samples (%d R / %d NR), %5d FDR<0.05\n",
              base, nrow(y), length(keep), ncol(y),
              sum(group == "R"), sum(group == "NR"), sum(out$pvals_adj < 0.05)))
  flush(stdout())
}
cat("limma-voom done\n")
