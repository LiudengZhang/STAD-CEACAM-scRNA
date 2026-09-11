#!/usr/bin/env Rscript
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# Extract epithelial-specific and CD8+ T cell-specific expression from BayesPrism result

suppressPackageStartupMessages({
    library(BayesPrism)
})

# One change from step4_extract_epi_cd8.R, and it is here to stop the script
# destroying its own control. As found, WORK_DIR is an absolute path into the
# deposited working directory and every write below lands there - so a re-run
# meant to be compared against those files would have overwritten them first.
# The script runs where it is invoked instead.
setwd(dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))))

cat("Loading BayesPrism result...\n")
bp_result <- readRDS("tcga_bayesprism_result.rds")

# 1. Epithelial-specific expression
cat("Extracting Epithelial cells expression...\n")
epi_expr <- get.exp(bp = bp_result, state.or.type = "type", cell.name = "Epithelial cells")
cat(sprintf("  Epithelial: %d samples x %d genes\n", nrow(epi_expr), ncol(epi_expr)))
write.table(epi_expr, "tcga_bayesprism_epithelial_expression.tsv", sep="\t", quote=FALSE)

# Check key genes
key_genes <- c("CEACAM5", "CEACAM6", "CD274")
for (g in key_genes) {
    if (g %in% colnames(epi_expr)) {
        cat(sprintf("  %s: mean=%.4f, median=%.4f\n", g, mean(epi_expr[, g]), median(epi_expr[, g])))
    } else {
        cat(sprintf("  %s: NOT FOUND\n", g))
    }
}

# 2. CD8+ T cell-specific expression
cat("\nExtracting CD8+ T cells expression...\n")
cd8_expr <- get.exp(bp = bp_result, state.or.type = "type", cell.name = "CD8+ T cells")
cat(sprintf("  CD8+ T: %d samples x %d genes\n", nrow(cd8_expr), ncol(cd8_expr)))
write.table(cd8_expr, "tcga_bayesprism_cd8_expression.tsv", sep="\t", quote=FALSE)

# Check Tex markers
tex_genes <- c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1")
cat("Tex markers in CD8 expression:\n")
for (g in tex_genes) {
    if (g %in% colnames(cd8_expr)) {
        cat(sprintf("  %s: mean=%.4f, median=%.4f\n", g, mean(cd8_expr[, g]), median(cd8_expr[, g])))
    } else {
        cat(sprintf("  %s: NOT FOUND\n", g))
    }
}

cat("\nDone!\n")
