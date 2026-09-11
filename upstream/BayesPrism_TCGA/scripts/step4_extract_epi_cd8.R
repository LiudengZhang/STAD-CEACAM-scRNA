#!/usr/bin/env Rscript
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# Extract epithelial-specific and CD8+ T cell-specific expression from BayesPrism result

suppressPackageStartupMessages({
    library(BayesPrism)
})

# Migrated 2026-09-03 on the author's ruling. WORK_DIR was
# submission-tree/temp-workspace/02162026_BayesPrism, a temp workspace that can be
# cleared at any time, and it served as both the input and the output
# directory - so a re-run wrote its own output over the very file it was being
# compared against. The posterior now comes from the pipeline's declared home,
# paths.py TCGA_BAYESPRISM_DIR; the md5 was checked on both sides of the copy
# and is unchanged (52ed935c08057299e1688712a4256517), and the temp workspace
# was copied from, never emptied. Writes land where the script is invoked.
IN_DIR <- "/path/to/submission-tree/02_Preparation_for_Panels/BayesPrism_TCGA"
setwd(dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))))

cat("Loading BayesPrism result...\n")
bp_result <- readRDS(file.path(IN_DIR, "tcga_bayesprism_result.rds"))

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
