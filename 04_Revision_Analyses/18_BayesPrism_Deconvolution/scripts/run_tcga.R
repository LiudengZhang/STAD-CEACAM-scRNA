#!/usr/bin/env Rscript
# BayesPrism on TCGA-STAD: raw-count reference with cell states, gene cleanup,
# protein-coding marker genes, epithelial cells as the tumour component.
# Writes fractions.tsv and epithelial_expression.tsv into the work directory.
# Usage: Rscript run_tcga.R ../work/tcga
suppressPackageStartupMessages({ library(BayesPrism); library(data.table) })
setwd(commandArgs(trailingOnly = TRUE)[1])
rd <- function(f) { d <- fread(f, sep = "\t", header = TRUE); r <- d[[1]]; d[, 1 := NULL]
                    m <- as.matrix(d); rownames(m) <- r; m }
sc <- rd("sc_counts.tsv"); bk <- rd("bulk_counts.tsv")
lab <- fread("../sc_cell_types.tsv", sep = "\t", header = TRUE)
stopifnot(identical(lab[[1]], rownames(sc)))
# Kept in the reference whether or not marker selection picks them.
KEEP <- c("CEACAM5", "CEACAM6", "CD274", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1",
          "CXCL13", "LAYN", "CD38", "BATF", "IRF4", "PRDM1", "TOX2", "ITGAE", "NR4A1", "NR4A2", "NR4A3")

sc <- cleanup.genes(input = sc, input.type = "count.matrix", species = "hs",
                    gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"), exp.cells = 5)
sc <- select.gene.type(sc, gene.type = "protein_coding")
sc <- sc[, colnames(sc) %in% colnames(bk)]
st <- get.exp.stat(sc.dat = sc, cell.type.labels = lab$cell_type, cell.state.labels = lab$cell_state,
                   pseudo.count = 0.1, cell.count.cutoff = 50, n.cores = 8)
mk <- select.marker(sc.dat = sc, stat = st, pval.max = 0.01, lfc.min = 0.1)
sc <- sc[, union(colnames(mk), intersect(KEEP, colnames(sc)))]

p <- new.prism(reference = sc, mixture = bk, input.type = "count.matrix",
               cell.type.labels = lab$cell_type, cell.state.labels = lab$cell_state,
               key = "Epithelial cells", outlier.cut = 0.01, outlier.fraction = 0.1)
bp <- run.prism(prism = p, n.cores = 8)
write.table(get.fraction(bp, which.theta = "final", state.or.type = "type"),
            "fractions.tsv", sep = "\t", quote = FALSE)
write.table(get.exp(bp = bp, state.or.type = "type", cell.name = "Epithelial cells"),
            "epithelial_expression.tsv", sep = "\t", quote = FALSE)
