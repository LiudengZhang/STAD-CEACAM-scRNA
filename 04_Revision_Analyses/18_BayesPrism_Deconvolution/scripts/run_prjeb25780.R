#!/usr/bin/env Rscript
# BayesPrism on PRJEB25780: reference of cell-type expression profiles
# (input.type = "GEP"), cell types as states, no tumour key.
# Writes bayesprism_fractions.tsv and bayesprism_epithelial_expression.tsv.
# Usage: Rscript run_prjeb25780.R ../work/prjeb25780
suppressPackageStartupMessages({ library(BayesPrism); library(data.table) })
setwd(commandArgs(trailingOnly = TRUE)[1])
rd <- function(f) { d <- fread(f, sep = "\t", header = TRUE); r <- d[[1]]; d[, 1 := NULL]
                    m <- as.matrix(d); rownames(m) <- r; m }
sc <- rd("sc_counts.tsv"); bk <- rd("bulk_counts.tsv")
lab <- fread("sc_cell_types.tsv", sep = "\t", header = TRUE)
stopifnot(identical(lab[[1]], rownames(sc)))

p <- new.prism(reference = sc, mixture = bk, input.type = "GEP",
               cell.type.labels = lab$cell_type, cell.state.labels = lab$cell_type,
               key = NULL, outlier.cut = 0.01, outlier.fraction = 0.1)
bp <- run.prism(prism = p, n.cores = 4)
write.table(get.fraction(bp, which.theta = "final", state.or.type = "type"),
            "bayesprism_fractions.tsv", sep = "\t", quote = FALSE)
write.table(get.exp(bp = bp, state.or.type = "type", cell.name = "Epithelial cells"),
            "bayesprism_epithelial_expression.tsv", sep = "\t", quote = FALSE)
