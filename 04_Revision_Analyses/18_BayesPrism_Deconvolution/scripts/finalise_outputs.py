"""Collect the BayesPrism results into outputs/.

TCGA epithelial expression is written as per-sample counts per million, with
the CD8A and CD8B counts and library size of each tumour (Figure 3E).
Run: python finalise_outputs.py   (after run_tcga.R and run_prjeb25780.R)
"""
from pathlib import Path
import shutil

import pandas as pd

MOD = Path(__file__).resolve().parents[1]
WORK, OUT = MOD / "work", MOD / "outputs"

for d in ("tcga", "prjeb25780"):
    (OUT / d).mkdir(parents=True, exist_ok=True)
epi = pd.read_csv(WORK / "tcga" / "epithelial_expression.tsv", sep="\t", index_col=0)
(epi.div(epi.sum(axis=1), axis=0) * 1e6).to_csv(OUT / "tcga" / "epithelial_cpm.tsv", sep="\t")
shutil.copy2(WORK / "tcga" / "fractions.tsv", OUT / "tcga" / "fractions.tsv")
bulk = pd.read_csv(WORK / "tcga" / "bulk_counts.tsv", sep="\t", index_col=0)
pd.DataFrame({"library_size": bulk.sum(axis=1), "CD8A": bulk["CD8A"], "CD8B": bulk["CD8B"]}).to_csv(
    OUT / "tcga" / "tumour_cd8_counts.tsv", sep="\t")
for name in ("bayesprism_epithelial_expression.tsv", "bayesprism_fractions.tsv"):
    shutil.copy2(WORK / "prjeb25780" / name, OUT / "prjeb25780" / name)
