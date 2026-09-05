"""
Shared readers for the claim ledger. Copied from
04_Revision_Analyses/16_GSEA_Metric_Sensitivity/scripts/metric_sensitivity.py
(its `load_rank` / `run_one` / gmt reader) so that the numbers this module
measures are produced by the same code path the sibling module used, and any
difference between the two is a difference in what was asked, not in how it
was read.

Everything here is read-only. Nothing outside 17_NFkB_Claim_Ledger/ is written.
"""

from pathlib import Path
import os

import numpy as np
import pandas as pd

os.environ.setdefault("PYTHONHASHSEED", "0")

HERE = Path(__file__).resolve().parent
MOD = HERE.parent
ROOT = MOD.parents[1]
OUT = MOD / "outputs"

HALLMARK = ROOT / "00_Reference" / "MSigDB_Hallmark_2020.gmt"

# `live`   - module 12's own outputs, from full_dataset.h5ad. These are the
#            tables verify_numbers.py checks the main text against.
# `sound13`- the 2026-08-31 recompute on sound per-cell-type inputs: twelve
#            types from the archive, neutrophils from module 13.
LIVE = ROOT / "04_Revision_Analyses" / "12_R1.8_DEG_Recompute" / "outputs" / "deg"
SOUND12 = (ROOT / "07_Archive"
           / "2026-08-31_deg_recompute_on_sound_per_cell_type_inputs"
           / "04_Revision_Analyses" / "12_R1.8_DEG_Recompute" / "outputs" / "deg")
NEUT = (ROOT / "04_Revision_Analyses" / "13_R1.8_Neutrophil_Rebuilt_Recompute"
        / "outputs" / "deg")

PUB_SOUND13 = (ROOT / "04_Revision_Analyses" / "13_R1.8_Neutrophil_Rebuilt_Recompute"
               / "outputs" / "nfkb_per_celltype_sound13.csv")
PUB_LIVE = (ROOT / "04_Revision_Analyses" / "13_R1.8_Neutrophil_Rebuilt_Recompute"
            / "outputs" / "nfkb_per_celltype_live.csv")
# the table verify_numbers.py actually reads
PUB_MANUSCRIPT = (ROOT / "04_Revision_Analyses" / "07_R1.8_NFkB_Specificity"
                  / "outputs" / "nfkb_per_celltype.csv")

CELLS = ["B_cells", "DC_cells", "Endothelial_cells", "Epithelial", "Fibroblast",
         "Mast_cells", "MoMac", "Neutrophils", "NK_cells", "Pericyte",
         "Plasma_cells", "TCD4_cells", "TCD8_cells"]
PHASES = ["pre", "post"]
NFKB_TERM = "TNF-alpha Signaling via NF-kB"
PUB_SEED = 42
SEEDS = [42, 1, 7, 13, 101, 2026, 777]


def read_gmt():
    sets = {}
    for line in HALLMARK.read_text().splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            sets[f[0]] = {g for g in f[2:] if g}
    return sets


GMT = read_gmt()


def deg_path(cell, phase, degset):
    if degset == "live":
        p = LIVE
    else:
        p = NEUT if cell == "Neutrophils" else SOUND12
    return p / f"{cell}_{phase}_ttest.csv"


def load_rank(cell, phase, degset):
    """Published ranking metric: logFC * -log10(P). recompute_deg.py:343."""
    d = pd.read_csv(deg_path(cell, phase, degset))
    d = d.dropna(subset=["logfoldchanges", "pvals", "scores"]).copy()
    d["metric"] = d["logfoldchanges"] * -np.log10(d["pvals"].clip(lower=1e-300))
    d = d[np.isfinite(d["metric"])]
    return (d.set_index("gene")["metric"].groupby(level=0).first()
            .sort_values(ascending=False))


def run_one(rnk, seed, perm=1000, min_size=15, max_size=500, threads=4):
    import gseapy as gp
    res = gp.prerank(rnk=rnk, gene_sets=str(HALLMARK), permutation_num=perm,
                     min_size=min_size, max_size=max_size, seed=seed,
                     no_plot=True, outdir=None, threads=threads)
    out = res.res2d.copy()
    out["Term"] = out["Term"].astype(str).str.replace(r"^.*__", "", regex=True)
    for c in ("ES", "NES", "NOM p-val", "FDR q-val"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out = out.sort_values("NES", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out
