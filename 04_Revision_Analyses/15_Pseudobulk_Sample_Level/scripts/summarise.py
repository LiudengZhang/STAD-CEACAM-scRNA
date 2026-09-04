"""
The sample-level result beside the published per-cell one.

The Hallmark tables are parsed with 07_R1.8_NFkB_Specificity's own `load_gsea`,
pointed at this module's gsea/ directory, so the NES, nominal P, FDR q and rank
are computed exactly as the shipped analysis computes them and a difference
here is a difference in the statistics, not in the reading of them. The counted
claims are 13_R1.8's own `claims()` for the same reason.

Published comparator: 13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/
nfkb_per_celltype_sound13.csv, method `ttest` - the per-cell Welch t-test on
the sound per-cell-type inputs with the rebuilt neutrophils.

Outputs
    nfkb_pseudobulk.csv          NF-kB row per cell type, phase and method
    nfkb_comparison.csv          the thirteen-by-two table, all three methods
    method_concordance.csv       sign and significance agreement between the
                                 three, contrast by contrast
    counted_claims.csv           the counted claims under each method
    n_samples.csv                samples and cells behind every contrast
    de_gene_counts.csv           genes tested and genes at FDR<0.05, per
                                 contrast and method

Run: python summarise.py
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
OUT = MODULE / "outputs"
NEW = MODULE.parent

sys.path.insert(0, str(NEW / "07_R1.8_NFkB_Specificity" / "scripts"))
sys.path.insert(0, str(NEW / "13_R1.8_Neutrophil_Rebuilt_Recompute" / "scripts"))
import nfkb_specificity as nf                                        # noqa: E402
import count_nfkb_claims as cc                                       # noqa: E402

PUBLISHED = (NEW / "13_R1.8_Neutrophil_Rebuilt_Recompute" / "outputs"
             / "nfkb_per_celltype_sound13.csv")
CELLS = ["B_cells", "DC_cells", "Endothelial_cells", "Epithelial",
         "Fibroblast", "Mast_cells", "MoMac", "Neutrophils", "NK_cells",
         "Pericyte", "Plasma_cells", "TCD4_cells", "TCD8_cells"]


def main():
    nf.RECOMPUTE_GSEA = OUT / "gsea"
    pb = nf.load_gsea()
    pb.to_csv(OUT / "nfkb_pseudobulk.csv", index=False)

    pub = pd.read_csv(PUBLISHED)
    pub = pub[pub["method"] == "ttest"].copy()
    pub["method"] = "ttest_percell"

    both = pd.concat([pb, pub], ignore_index=True)
    build = pd.read_csv(OUT / "build_report.csv")

    rows = []
    for cell in CELLS:
        for phase in ("pre", "post"):
            b = build[(build["cell"] == cell) & (build["phase"] == phase)]
            r = dict(cell_type=cell, phase=phase)
            if len(b):
                b = b.iloc[0]
                r.update(n_samples=int(b["n_samples_kept"]),
                         n_R=int(b["n_R"]), n_NR=int(b["n_NR"]),
                         n_cells=int(b["n_cells"]),
                         median_cells_per_sample=b["median_cells_per_sample"],
                         samples_dropped=int(b["n_samples_dropped"]),
                         testable=bool(b["testable"]))
            for m in ("limma", "deseq2", "ttest_percell"):
                s = both[(both["method"] == m) & (both["cell_type"] == cell)
                         & (both["phase"] == phase)]
                if len(s):
                    s = s.iloc[0]
                    r[f"{m}_nes"] = round(float(s["nes"]), 4)
                    r[f"{m}_p"] = float(s["nom_p"])
                    r[f"{m}_q"] = float(s["fdr_q"])
                    r[f"{m}_rank"] = int(s["rank"])
                    r[f"{m}_nsets"] = int(s["n_sets"])
            rows.append(r)
    comp = pd.DataFrame(rows)
    comp.to_csv(OUT / "nfkb_comparison.csv", index=False)

    claims = []
    for m in ("limma", "deseq2", "ttest_percell"):
        cc.PRIMARY = m
        claims.append(cc.claims(both, m))
    cl = pd.DataFrame(claims)
    cl.to_csv(OUT / "counted_claims.csv", index=False)

    # Agreement, contrast by contrast. Two methods that share an input and a
    # gene universe but not a model either agree, in which case the enrichment
    # is not an artefact of either, or they do not.
    conc = []
    for a, b in (("limma", "deseq2"), ("limma", "ttest_percell"),
                 ("deseq2", "ttest_percell")):
        for phase in ("pre", "post", "both"):
            s_ = comp if phase == "both" else comp[comp["phase"] == phase]
            s_ = s_.dropna(subset=[f"{a}_nes", f"{b}_nes"])
            same = np.sign(s_[f"{a}_nes"]) == np.sign(s_[f"{b}_nes"])
            conc.append(dict(
                method_a=a, method_b=b, phase=phase, n_contrasts=len(s_),
                same_sign=int(same.sum()),
                both_pos_q05=int(((s_[f"{a}_nes"] > 0) & (s_[f"{a}_q"] < 0.05)
                                  & (s_[f"{b}_nes"] > 0)
                                  & (s_[f"{b}_q"] < 0.05)).sum()),
                nes_pearson=round(float(s_[f"{a}_nes"].corr(s_[f"{b}_nes"])), 4),
                nes_spearman=round(float(s_[f"{a}_nes"].corr(
                    s_[f"{b}_nes"], method="spearman")), 4),
                disagreeing=", ".join(
                    s_.loc[~same, "cell_type"] + " " + s_.loc[~same, "phase"])))
    pd.DataFrame(conc).to_csv(OUT / "method_concordance.csv", index=False)

    build[["cell", "phase", "n_cells", "n_samples_all", "n_samples_kept",
           "n_samples_dropped", "n_R", "n_NR", "median_cells_per_sample",
           "min_cells_per_sample", "max_cells_per_sample", "n_genes_nonzero",
           "counts_source", "testable", "status"]].to_csv(
        OUT / "n_samples.csv", index=False)

    # How much gene-level signal a sample-level test finds at all. This is the
    # number the GSEA result has to be read against: a prerank enrichment can
    # be strong on a ranking whose individual genes are all non-significant.
    deg = []
    for f in sorted((OUT / "de").glob("*_limma.csv")):
        base = f.name[:-len("_limma.csv")]
        cell, phase = base.rsplit("_", 1)
        row = dict(cell_type=cell, phase=phase)
        for m in ("limma", "deseq2"):
            g = OUT / "de" / f"{base}_{m}.csv"
            if not g.exists():
                continue
            d = pd.read_csv(g)
            row[f"{m}_n_genes"] = len(d)
            row[f"{m}_n_p_na"] = int(d["pvals"].isna().sum())
            row[f"{m}_n_fdr05"] = int((d["pvals_adj"] < 0.05).sum())
            row[f"{m}_n_fdr10"] = int((d["pvals_adj"] < 0.10).sum())
        deg.append(row)
    pd.DataFrame(deg).to_csv(OUT / "de_gene_counts.csv", index=False)

    pd.set_option("display.width", 250)
    show = ["cell_type", "phase", "n_samples", "n_R", "n_NR",
            "limma_nes", "limma_q", "limma_rank",
            "deseq2_nes", "deseq2_q", "deseq2_rank",
            "ttest_percell_nes", "ttest_percell_q", "ttest_percell_rank"]
    print(comp[[c for c in show if c in comp]].to_string(index=False))
    print()
    print(cl.to_string(index=False))
    print()
    print(pd.DataFrame(conc).to_string(index=False))


if __name__ == "__main__":
    main()
