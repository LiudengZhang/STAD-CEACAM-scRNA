"""
Collect every specification's per-contrast summary into one table, add the
TNF-alpha/NF-kB row from each Hallmark result, and tabulate the agreement
between specifications gene by gene.
"""
from pathlib import Path
import itertools
import sys

import numpy as np
import pandas as pd
from scipy import stats

MOD = Path(__file__).resolve().parents[1]
OUT = MOD / "outputs"
CELLS = ["B_cells", "DC_cells", "Endothelial_cells", "Epithelial", "Fibroblast",
         "Mast_cells", "MoMac", "Neutrophils", "NK_cells", "Pericyte",
         "Plasma_cells", "TCD4_cells", "TCD8_cells"]
SPECS = ["dep", "A0", "A", "B", "T"]

rows = [pd.read_csv(f) for f in sorted((OUT / "summaries").glob("*.csv"))]
if not rows:
    sys.exit("no summaries yet")
summ = pd.concat(rows, ignore_index=True)
summ = summ.sort_values(["cell", "phase", "spec"])
summ.to_csv(OUT / "all_specifications_summary.csv", index=False)

# --- the NF-kB table the brief asks for -------------------------------------
nf = summ[["cell", "phase", "spec", "formula", "n_ref_cells", "n_test_cells",
           "n_samples", "n_genes_tested", "nfkb_nes", "nfkb_nom_p",
           "nfkb_fdr_q", "nfkb_rank", "n_sets"]].copy()
nf.to_csv(OUT / "nfkb_by_specification.csv", index=False)

wide = nf.pivot_table(index=["cell", "phase"], columns="spec",
                      values="nfkb_nes")
wide.columns = [f"NES_{c}" for c in wide.columns]
wq = nf.pivot_table(index=["cell", "phase"], columns="spec", values="nfkb_fdr_q")
wq.columns = [f"FDRq_{c}" for c in wq.columns]
w = wide.join(wq).reset_index()
w.to_csv(OUT / "nfkb_nes_wide.csv", index=False)
print(w.to_string(index=False))

# --- gene-level concordance between specifications --------------------------
con = []
for cell in CELLS:
    for phase in ("pre", "post"):
        d = {}
        for s in SPECS:
            f = OUT / "deg" / f"{cell}_{phase}_{s}.csv"
            if f.exists():
                t = pd.read_csv(f)
                t = t.dropna(subset=["logfoldchanges", "pvals"])
                t["metric"] = (t["logfoldchanges"]
                               * -np.log10(t["pvals"].clip(lower=1e-300)))
                d[s] = t.set_index("gene")
        for a, b in itertools.combinations([s for s in SPECS if s in d], 2):
            ga = d[a]["logfoldchanges"]
            gb = d[b]["logfoldchanges"]
            common = ga.index.intersection(gb.index)
            if len(common) < 20:
                continue
            va, vb = ga.loc[common], gb.loc[common]
            ma = d[a]["metric"].loc[common]
            mb = d[b]["metric"].loc[common]
            con.append(dict(
                cell=cell, phase=phase, spec_a=a, spec_b=b, n_genes=len(common),
                spearman_logFC=float(stats.spearmanr(va, vb).statistic),
                pearson_logFC=float(stats.pearsonr(va, vb).statistic),
                spearman_rankmetric=float(stats.spearmanr(ma, mb).statistic),
                pct_same_sign_logFC=float(100 * (np.sign(va) == np.sign(vb)).mean()),
            ))
if con:
    c = pd.DataFrame(con)
    c.to_csv(OUT / "gene_level_concordance.csv", index=False)
    print()
    print(c.groupby(["spec_a", "spec_b"])[
        ["spearman_logFC", "spearman_rankmetric", "pct_same_sign_logFC"]
    ].median().round(3).to_string())

# --- NES sign agreement, per pair of specifications -------------------------
agr = []
for a, b in itertools.combinations(SPECS, 2):
    ca, cb = f"NES_{a}", f"NES_{b}"
    if ca not in w.columns or cb not in w.columns:
        continue
    d = w[["cell", "phase", ca, cb]].dropna()
    if not len(d):
        continue
    agr.append(dict(spec_a=a, spec_b=b, n_contrasts=len(d),
                    n_same_sign=int((np.sign(d[ca]) == np.sign(d[cb])).sum()),
                    spearman_NES=float(stats.spearmanr(d[ca], d[cb]).statistic)))
if agr:
    ag = pd.DataFrame(agr)
    ag.to_csv(OUT / "nfkb_nes_agreement.csv", index=False)
    print()
    print(ag.to_string(index=False))
