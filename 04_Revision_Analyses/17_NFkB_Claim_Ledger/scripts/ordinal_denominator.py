"""
The printed Methods say the NF-kB prerank used `min_size = 5` (Manuscript
Methods, "Cell-Type-Specific NF-kB Pathway Enrichment"). The tables the numbers
are actually read from were produced at `min_size = 15`
(12_R1.8_DEG_Recompute/scripts/recompute_deg.py:351-352, and the same in
13_.../recompute_deg_neutrophils.py:375-376).

min_size does not change NES at all - 16_GSEA_Metric_Sensitivity measured
max |dNES| = 0.000 across min_size 5/10/25 - but it sets the denominator of
every printed ordinal. This measures, for the two printed ordinals and for the
"top-ranked in five" count, what they would be under the min_size the Methods
declare.

Both DE sets, published seed 42, everything else fixed.
"""
import sys
import numpy as np
import pandas as pd
from gsea_common import load_rank, run_one, CELLS, NFKB_TERM, OUT

rows = []
for degset in ("live", "sound13"):
    for mn in (5, 10, 15, 25):
        rank1 = []
        for cell in CELLS:
            rnk = load_rank(cell, "post", degset)
            tab = run_one(rnk, 42, min_size=mn)
            r = tab[tab["Term"] == NFKB_TERM]
            if r.empty:
                rows.append(dict(degset=degset, min_size=mn, cell=cell,
                                 n_sets=len(tab), rank=np.nan, nes=np.nan))
                continue
            r = r.iloc[0]
            rows.append(dict(degset=degset, min_size=mn, cell=cell,
                             n_sets=len(tab), rank=int(r["rank"]),
                             nes=float(r["NES"]), fdr_q=float(r["FDR q-val"])))
            if int(r["rank"]) == 1:
                rank1.append(cell)
        print(f"{degset:<8} min_size={mn:<3} top-ranked in {len(rank1)}: "
              f"{', '.join(rank1)}", flush=True)
d = pd.DataFrame(rows)
d.to_csv(OUT / "ordinal_denominator.csv", index=False)
print()
print(d[d["cell"].isin(["Pericyte", "B_cells"])].to_string(index=False))
sys.exit(0)
