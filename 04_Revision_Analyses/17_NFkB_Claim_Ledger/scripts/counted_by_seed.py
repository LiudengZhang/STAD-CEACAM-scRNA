"""
The counted claims, recomputed at every seed. This is the direct measurement of
how stable an ordinal or a count is when nothing but the permutation draw moves.

Reads outputs/signal_runs.csv only.
"""
import sys
import pandas as pd
from gsea_common import OUT

runs = pd.read_csv(OUT / "signal_runs.csv")
rows = []
for (degset, phase, seed), g in runs.groupby(["degset", "phase", "seed"]):
    pos = g["nfkb_nes"] > 0
    rows.append(dict(
        degset=degset, phase=phase, seed=seed,
        n_positive=int(pos.sum()),
        n_positive_q05=int((pos & (g["nfkb_fdr_q"] < 0.05)).sum()),
        n_positive_q25=int((pos & (g["nfkb_fdr_q"] < 0.25)).sum()),
        n_rank1=int((g["nfkb_rank"] == 1).sum()),
        types_rank1="; ".join(sorted(g.loc[g["nfkb_rank"] == 1, "cell"])),
        n_negative=int((~pos).sum()),
        pericyte_rank=int(g.loc[g["cell"] == "Pericyte", "nfkb_rank"].iloc[0]),
        pericyte_nsets=int(g.loc[g["cell"] == "Pericyte", "n_sets_tested"].iloc[0]),
        bcells_rank=int(g.loc[g["cell"] == "B_cells", "nfkb_rank"].iloc[0]),
        bcells_nsets=int(g.loc[g["cell"] == "B_cells", "n_sets_tested"].iloc[0]),
        momac_nes=round(float(g.loc[g["cell"] == "MoMac", "nfkb_nes"].iloc[0]), 4),
        momac_q=round(float(g.loc[g["cell"] == "MoMac", "nfkb_fdr_q"].iloc[0]), 4),
        momac_rank=int(g.loc[g["cell"] == "MoMac", "nfkb_rank"].iloc[0]),
        epithelial_nes=round(float(g.loc[g["cell"] == "Epithelial", "nfkb_nes"].iloc[0]), 4),
        epithelial_rank=int(g.loc[g["cell"] == "Epithelial", "nfkb_rank"].iloc[0]),
        strongest_post=g.loc[g["nfkb_nes"].idxmax(), "cell"]))
d = pd.DataFrame(rows).sort_values(["degset", "phase", "seed"])
d.to_csv(OUT / "counted_claims_by_seed.csv", index=False)
pd.set_option("display.width", 250)
print(d.to_string(index=False))
sys.exit(0)
