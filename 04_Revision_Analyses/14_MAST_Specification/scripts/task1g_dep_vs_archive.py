"""
How far the deposited model, re-run here, lands from the deposit's own archived
output - per contrast, gene by gene.

Same model, same cells, same genes, same MAST. The one difference is the sample
labelling (positional labels in Python's sorted() order here; the real
identifiers in R's collation order there). Twelve of the thirteen cell types
have an archived sound-run table; neutrophils do not, because that run failed
on the damaged input and the rebuilt object post-dates it.
"""
from pathlib import Path

import numpy as np
import pandas as pd

MOD = Path(__file__).resolve().parents[1]
OUT = MOD / "outputs"
ARCH = (MOD.parents[1] / "07_Archive"
        / "2026-08-31_deg_recompute_on_sound_per_cell_type_inputs"
        / "04_Revision_Analyses" / "12_R1.8_DEG_Recompute" / "outputs" / "deg")
CELLS = ["B_cells", "DC_cells", "Endothelial_cells", "Epithelial", "Fibroblast",
         "Mast_cells", "MoMac", "Neutrophils", "NK_cells", "Pericyte",
         "Plasma_cells", "TCD4_cells", "TCD8_cells"]

rows = []
for cell in CELLS:
    for phase in ("pre", "post"):
        for method, mine_spec in (("mast", "dep"), ("ttest", "T")):
            a = ARCH / f"{cell}_{phase}_{method}.csv"
            b = OUT / "deg" / f"{cell}_{phase}_{mine_spec}.csv"
            if not a.exists() or not b.exists():
                rows.append(dict(cell=cell, phase=phase, method=method,
                                 status="no archived table" if not a.exists()
                                 else "not run here"))
                continue
            da = pd.read_csv(a).set_index("gene")
            db = pd.read_csv(b).set_index("gene")
            c = da.index.intersection(db.index)
            dl = (da.loc[c, "logfoldchanges"] - db.loc[c, "logfoldchanges"]).abs()
            dp = (np.log10(da.loc[c, "pvals"].clip(lower=1e-300))
                  - np.log10(db.loc[c, "pvals"].clip(lower=1e-300))).abs()
            rows.append(dict(cell=cell, phase=phase, method=method,
                             status="ok", n_archive=len(da), n_here=len(db),
                             n_common=len(c),
                             max_abs_dlogFC=float(dl.max()),
                             median_abs_dlogFC=float(dl.median()),
                             max_abs_dlog10P=float(dp.max()),
                             reproduces=bool(dl.max() < 1e-9)))
df = pd.DataFrame(rows)
df.to_csv(OUT / "task1_dep_vs_archive_by_contrast.csv", index=False)
ok = df[df.status == "ok"]
print(ok.to_string(index=False))
print()
print(ok.groupby(["method", "phase"])["reproduces"].agg(["sum", "count"]).to_string())
