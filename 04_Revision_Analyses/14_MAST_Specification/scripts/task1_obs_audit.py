"""
Task 1, part 1: is `condition` constant within `sample`, for every cell type
and phase?

Reads obs only (backed mode), so no expression matrix is touched. Sample
identifiers are never written out - only counts and the anonymised index of a
sample within its own contrast.
"""
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

from sources import CELL_SOURCES, PHASES

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

rows = []
for cell, path in CELL_SOURCES.items():
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs
    for phase, cfg in PHASES.items():
        col = cfg["column"]
        if col not in obs.columns:
            rows.append(dict(cell=cell, phase=phase, status=f"no column {col}"))
            continue
        g = obs[col].astype(str)
        m = g.isin([cfg["ref"], cfg["test"]])
        sub = obs[m]
        if len(sub) == 0:
            rows.append(dict(cell=cell, phase=phase, status="no cells"))
            continue
        cond = sub[col].astype(str)
        samp = sub["sample"].astype(str)
        tab = pd.crosstab(samp, cond)
        n_mixed = int((tab > 0).sum(axis=1).gt(1).sum())
        n_ref_s = int((tab.get(cfg["ref"], 0) > 0).sum())
        n_test_s = int((tab.get(cfg["test"], 0) > 0).sum())
        rows.append(dict(
            cell=cell, phase=phase, status="ok",
            n_cells=int(len(sub)),
            n_ref_cells=int((cond == cfg["ref"]).sum()),
            n_test_cells=int((cond == cfg["test"]).sum()),
            n_samples=int(tab.shape[0]),
            n_R_samples=n_ref_s, n_NR_samples=n_test_s,
            n_samples_with_both_conditions=n_mixed,
            condition_constant_within_sample=(n_mixed == 0),
            min_cells_per_sample=int(tab.sum(axis=1).min()),
            median_cells_per_sample=int(tab.sum(axis=1).median()),
            max_cells_per_sample=int(tab.sum(axis=1).max()),
        ))
        print(f"[{cell} {phase}] samples={tab.shape[0]} R={n_ref_s} NR={n_test_s} "
              f"mixed={n_mixed} cells={len(sub)}", flush=True)
    a.file.close()

df = pd.DataFrame(rows)
df.to_csv(OUT / "task1_sample_condition_structure.csv", index=False)
print()
print(df.to_string(index=False))
