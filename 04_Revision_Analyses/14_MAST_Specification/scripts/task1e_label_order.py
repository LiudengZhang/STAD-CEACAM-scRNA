"""
Why the deposited MAST result changes when the samples are renamed.

Running the deposited model through this module's harness reproduces the
archived 2026-08-31 sound run to 1e-15 for every `post` contrast and NOT for
any `pre` contrast, on a matrix that is identical (the Welch t-test on the same
cells reproduces exactly, and the cell and gene counts match).

The only thing this harness changes is the sample labels: it replaces the
specimen identifiers with positional labels S01, S02, ... assigned in Python's
sorted() order. That is meant to be a no-op. It is not, because the factor
level order decides which sample dummy MAST drops as "never estimible", and
Python's sorted() (codepoint order) and R's factor() (locale collation) do not
always agree on identifiers containing punctuation.

This script measures that, per contrast, without writing any identifier: it
reports only whether the two orders agree and how many positions differ.
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from sources import CELL_SOURCES, PHASES

OUT = Path(__file__).resolve().parents[1] / "outputs"

import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter

rows = []
for cell, path in CELL_SOURCES.items():
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs
    for phase, cfg in PHASES.items():
        g = obs[cfg["column"]].astype(str)
        m = g.isin([cfg["ref"], cfg["test"]])
        if m.sum() == 0:
            continue
        samples = obs.loc[m, "sample"].astype(str).unique()
        py_order = sorted(samples)
        with localconverter(ro.default_converter + pandas2ri.converter
                            + numpy2ri.converter):
            ro.r.assign("s", np.asarray([str(x) for x in samples], dtype=object))
        r_order = list(ro.r('levels(factor(as.character(s)))'))
        n_diff = int(sum(1 for x, y in zip(py_order, r_order) if x != y))
        rows.append(dict(cell=cell, phase=phase, n_samples=len(samples),
                         python_sorted_equals_R_factor_levels=(py_order == r_order),
                         n_positions_differing=n_diff))
        print(f"[{cell} {phase}] n={len(samples)} orders agree="
              f"{py_order == r_order} positions differing={n_diff}", flush=True)
    a.file.close()

df = pd.DataFrame(rows)
df.to_csv(OUT / "task1_label_order_agreement.csv", index=False)
print()
print(df.to_string(index=False))
print("\nby phase:")
print(df.groupby("phase")["python_sorted_equals_R_factor_levels"]
      .agg(["sum", "count"]).to_string())
print("\ndone", flush=True)
