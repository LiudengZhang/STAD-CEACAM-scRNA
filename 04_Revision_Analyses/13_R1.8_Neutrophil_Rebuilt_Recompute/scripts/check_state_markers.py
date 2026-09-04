"""
Do the eight shipped neutrophil states still hold on the corrected matrix?

The labels `C0_Neu_S100A12` ... `C7_Neu_IGKC` were assigned by a leiden run
(res 0.8) computed on the doubly normalised matrix, and each is named after the
marker that topped its cluster in that run. The manuscript's neutrophil subsets
are these eight, so they are carried across the rebuild verbatim and nothing
here relabels anything.

This is a sanity check and only that: it re-ranks markers on the sound matrix,
with the same cells in the same groups, and reports whether each state's naming
marker is still its top marker. A marker that no longer holds is a finding to
put in front of the author, not a reason to recluster.

Run:
    python check_state_markers.py [--h5ad PATH] [--out DIR]
"""

from pathlib import Path
import argparse
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

HERE = Path(__file__).resolve().parent
DEFAULT = (HERE.parents[2] / "06_Clean_Data" / "02_Rebuilt"
           / "Neutrophils_sound.h5ad")
KEY = "minor_cell_state"
N_TOP = 200


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=str(DEFAULT))
    ap.add_argument("--out", default=str(HERE.parents[0] / "outputs"))
    args = ap.parse_args()

    # Only .X and the one obs column are read. Pulling the whole object in
    # brings .raw, layers['counts'] and the two 61,167 x 61,167 graphs with it
    # for no benefit.
    import anndata as ad
    import h5py
    from anndata.experimental import read_elem
    with h5py.File(args.h5ad, "r") as f:
        X = read_elem(f["X"])
        obs = read_elem(f["obs"])[[KEY]]
        var = read_elem(f["var"])[[]]
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.obs[KEY] = a.obs[KEY].astype(str)
    print(f"{args.h5ad}\n{a.n_obs:,} cells x {a.n_vars:,} genes")

    sc.tl.rank_genes_groups(a, groupby=KEY, method="wilcoxon", use_raw=False,
                            n_genes=N_TOP)
    rows = []
    for state in sorted(a.obs[KEY].unique()):
        marker = state.rsplit("_", 1)[-1]
        d = sc.get.rank_genes_groups_df(a, group=state)
        top = d["names"].tolist()
        rank = top.index(marker) + 1 if marker in top else None
        rows.append(dict(
            state=state, n_cells=int((a.obs[KEY] == state).sum()),
            naming_marker=marker,
            marker_in_var=bool(marker in a.var_names),
            rank_of_naming_marker=rank,
            still_top_marker=bool(rank == 1),
            top5=", ".join(top[:5])))
    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "state_marker_check.csv", index=False)
    print(df.to_string(index=False))
    held = int(df["still_top_marker"].sum())
    print(f"\nnaming marker is still the top marker in {held} of {len(df)} "
          f"states (Wilcoxon, one state against the other seven, on the "
          f"singly normalised matrix)")
    print(f"wrote {out / 'state_marker_check.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
