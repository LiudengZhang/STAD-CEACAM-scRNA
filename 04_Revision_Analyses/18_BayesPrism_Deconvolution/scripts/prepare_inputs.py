"""Inputs for the BayesPrism runs, from the deposited single-cell matrix.

Reference: 6,834 stomach cells (up to 500 per cell type). TCGA-STAD: their raw counts
with cell states (states under 20 cells merged per type), and the unstranded
STAR counts of each patient's primary tumour. PRJEB25780: their expression on
the 5,020 genes in outputs/prjeb25780/reference_genes.txt, and the cohort's expression
matrix (fetched; see TIGER_EXPR) for the 45 runs in samples.txt.

Writes ../work/. Run: python prepare_inputs.py
"""
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
import paths  # noqa: E402

WORK = Path(__file__).resolve().parents[1] / "work"
H5AD = paths.FULL_DATASET_DEPOSIT_H5AD


def h5_index(group):
    key = group.attrs["_index"]
    key = key.decode() if isinstance(key, bytes) else key
    return pd.Index([x.decode() if isinstance(x, bytes) else x for x in group[key][:]])


def read_csr(f, key, rows, n_cols):
    """Dense rows of a CSR matrix in an open h5ad, in the order given."""
    g = f[key]
    ip = g["indptr"][:]
    out = np.zeros((len(rows), n_cols), dtype=np.float64)
    for i, r in enumerate(rows):
        out[i, g["indices"][ip[r]:ip[r + 1]]] = g["data"][ip[r]:ip[r + 1]]
    return out


def cell_labels():
    ct = pd.read_csv(paths.TCGA_BAYESPRISM_DIR / "sc_cell_types.tsv", sep="\t", index_col=0)
    n = ct.groupby("cell_state").size()
    small = ct.cell_state.map(n) < 20
    ct.loc[small, "cell_state"] = ct.loc[small, "cell_type"] + "_other"
    # Hepatocytes carry no subcluster: one state, named after the type.
    ct["cell_state"] = ct.cell_state.fillna(ct.cell_type)
    return ct


def tcga_tumour_counts():
    """Unstranded counts, one primary-tumour aliquot per patient; duplicate
    symbols keep the highest-mean gene."""
    raw = Path(paths.TCGA_RAW_DIR)
    meta = pd.read_csv(raw.parent / "file_metadata.tsv", sep="\t")
    meta = meta[meta.sample_type == "Primary Tumor"]
    fpkm = pd.read_csv(paths.TCGA_BULK_TUMOR_ONLY, sep="\t", index_col=0)
    fpkm = fpkm[fpkm.index.str.startswith("ENSG")]
    patients = pd.read_csv(paths.TCGA_BAYESPRISM_DIR / "tcga_bulk_counts.tsv",
                           sep="\t", index_col=0, usecols=[0]).index
    files = {p.name: p for p in raw.rglob("*.augmented_star_gene_counts.tsv")}
    cols, sym = {}, None
    for s in patients:
        names = sorted(meta.loc[meta.submitter_id == s, "file_name"])
        if not names:
            continue
        tabs = [pd.read_csv(files[f], sep="\t", comment="#", index_col=0) for f in names]
        tabs = [t[t.index.str.startswith("ENSG")] for t in tabs]
        pick = tabs[0]
        if len(tabs) > 1:
            match = [t for t in tabs if s in fpkm.columns and np.allclose(
                t.fpkm_unstranded.reindex(fpkm.index).values, fpkm[s].values,
                atol=1e-3, equal_nan=True)]
            if len(match) != 1:
                raise SystemExit(f"{s}: {len(match)} tumour aliquots match the tumour-only matrix")
            pick = match[0]
        cols[s], sym = pick.unstranded, pick.gene_name
    bk = pd.DataFrame(cols)
    bk["sym"] = sym
    bk["m"] = bk[list(cols)].mean(axis=1)
    bk = bk.sort_values("m", ascending=False).drop_duplicates("sym").set_index("sym")
    return bk.drop(columns="m").T


def main():
    ct = cell_labels()
    with h5py.File(H5AD, "r") as f:
        obs, var = h5_index(f["obs"]), h5_index(f["var"])
        rows = obs.get_indexer(ct.index)
        if (rows < 0).any():
            raise SystemExit("reference cells missing from the deposited h5ad")
        counts = read_csr(f, "layers/counts", rows, len(var))
        if not np.allclose(counts, np.round(counts)):
            raise SystemExit("layers['counts'] is not integer")
        lognorm = read_csr(f, "X", rows, len(var))

    for d in ("tcga", "prjeb25780"):
        (WORK / d).mkdir(parents=True, exist_ok=True)
    ct.to_csv(WORK / "sc_cell_types.tsv", sep="\t")

    sc = pd.DataFrame(np.round(counts).astype(np.int32), index=ct.index, columns=var)
    sc = sc.loc[:, ~sc.columns.duplicated()]
    sc.to_csv(WORK / "tcga" / "sc_counts.tsv", sep="\t")
    tcga_tumour_counts().to_csv(WORK / "tcga" / "bulk_counts.tsv", sep="\t")

    lists = Path(__file__).resolve().parents[1] / "outputs" / "prjeb25780"
    genes = (lists / "reference_genes.txt").read_text().split()
    samples = (lists / "samples.txt").read_text().split()
    ref = pd.DataFrame(np.expm1(lognorm), index=ct.index, columns=var)
    ref = ref.loc[:, ~ref.columns.duplicated()]
    ref[genes].to_csv(WORK / "prjeb25780" / "sc_counts.tsv", sep="\t")
    pd.read_csv(paths.TCGA_BAYESPRISM_DIR / "sc_cell_types.tsv", sep="\t").to_csv(
        WORK / "prjeb25780" / "sc_cell_types.tsv", sep="\t", index=False)
    bulk = pd.read_csv(paths.TIGER_EXPR, sep="\t", index_col=0)
    bulk.loc[genes, samples].T.to_csv(WORK / "prjeb25780" / "bulk_counts.tsv", sep="\t")
    print("written:", WORK)


if __name__ == "__main__":
    main()
