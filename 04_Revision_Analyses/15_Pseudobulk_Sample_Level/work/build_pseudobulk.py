"""
Sample-level (pseudobulk) matrices for the NF-kB recompute.

12_R1.8_DEG_Recompute/scripts/recompute_deg.py tests cells as if they were
independent replicates: `rank_genes_groups(..., method='t-test')` over tens of
thousands of cells drawn from about a dozen patients. This module changes the
unit of replication to the sample and nothing else. The cell-type sources, the
R/NR labels, the pre/post split and the cell selection are recompute_deg.py's,
copied rather than rewritten, so that the statistics are the only difference.

What is summed
--------------
Integer UMI counts, from `layers['counts']`, summed over the cells of one
sample. Summing log-normalised values is not a pseudobulk and is not done here;
neither is back-calculating counts from log values, which would be fabricated
data.

  eleven cell types   06_Clean_Data/01_H5AD/<type>.h5ad, layers['counts']
                      attached 31 Aug 2026 by attach_counts_layer.py and
                      verified per file in counts_layer_report.csv: integer,
                      100% of cells matched to their per-sample droplet
                      matrix, zero genes zero-filled, .X byte-identical after.
  Neutrophils         06_Clean_Data/02_Rebuilt/Neutrophils_sound.h5ad, whose
                      layers['counts'] is the same recovery. The shipped
                      Neutrophils.h5ad is doubly normalised in both .X and
                      .raw; the rebuilt object is the one 13_R1.8 used.
  Mast, Plasma        no counts layer exists: their objects are the upstream-pipeline
                      *_integrated.h5ad files, which carry .raw (log1p CP10K,
                      sound) and no counts anywhere. The counts are recovered
                      here the same way attach_counts_layer.py recovered the
                      other fourteen - from the 70 per-sample droplet
                      matrices, matched on Ensembl gene ID and on the
                      24-character barcode - by importing that script and
                      calling its own functions. Nothing is written back into
                      the upstream-pipeline files.

Sample identifiers
------------------
Never written out. Each sample is given a code SMP01.. from the sorted order of
the samples seen across the whole run; the mapping is not saved.

Outputs (all under ../outputs/pseudobulk/)
    <cell>_<phase>_counts.tsv.gz    genes x kept samples, integer
    <cell>_<phase>_samples.tsv      code, group, n_cells for the kept samples
    sample_census.csv               every sample of every contrast, kept or not
    build_report.csv                one row per contrast

Run: python build_pseudobulk.py [--cell CELL] [--phase pre|post]
"""

from pathlib import Path
import argparse
import gc
import os
import sys
import warnings
from contextlib import redirect_stdout

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                       # Round_7_major_revision
OUT = HERE.parents[0] / "outputs" / "pseudobulk"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT / "00_Config"))
from paths import PROJECT_ROOT                                       # noqa: E402

CLEAN = ROOT / "06_Clean_Data" / "01_H5AD"
REBUILT = ROOT / "06_Clean_Data" / "02_Rebuilt" / "Neutrophils_sound.h5ad"
ROUND_4_MAJOR = (PROJECT_ROOT.parent / "upstream-pipeline" / "04_Final_Panels"
                 / "00_Set_Ups" / "00_Data" / "01_Major_Cell_Types")

# recompute_deg.py's thirteen populations, on the prepared inputs. Every file
# holds exactly one population, so there is no obs label to filter on.
CELL_SOURCES = {
    "B_cells": CLEAN / "B_cells.h5ad",
    "DC_cells": CLEAN / "DC_cells.h5ad",
    "Endothelial_cells": CLEAN / "Endothelial.h5ad",
    "Epithelial": CLEAN / "Epithelial.h5ad",
    "Fibroblast": CLEAN / "Fibroblast.h5ad",
    "Mast_cells": ROUND_4_MAJOR / "mast_cell_integrated.h5ad",
    "MoMac": CLEAN / "MoMac.h5ad",
    "Neutrophils": REBUILT,
    "NK_cells": CLEAN / "NK_cells.h5ad",
    "Pericyte": CLEAN / "Pericyte.h5ad",
    "Plasma_cells": ROUND_4_MAJOR / "plasma_cell_integrated.h5ad",
    "TCD4_cells": CLEAN / "TCD4.h5ad",
    "TCD8_cells": CLEAN / "TCD8.h5ad",
}
NO_COUNTS_LAYER = {"Mast_cells", "Plasma_cells"}

PHASES = {
    "pre": dict(column="stomach_pre_grouping", ref="Responsed", test="No-response"),
    "post": dict(column="stomach_post_grouping", ref="Responsed", test="No-response"),
}

# A sample contributing fewer than this many cells to a contrast is dropped:
# its summed profile is too shallow to be a library. 10 is the threshold used
# throughout the pseudobulk literature and is stated in FINDINGS.md.
MIN_CELLS_PER_SAMPLE = 10
# Below this a contrast cannot be tested at all.
MIN_SAMPLES_PER_GROUP = 3

ROW_CHUNK = 20000


def read_str(f, path):
    d = f[path]
    if isinstance(d, h5py.Group) and "categories" in d:
        cats = np.asarray(d["categories"].asstr()[:])
        return cats[d["codes"][:]]
    return np.asarray(d.asstr()[:])


def obs_frame(path, columns):
    with h5py.File(path, "r") as f:
        idx = f["obs"].attrs.get("_index", "_index")
        out = {"_name": read_str(f, f"obs/{idx}")}
        for c in columns:
            if c not in f["obs"]:
                raise SystemExit(f"{path.name}: no obs['{c}']")
            out[c] = read_str(f, f"obs/{c}")
    return pd.DataFrame(out)


def var_names(path, group="var"):
    with h5py.File(path, "r") as f:
        idx = f[group].attrs.get("_index", "_index")
        return read_str(f, f"{group}/{idx}")


def sum_counts_layer(path, rows, sample_codes, n_codes):
    """Sum layers['counts'] over the given rows, grouped by sample.

    Read in row blocks straight off the HDF5 csr rather than through anndata:
    Epithelial is 149,373 x 56,034 with 204 million non-zeros and there is no
    reason to hold it.
    """
    with h5py.File(path, "r") as f:
        g = f["layers/counts"]
        if g.attrs.get("encoding-type") != "csr_matrix":
            raise SystemExit(f"{path.name}: layers['counts'] is not csr")
        n_genes = int(g.attrs["shape"][1])
        indptr = g["indptr"][:]
        total = np.zeros((n_codes, n_genes), dtype=np.float64)
        rows = np.asarray(rows)
        for start in range(0, len(rows), ROW_CHUNK):
            blk = rows[start:start + ROW_CHUNK]
            lo, hi = int(indptr[blk.min()]), int(indptr[blk.max() + 1])
            data = g["data"][lo:hi]
            indices = g["indices"][lo:hi]
            ptr = indptr[blk.min():blk.max() + 2] - lo
            M = sparse.csr_matrix((data, indices, ptr),
                                  shape=(blk.max() - blk.min() + 1, n_genes))
            M = M[blk - blk.min()]
            ind = sparse.csr_matrix(
                (np.ones(len(blk)), (sample_codes[start:start + len(blk)],
                                     np.arange(len(blk)))),
                shape=(n_codes, len(blk)))
            total += np.asarray((ind @ M).todense())
            del M, data, indices
            gc.collect()
    return total


def sum_from_droplet_cache(path, rows, sample_codes, n_codes):
    """Counts for Mast and Plasma, recovered from the 70 per-sample matrices.

    attach_counts_layer.py is imported and its own functions are called, so the
    barcode rule, the Ensembl-ID gene bridge and the "a barcode that cannot be
    matched raises" behaviour are that script's and not a second copy of them.
    """
    sys.path.insert(0, str(ROOT / "06_Clean_Data"))
    import attach_counts_layer as acl

    genes = sorted(set().union(*(set(acl.read_axes(p)[3])
                                 for p in acl.clean_files())))
    gene_pos = {g: i for i, g in enumerate(genes)}

    with h5py.File(path, "r") as f:
        idx = f["obs"].attrs.get("_index", "_index")
        names = read_str(f, f"obs/{idx}")
        sample = read_str(f, "obs/sample")
        gid = read_str(f, "raw/var/gene_id")
        vidx = f["raw/var"].attrs.get("_index", "_index")
        symbols = read_str(f, f"raw/var/{vidx}")
    barcode = np.array([n[:acl.BARCODE_LEN] for n in names])
    tail = np.array([n[acl.BARCODE_LEN:] for n in names])
    bad = [(t, s) for t, s in zip(tail, sample)
           if not (t == "-" + s or t.startswith("-" + s + "-"))]
    if bad:
        raise SystemExit(f"{path.name}: {len(bad)} cell names do not read as "
                         f"<{acl.BARCODE_LEN}-char barcode>-<sample>")

    colmap = np.array([gene_pos[g] for g in gid])   # KeyError if a gene is new
    n_genes = len(gid)
    total = np.zeros((n_codes, n_genes), dtype=np.float64)
    rows = np.asarray(rows)
    code_of = dict(zip(rows, sample_codes))
    inv = np.full(len(genes), -1, dtype=np.int64)
    inv[colmap] = np.arange(n_genes)

    # The check that binds, and the one attach_counts_layer.py used: a
    # log-normalisation is strictly increasing and maps zero to zero, so a
    # correctly matched cell has exactly the same non-zero genes in the
    # recovered counts as in .raw. Anything less than 100% means the barcode
    # bridge put another cell's counts on this row.
    spot_rows, spot_sets = [], []
    rng = np.random.default_rng(42)
    want_spot = set(rng.choice(rows, min(200, len(rows)), replace=False).tolist())

    cache = {}
    for s in sorted(set(sample[rows])):
        idx_s = rows[sample[rows] == s]
        # acl.load_block prints the sample name; nothing here may write an
        # internal specimen identifier to a log.
        with open(os.devnull, "w") as devnull, redirect_stdout(devnull):
            block, bpos = acl.load_block(cache, s)
        miss = [b for b in barcode[idx_s] if b not in bpos]
        if miss:
            raise SystemExit(
                f"{path.name}: {len(miss)} barcodes of one sample are not in "
                f"its droplet matrix; counts are not invented.")
        local = np.array([bpos[b] for b in barcode[idx_s]])
        sub = acl.reindex_columns(block[local], inv, n_genes)
        codes = np.array([code_of[i] for i in idx_s])
        ind = sparse.csr_matrix(
            (np.ones(len(idx_s)), (codes, np.arange(len(idx_s)))),
            shape=(n_codes, len(idx_s)))
        total += np.asarray((ind @ sub).todense())
        for j, r in enumerate(idx_s):
            if r in want_spot:
                spot_rows.append(int(r))
                spot_sets.append(set(sub.indices[sub.indptr[j]:sub.indptr[j + 1]]))
        cache.clear()
        gc.collect()

    exact = 0
    with h5py.File(path, "r") as f:
        rp = f["raw/X"]
        rptr = rp["indptr"][:]
        for r, cs in zip(spot_rows, spot_sets):
            lo, hi = int(rptr[r]), int(rptr[r + 1])
            if set(rp["indices"][lo:hi]) == cs:
                exact += 1
    frac = exact / max(len(spot_rows), 1)
    print(f"    droplet recovery support check: {exact}/{len(spot_rows)} cells "
          f"match .raw exactly ({frac:.1%})", flush=True)
    if frac < 1.0:
        raise SystemExit(
            f"{path.name}: recovered counts do not reproduce .raw's non-zero "
            f"support on {len(spot_rows) - exact} of {len(spot_rows)} cells. "
            f"The barcode bridge is wrong; refusing to sum them.")
    return total, symbols


def one(cell, phase, code_of_sample):
    cfg = PHASES[phase]
    path = CELL_SOURCES[cell]
    obs = obs_frame(path, ["sample", cfg["column"]])
    keep = obs[cfg["column"]].isin([cfg["ref"], cfg["test"]]).values
    rows = np.flatnonzero(keep)
    if len(rows) == 0:
        return dict(cell=cell, phase=phase, status="no cells"), None, None

    samples = obs["sample"].values[rows]
    groups = obs[cfg["column"]].values[rows]
    uniq = sorted(set(samples))
    local = {s: i for i, s in enumerate(uniq)}
    codes = np.array([local[s] for s in samples])

    # A sample must sit in exactly one response group or the label is not a
    # sample-level variable and a sample-level test is meaningless.
    g_of = {}
    for s, g in zip(samples, groups):
        g_of.setdefault(s, set()).add(g)
    mixed = [s for s, v in g_of.items() if len(v) > 1]
    if mixed:
        raise SystemExit(f"{cell} {phase}: {len(mixed)} samples carry both "
                         f"response labels; the contrast is not sample-level")

    if cell in NO_COUNTS_LAYER:
        total, genes = sum_from_droplet_cache(path, rows, codes, len(uniq))
    else:
        total = sum_counts_layer(path, rows, codes, len(uniq))
        genes = var_names(path)

    n_cells = np.bincount(codes, minlength=len(uniq))
    census = pd.DataFrame(dict(
        cell=cell, phase=phase,
        sample_code=[code_of_sample(s) for s in uniq],
        group=[("NR" if next(iter(g_of[s])) == cfg["test"] else "R") for s in uniq],
        n_cells=n_cells,
        library=total.sum(axis=1).astype(np.int64),
        kept=n_cells >= MIN_CELLS_PER_SAMPLE))

    kept = census["kept"].values
    mat = pd.DataFrame(np.rint(total[kept]).astype(np.int64).T,
                       index=genes, columns=census.loc[kept, "sample_code"])
    mat = mat[mat.sum(axis=1) > 0]
    meta = census.loc[kept, ["sample_code", "group", "n_cells", "library"]]

    n_r = int((meta["group"] == "R").sum())
    n_nr = int((meta["group"] == "NR").sum())
    testable = n_r >= MIN_SAMPLES_PER_GROUP and n_nr >= MIN_SAMPLES_PER_GROUP
    if testable:
        mat.to_csv(OUT / f"{cell}_{phase}_counts.tsv.gz", sep="\t")
        meta.to_csv(OUT / f"{cell}_{phase}_samples.tsv", sep="\t", index=False)

    row = dict(
        cell=cell, phase=phase,
        status="ok" if testable else
        f"untestable {n_r}R/{n_nr}NR after the {MIN_CELLS_PER_SAMPLE}-cell rule",
        source=path.name, counts_source=("droplet matrices"
                                         if cell in NO_COUNTS_LAYER
                                         else "layers['counts']"),
        n_cells=int(len(rows)), n_samples_all=len(uniq),
        n_samples_kept=int(kept.sum()),
        n_samples_dropped=int((~kept).sum()),
        n_R=n_r, n_NR=n_nr,
        median_cells_per_sample=float(np.median(n_cells[kept])) if kept.any() else 0,
        min_cells_per_sample=int(n_cells[kept].min()) if kept.any() else 0,
        max_cells_per_sample=int(n_cells[kept].max()) if kept.any() else 0,
        n_genes_nonzero=int(len(mat)), testable=testable)
    return row, census, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", default=None)
    ap.add_argument("--phase", default=None, choices=list(PHASES))
    args = ap.parse_args()

    cells = [args.cell] if args.cell else list(CELL_SOURCES)
    phases = [args.phase] if args.phase else list(PHASES)
    missing = [c for c in cells if not CELL_SOURCES[c].exists()]
    if missing:
        raise SystemExit("no input file for: " + ", ".join(missing))

    # One code table for the whole run, from the sorted union of samples. The
    # mapping is never written.
    allsamp = set()
    for c in cells:
        allsamp |= set(obs_frame(CELL_SOURCES[c], ["sample"])["sample"])
    code = {s: f"SMP{i:02d}" for i, s in enumerate(sorted(allsamp), 1)}

    rows, cens = [], []
    for c in cells:
        for p in phases:
            print(f"[{c} {p}] start", flush=True)
            r, cn, _ = one(c, p, lambda s: code[s])
            rows.append(r)
            if cn is not None:
                cens.append(cn)
            print(f"[{c} {p}] {r['status']}  "
                  f"{r.get('n_samples_kept','?')} samples "
                  f"({r.get('n_R','?')}R/{r.get('n_NR','?')}NR) "
                  f"{r.get('n_cells','?')} cells", flush=True)

    tag = f"_{args.cell}_{args.phase}" if args.cell else ""
    pd.DataFrame(rows).to_csv(OUT.parent / f"build_report{tag}.csv", index=False)
    if cens:
        pd.concat(cens).to_csv(OUT.parent / f"sample_census{tag}.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
