#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""Compare two GraphST output h5ads: bytes first, then every stored value.

Byte equality is the smoke-test criterion.  The value-level pass exists so a
byte difference can be ATTRIBUTED -- HDF5 container noise vs numbers actually
moving -- rather than merely reported.

Large datasets (map_matrix is n_spot x n_cell, ~3.7 GB) are compared in row
chunks; loading two copies as float64 would need ~30 GB.

usage: compare_h5ad.py A.h5ad B.h5ad [--chunk-rows N]
"""
import sys, os, hashlib, filecmp
import numpy as np, h5py

A, B = sys.argv[1], sys.argv[2]
CHUNK = 64
if "--chunk-rows" in sys.argv:
    CHUNK = int(sys.argv[sys.argv.index("--chunk-rows") + 1])

def md5(p, chunk=1 << 24):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

print(f"A {A}\n  bytes={os.path.getsize(A)} md5={md5(A)}")
print(f"B {B}\n  bytes={os.path.getsize(B)} md5={md5(B)}")
print(f"\nBYTE_IDENTICAL {'yes' if filecmp.cmp(A, B, shallow=False) else 'no'}")

def collect(g, out, pre=""):
    for k, v in g.items():
        p = f"{pre}/{k}"
        collect(v, out, p) if isinstance(v, h5py.Group) else out.append(p)

def maxdiff(da, db):
    """max |a-b| and max |a|, streamed in row chunks."""
    m = 0.0; scale = 0.0
    if da.ndim == 0 or da.size == 0:
        va = np.asarray(da[()], float); vb = np.asarray(db[()], float)
        if va.size:
            m = float(np.nanmax(np.abs(va - vb))); scale = float(np.nanmax(np.abs(va)))
        return m, scale
    n = da.shape[0]
    step = max(1, CHUNK) if da.size > (1 << 24) else n
    for i in range(0, n, step):
        va = np.asarray(da[i:i + step], dtype=np.float64)
        vb = np.asarray(db[i:i + step], dtype=np.float64)
        d = np.abs(va - vb)
        if d.size:
            m = max(m, float(np.nanmax(d))); scale = max(scale, float(np.nanmax(np.abs(va))))
    return m, scale

pa, pb = [], []
with h5py.File(A, "r") as fa, h5py.File(B, "r") as fb:
    collect(fa, pa); collect(fb, pb)
    sa, sb = set(pa), set(pb)
    if sa - sb: print(f"DATASETS_ONLY_IN_A {len(sa-sb)}: {sorted(sa-sb)[:10]}")
    if sb - sa: print(f"DATASETS_ONLY_IN_B {len(sb-sa)}: {sorted(sb-sa)[:10]}")

    worst, n_diff, n_float = [], 0, 0
    for p in sorted(sa & sb):
        da, db = fa[p], fb[p]
        if da.shape != db.shape:
            print(f"SHAPE_DIFF {p}: {da.shape} vs {db.shape}"); n_diff += 1; continue
        if getattr(da, "dtype", None) is not None and da.dtype.kind in "fc":
            n_float += 1
            m, scale = maxdiff(da, db)
            worst.append((m, m / scale if scale else m, p))
            if m > 0: n_diff += 1
        else:
            if not np.array_equal(da[()], db[()]):
                print(f"VALUE_DIFF (non-float) {p}"); n_diff += 1

worst.sort(reverse=True)
print(f"\nFLOAT_DATASETS_COMPARED {n_float}")
print(f"DATASETS_WITH_ANY_DIFFERENCE {n_diff}")
print("\nlargest absolute differences:")
for m, rel, p in worst[:15]:
    print(f"  {m:.6e}  (rel {rel:.3e})  {p}")
if worst:
    print(f"\nMAX_ABS_DIFF_OVERALL {worst[0][0]:.6e}")
    print(f"MAX_REL_DIFF_OVERALL {max(w[1] for w in worst):.6e}")
