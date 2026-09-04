#!/usr/bin/env python3
"""Preflight for the GraphST reproduction.

1. print the installed versions and diff them against the versions measured in
   the original environment (Liudeng_Python_310);
2. exercise the exact GraphST call chain the runner uses -- GraphSTModel(...,
   deconvolution=True).train_map() and project_cell_to_spot -- on small
   synthetic arrays, with the determinism pinning active, so that a refusal
   from torch.use_deterministic_algorithms(True) shows up in minutes instead of
   two hours into a real sample.

No project data is read and nothing is written.
"""
import importlib.metadata as md
import sys, traceback

PINNED = {
    "torch": "2.8.0", "GraphST": "1.1.1", "scanpy": "1.9.6", "anndata": "0.10.9",
    "numpy": "1.26.4", "scipy": "1.15.3", "scikit-learn": "1.7.2",
    "pandas": "2.3.3", "POT": "0.9.6.post1", "leidenalg": "0.10.2",
    "igraph": "0.11.8", "matplotlib": "3.10.5", "PyYAML": "6.0.3",
    # transitive deps that the first smoke attempt exposed: scikit-misc was
    # missing entirely (scanpy seurat_v3 HVG, used by GraphST.preprocess) and
    # numba/llvmlite/pynndescent/umap-learn had drifted -- they drive
    # sc.pp.neighbors -> leiden -> spatial_domain.
    "scikit-misc": "0.1.4", "numba": "0.65.1", "llvmlite": "0.47.0",
    "pynndescent": "0.5.13", "umap-learn": "0.5.7", "joblib": "1.5.3",
    "h5py": "3.16.0",
}

print("=" * 72); print("VERSION DIFF vs the original environment"); print("=" * 72)
mismatch = 0
for pkg, want in PINNED.items():
    try:
        got = md.version(pkg)
    except Exception as exc:
        got = f"<{type(exc).__name__}>"
    ok = (got == want)
    mismatch += (not ok)
    print(f"  {pkg:<16} want {want:<14} got {got:<14} {'OK' if ok else 'MISMATCH'}")

import torch
print(f"\n  torch.__version__ (full build tag): {torch.__version__}")
print(f"  torch.version.cuda: {torch.version.cuda}")
print(f"  torch.cuda.is_available(): {torch.cuda.is_available()}")
print(f"  torch.get_num_threads(): {torch.get_num_threads()}")
print(f"  torch.are_deterministic_algorithms_enabled(): "
      f"{torch.are_deterministic_algorithms_enabled()}")
print(f"\nVERSION_MISMATCHES {mismatch}")

print("\n" + "=" * 72); print("GraphST call-chain probe (synthetic, deconvolution=True)"); print("=" * 72)
import numpy as np, anndata as ad, scanpy as sc
from GraphST.GraphST import GraphST as GraphSTModel
from GraphST.utils import project_cell_to_spot

rng = np.random.default_rng(0)
n_spot, n_cell, n_gene = 60, 200, 120
genes = [f"g{i}" for i in range(n_gene)]

sp = ad.AnnData(rng.poisson(2.0, (n_spot, n_gene)).astype(np.float32),
                var=__import__("pandas").DataFrame(index=genes))
sp.obsm["spatial"] = rng.integers(0, 500, (n_spot, 2)).astype(float)

pd = __import__("pandas")
scd = ad.AnnData(rng.poisson(2.0, (n_cell, n_gene)).astype(np.float32),
                 var=pd.DataFrame(index=genes))
scd.obs["cell_type"] = pd.Categorical(rng.choice(["A", "B", "C"], n_cell))

for a in (sp, scd):
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sp.var["highly_variable"] = True
scd.var["highly_variable"] = True

rc = 0
try:
    m = GraphSTModel(sp, scd, random_seed=42, epochs=5, learning_rate=0.001,
                     deconvolution=True, device="cpu")
    out_sp, out_sc = m.train_map()
    project_cell_to_spot(out_sp, out_sc, retain_percent=0.15)
    got = [c for c in ("A", "B", "C") if c in out_sp.obs.columns]
    print(f"\n  train_map + project_cell_to_spot OK; cell-type columns written: {got}")
    print("  PROBE_OK")
except Exception:
    rc = 1
    print("\n  PROBE_FAILED -- traceback verbatim:")
    traceback.print_exc()

sys.exit(rc)
