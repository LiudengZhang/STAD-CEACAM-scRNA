"""Decide which matrix of an AnnData holds the log1p CP10K expression.

The analysis inputs keep the expression in `.raw` (their `.X` is scaled, with
NaN rows). The deposited objects have no `.raw` and hold the expression in
`.X`, with integer counts in layers['counts']. expression_source() returns
`.raw` where it exists and passes the checks below, otherwise the object
itself if its `.X` is identified by digest, and raises in every other case.

Checks:
  - `.raw` must pass the integer-ladder test (log1p_cp10k_problems) and not
    be listed in DAMAGED_RAW (full_dataset.h5ad and Neutrophils.h5ad carry a
    `.raw` that was normalised twice).
  - Without `.raw`, `.X` must be sparse, finite and non-negative with a
    counts layer (promotion_problems), and the md5 of the file's /X/data must
    match PROVEN_X or SOUND_NOT_RAW for that file name. `path` is therefore
    required when `.raw` is absent.

Use:
    from shared.expression import expression_source, expression_adata

    src = expression_source(adata, path)    # has .X and .var_names
    i = list(src.var_names).index("CEACAM5")
    x = src.X[:, i]

    ad = expression_adata(adata, path)      # a full AnnData, obs carried
"""

from __future__ import annotations

__all__ = ["expression_source", "expression_adata", "promotion_problems",
           "x_data_md5", "file_x_data_md5", "assert_x_stands_in_for_raw",
           "log1p_cp10k_problems", "raw_problems", "assert_raw_is_sound",
           "PROVEN_X", "DIFFERENT_X", "DAMAGED_RAW", "SOUND_NOT_RAW",
           "LADDER_CELLS", "LADDER_TOL", "LADDER_SEED"]


def _describe(adata, name):
    return name or getattr(adata, "filename", None) or "this object"


def promotion_problems(adata):
    """Why `.X` may not be read as expression. Empty means it may."""
    import numpy as np
    from scipy import sparse

    problems = []

    layers = getattr(adata, "layers", None)
    if layers is None or "counts" not in layers:
        problems.append(
            "layers['counts'] is absent - build_clean_h5ad.py and "
            "attach_counts_layer.py put the counts there, so an object without "
            "it is not the promoted deposit")

    X = adata.X
    if not sparse.issparse(X):
        problems.append(
            "`.X` is dense - the promoted matrix is csr float32; the dense one "
            "is the scaled matrix the double normalisation damaged")
    else:
        data = X.data
        if data.size:
            if not bool(np.isfinite(data).all()):
                problems.append(
                    "`.X` contains NaN or inf - a log1p CP10K matrix cannot; "
                    "see 00_Data_Audit/FINDINGS.md sections 1 and 7")
            if float(data.min()) < 0.0:
                problems.append(
                    f"`.X` has negative values (min {float(data.min()):.4g}) - "
                    f"a log1p CP10K matrix cannot")
    return problems


# Identity tables: md5 of each deposited object's /X/data, keyed by file name.

# File name -> md5 of `.X.data` for deposited objects whose `.X` is
# byte-identical to the input's `.raw`.
PROVEN_X = {
    "B_cells.h5ad":     "c9e74502a892e5b1f2140232819f4392",
    "DC_cells.h5ad":    "43cbdd42f8aedf93a2b2546967cb912f",
    "Endothelial.h5ad": "f39499d56c7bdfc14aa1b2b74f49af0e",
    "Epithelial.h5ad":  "8f8d3360369ac3e36a4b3724419302b9",
    "Fibroblast.h5ad":  "3c3131b2d6fffd14947cbac016a9a387",
    "MoMac.h5ad":       "9ad0d7e28d4434e1abb72410a11b1f10",
    "NK_cells.h5ad":    "af714cea22320439e867b12902e4db74",
    "Pericyte.h5ad":    "aae832c2ba56bdc269ef2a1bed434a7c",
    "TCD4.h5ad":        "5ddf48a7e0c7cbbec0350a10d1f25f94",
    "TCD8.h5ad":        "b10a3df1c9f860d2848590662c5e1bb8",
}

# File name -> (md5 of the deposited `.X.data`, why it differs from the
# input's `.raw`). SOUND_NOT_RAW is checked first.
DIFFERENT_X = {
    "Neutrophils.h5ad": (
        "fb1b5c027960dd617d288dbc015befab",
        "the deposited object is Neutrophils_sound.h5ad, the rebuild that "
        "repaired the double normalisation; same shape, same nnz, same "
        "sparsity, different values. The published neutrophil numbers came "
        "from the author's `.raw`, which is the damaged matrix. Reading this "
        "one instead is what moved Supplementary Figure S7 E's P from 0.200 "
        "to 0.583"),
    "full_dataset.h5ad": (
        "226abe3ca32fb22ccfdf23b5bdc9541d",
        "667,573,017 non-zeros against the author's 667,390,633 at the same "
        "542,121 x 57,058 shape - a different sparsity pattern, so not a "
        "rescaling of the same matrix"),
}


# Is `.raw` log1p CP10K? Within a cell the implied counts expm1(v) / expm1(min v)
# must be 1, 2, 3, ... A second normalisation keeps the ranks and the 1e4 row
# sum but puts the implied counts about halfway between integers (gap ~0.5,
# against <4e-4 for a correct matrix).

# Cells sampled for the ladder, the tolerance, and the seed.
LADDER_CELLS = 40
LADDER_TOL = 0.01
LADDER_SEED = 42

# File name -> (measured ladder gap, explanation). These `.raw` matrices were
# normalised twice and may not be read as expression.
DAMAGED_RAW = {
    "full_dataset.h5ad": (
        0.499034,
        "its `.raw` holds the log-normalised values put through "
        "normalize_total and log1p a second time. The sound matrix is the "
        "deposited `06_Clean_Data/01_H5AD/full_dataset.h5ad` `.X` "
        "(paths.FULL_DATASET_DEPOSIT_H5AD), which reproduces "
        "log1p(counts/rowsum*1e4) from that file's own `layers['counts']` to "
        "1e-6. Figure 5H already reads the deposit and says so"),
    "Neutrophils.h5ad": (
        0.499817,
        "its `.raw` is log1p(CP10K(x)) applied to an already log1p CP10K "
        "matrix - the defect in "
        "upstream-pipeline/.../Neutrophils_filtered_v2.1/filter_and_visualize.py. The "
        "sound matrix is `06_Clean_Data/01_H5AD/Neutrophils_sound.h5ad` "
        "(paths.NEUTROPHILS_SOUND_H5AD), which keeps all 61,167 cells and the "
        "same 4 R / 5 NR split - a change of SOURCE, not of COHORT, proven in "
        "06_Clean_Data/verify_neutrophil_rebuild.py"),
}

# File name -> (md5 of the deposited `.X.data`, explanation). These deposited
# `.X` matrices replace a `.raw` that was normalised twice; reading them is
# allowed and reported on stderr.
SOUND_NOT_RAW = {
    "Neutrophils.h5ad": (
        "fb1b5c027960dd617d288dbc015befab",
        "the sound rebuild, byte-identical to Neutrophils_sound.h5ad on all "
        "three CSR arrays. Reading it rather than the author's damaged `.raw` "
        "is what moved Supplementary Figure S7 E's P from 0.200 to 0.583"),
    "full_dataset.h5ad": (
        "226abe3ca32fb22ccfdf23b5bdc9541d",
        "667,573,017 non-zeros against the author's 667,390,633 at the same "
        "542,121 x 57,058 shape. The extra non-zeros are the ~22,000 genes "
        "the symbol-matched reconstruction dropped and the Ensembl-matched "
        "one keeps; this `.X` reproduces its own integer counts to 1e-6 and "
        "the author's `.raw` does not reproduce them at all"),
}


def log1p_cp10k_problems(X, label="this matrix"):
    """Why `X` is not log1p CP10K. Empty means it is."""
    import numpy as np
    from scipy import sparse

    Xc = X.tocsr() if sparse.issparse(X) else np.asarray(X)
    n = Xc.shape[0]
    if n == 0:
        return [f"{label} has no cells, so the ladder cannot be run"]
    rng = np.random.default_rng(LADDER_SEED)
    rows = rng.choice(n, min(LADDER_CELLS, n), replace=False)
    worst, tested, integral = 0.0, 0, 0
    for i in rows:
        i = int(i)
        if sparse.issparse(Xc):
            v = Xc.data[Xc.indptr[i]:Xc.indptr[i + 1]]
        else:
            row = Xc[i]
            v = row[row > 0]
        v = np.asarray(v, dtype=np.float64)
        if v.size < 20:
            continue
        tested += 1
        if np.abs(v - np.round(v)).max() < 1e-6:
            integral += 1
            continue
        c = np.expm1(v) / np.expm1(v.min())
        worst = max(worst, float(np.abs(c - np.round(c)).max()))
    if tested == 0:
        return [f"{label}: too few detected genes in {LADDER_CELLS} sampled "
                f"cells to run the ladder, so nothing can be said about it"]
    if integral > tested / 2:
        return [f"{label} holds INTEGER COUNTS, not log1p CP10K "
                f"({integral} of {tested} sampled cells are integral)"]
    if worst >= LADDER_TOL:
        return [f"{label} is NOT log1p CP10K: the counts it implies miss the "
                f"integers by up to {worst:.4f}, against a tolerance of "
                f"{LADDER_TOL} and {3.2e-4:.1e} for the worst sound matrix in "
                f"this project. A matrix normalised twice keeps its ranks and "
                f"still sums to 1e4, so only this spacing gives it away"]
    return []


def raw_problems(adata, path=None):
    """Why this object's `.raw` may not be read as expression. Empty means it may."""
    from pathlib import Path as _Path
    raw = getattr(adata, "raw", None)
    if raw is None:
        return []
    problems = []
    if path is not None:
        name = _Path(path).name
        if name in DAMAGED_RAW:
            gap, why = DAMAGED_RAW[name]
            problems.append(
                f"{name}'s `.raw` was MEASURED on 2026-10-05 and is not "
                f"log1p CP10K - implied counts miss the integers by "
                f"{gap:.4f}: {why}")
    problems.extend(log1p_cp10k_problems(
        raw.X, f"{_Path(path).name if path else 'this object'}'s `.raw`"))
    return problems


def assert_raw_is_sound(adata, path=None):
    """Raise unless this object's `.raw` is a log1p CP10K matrix."""
    problems = raw_problems(adata, path)
    if problems:
        from pathlib import Path as _Path
        label = _describe(adata, _Path(path).name if path else None)
        raise SystemExit(
            f"{label}'s `.raw` may not be read as the log-normalised "
            f"expression matrix:\n  - " + "\n  - ".join(problems)
            + f"\nRefusing to read it. Point the "
              f"script at the correctly normalised object named above.")


def x_data_md5(adata):
    """md5 of `.X.data`'s bytes, for an object held WHOLE in memory."""
    import hashlib
    import numpy as np
    return hashlib.md5(
        np.ascontiguousarray(adata.X.data).tobytes()).hexdigest()


def file_x_data_md5(path):
    """md5 of `/X/data` read off the FILE, for callers holding a subset."""
    import hashlib
    import h5py
    import numpy as np
    h = hashlib.md5()
    with h5py.File(str(path), "r") as f:
        ds = f["X/data"]
        step = 1 << 24
        for i in range(0, ds.shape[0], step):
            h.update(np.ascontiguousarray(ds[i:i + step]).tobytes())
    return h.hexdigest()


def _verdict(name, digest):
    """Check a file name and its `.X.data` digest against the identity tables.
    Returns None if allowed, else the reason."""
    if not name:
        return ("the caller did not say which object this is. The identity "
                "table is keyed by file name, because two deposited objects "
                "hold the same bytes under different names and those bytes "
                "are right under one name and wrong under the other. Pass "
                "`name=<path>.name`")
    if name in SOUND_NOT_RAW:
        want, why = SOUND_NOT_RAW[name]
        if digest == want:
            # Allowed, but reported on stderr.
            import sys
            print(f"NOTE  {name}: reading the deposited `.X`, which is the "
                  f"correctly normalised matrix - {why}. The original object's `.raw` is not "
                  f"log1p CP10K (gap "
                  f"{DAMAGED_RAW.get(name, (float('nan'), ''))[0]:.4f}).",
                  file=sys.stderr)
            return None
        return (f"{name} was MEASURED on 2026-10-05 and this is not what the "
                f"measurement found: its `.X.data` digests {digest}, where the "
                f"sound deposited copy digests {want}. Neither the sound "
                f"deposit nor the author's `.raw`, under a name that has a "
                f"measurement - which is the one thing the table exists to "
                f"catch")
    if name in DIFFERENT_X:
        want, why = DIFFERENT_X[name]
        where = ("and its `.X` is the matrix that measurement found there"
                 if digest == want else
                 f"and its `.X` digest {digest} matches neither the author's "
                 f"`.raw` nor the deposited copy measured on 2026-10-05")
        return (f"{name} was MEASURED on 2026-10-05 and its deposited `.X` is "
                f"NOT the author's `.raw` {where}: {why}")
    if name not in PROVEN_X:
        return (f"{name} is not in PROVEN_X. No measurement says its `.X` is "
                f"the matrix the published numbers came from, and an "
                f"unmeasured object is refused rather than assumed. Measure "
                f"it against the author's `.raw` and record the digest")
    want = PROVEN_X[name]
    if digest != want:
        return (f"{name}'s `.X.data` digest is {digest}; the object proven "
                f"identical to the author's `.raw` on 2026-10-05 digests "
                f"{want}. This is a different matrix under a name that was "
                f"cleared, which is the one thing the table exists to catch")
    return None


def assert_x_stands_in_for_raw(adata, path):
    """Raise unless this object's `.X` may be read where `.raw` was asked for."""
    from pathlib import Path as _P
    if getattr(adata, "raw", None) is not None:
        assert_raw_is_sound(adata, path)
        return
    name = _P(path).name
    problems = promotion_problems(adata)
    if problems:
        raise SystemExit(
            f"{name} has no `.raw`, and its `.X` cannot be read as the "
            f"log-normalised expression matrix:\n  - " + "\n  - ".join(problems)
            + "\nRefusing to read it.")
    problem = _identity_problem(path)
    if problem:
        raise SystemExit(
            f"{name} has no `.raw`, and its `.X` may not be read in its "
            f"place:\n  - {problem}\nRefusing to read it. Where the two "
            f"matrices are not PROVEN to be the same, the script stops rather "
            f"than substitute one for the other.")


def _identity_problem(path):
    """Why this `.X` may not stand in for `.raw`. None means it may."""
    if path is None:
        return ("the caller did not say which object this is. The identity "
                "table is keyed by file name, because two deposited objects "
                "hold the same bytes under different names and those bytes "
                "are right under one name and wrong under the other. Pass "
                "`path=<the h5ad>`")
    from pathlib import Path as _Path
    p = _Path(path)
    if not p.is_file():
        return (f"`{path}` is not a file. The proof is a digest of the h5ad, "
                f"so the argument has to be the path it was read from and not "
                f"its name - a name cannot be resolved to bytes")
    return _verdict(p.name, file_x_data_md5(p))


def _check_promoted(adata, path):
    from pathlib import Path as _Path
    label = _describe(adata, _Path(path).name if path else None)
    problems = promotion_problems(adata)
    if problems:
        raise SystemExit(
            f"{label} has no `.raw`, and its `.X` cannot be "
            f"read as the log-normalised expression matrix:\n  - "
            + "\n  - ".join(problems)
            + "\nRefusing to read it. The analysis inputs keep the expression "
              "in `.raw`; the deposited clean h5ads promote it into `.X`. "
              "Nothing else is either.")
    # Structure says it is a deposit; the digest says which matrix it holds.
    problem = _identity_problem(path)
    if problem:
        raise SystemExit(
            f"{label} has no `.raw`, and its `.X` may not be "
            f"read in its place:\n  - {problem}\nRefusing to read it. Where "
            f"the two matrices are not PROVEN to be the same, the script stops "
            f"rather than substitute one for the other.")


def expression_source(adata, path=None):
    """The object whose `.X` is the log-normalised expression matrix.

    `.raw` where it passes raw_problems(); otherwise the object itself if
    it is a deposit whose `.X` is proven by digest. Raises in every other case.

    `path` is the h5ad `adata` was read from (it may since have been subset).
    Used for the `.raw` name check, and required when there is no `.raw`.
    """
    if getattr(adata, "raw", None) is not None:
        assert_raw_is_sound(adata, path)
        return adata.raw
    _check_promoted(adata, path)
    return adata


def expression_adata(adata, path=None):
    """The same, as a full AnnData with `obs` carried.

    With a `.raw`, returns `.raw.to_adata()` indexed by the parent's
    obs names.
    """
    if getattr(adata, "raw", None) is not None:
        assert_raw_is_sound(adata, path)           # see expression_source
        return adata.raw.to_adata()[adata.obs_names]
    _check_promoted(adata, path)
    return adata


if __name__ == "__main__":
    from expression_selftest import self_test
    raise SystemExit(self_test())
