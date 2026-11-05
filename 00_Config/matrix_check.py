"""Checks that an expression matrix really is log1p CP10K.

A matrix normalised twice keeps its ranks and still sums to 1e4 per cell, so
finiteness, non-negativity and rank-based statistics do not reveal it. Two
tests do, and each catches what the other misses:

  assert_log1p_cp10k    integer ladder: within a cell,
                        expm1(value) / expm1(min value) must be 1, 2, 3, ...
                        Needs only the matrix; blind to a wrong overall scale.

  assert_counts_ladder  count-stratum ladder: mean value must rise steeply with
                        the raw count behind each entry. Needs a `counts`
                        layer.

Run both where a counts layer exists, the first alone otherwise.

    python matrix_check.py --self-test
"""
import numpy as np

try:                                        # scipy is in every environment that
    from scipy import sparse                # reads an h5ad; guarded so that
except ImportError:                         # --self-test still runs without it
    sparse = None

# The ladders sample cells; the seed fixes which.
SEED = 42

# The integer ladder.
LADDER_CELLS = 40
LADDER_TOL = 0.01       # log1p CP10K matrices sit near 3e-4 (float32 precision)

# The count-stratum ladder: about 3.9 for log1p CP10K, about 1.1 for a doubly
# normalised matrix.
COUNT_STRATA = ((1, 1), (2, 2), (3, 5), (6, 20), (21, 100), (101, 10 ** 9))
STRATUM_CELLS = 200
MIN_STRATUM_ENTRIES = 50
COUNT_LADDER_MIN_RATIO = 2.0


class BadMatrix(RuntimeError):
    """The matrix is not what its container claims."""


def _rows(X, n):
    Xc = X.tocsr() if (sparse is not None and sparse.issparse(X)) else np.asarray(X)
    rng = np.random.default_rng(SEED)
    idx = rng.choice(Xc.shape[0], min(n, Xc.shape[0]), replace=False)
    return Xc, idx


def _row_values(Xc, i):
    if sparse is not None and sparse.issparse(Xc):
        return Xc.data[Xc.indptr[i]:Xc.indptr[i + 1]].astype(np.float64)
    row = Xc[i]
    return row[row > 0].astype(np.float64)


def assert_log1p_cp10k(X, source):
    """Refuse a matrix that is not log1p CP10K, by the integer-ladder test.

    log1p CP10K means

        value = log1p(count * 1e4 / library)

    so within one cell `expm1(value) / expm1(smallest value)` must be
    1, 2, 3, 4 ... The scale is taken from the matrix, not from
    obs["total_counts"], which may predate gene filtering.
    """
    Xc, rows = _rows(X, LADDER_CELLS)
    worst, tested, integral = 0.0, 0, 0
    for i in rows:
        v = _row_values(Xc, i)
        if v.size < 20:
            continue
        tested += 1
        if np.abs(v - np.round(v)).max() < 1e-6:
            integral += 1
            continue
        c = np.expm1(v) / np.expm1(v.min())      # one count is the smallest step
        worst = max(worst, float(np.abs(c - np.round(c)).max()))
    if tested == 0:
        raise BadMatrix(f"{source}: too few detected genes to test the matrix")
    if integral > tested / 2:
        raise BadMatrix(
            f"{source}: .raw holds integer counts, not log1p CP10K. "
            "Normalising it here would hide which file was wrong; fix the input.")
    if worst >= LADDER_TOL:
        raise BadMatrix(
            f"{source}: .raw is NOT log1p CP10K - implied counts miss the "
            f"integers by up to {worst:.3f} (correct matrices sit at 3e-4). This "
            "is the signature of a second normalisation on already normalised "
            "values. "
            "Refusing to test differential expression on it.")
    return worst


def count_ladder(X, counts, cells=STRATUM_CELLS):
    """Mean value of X per count stratum. -> [(label, n_entries, mean), ...]."""
    if X.shape != counts.shape:
        raise BadMatrix(f"the matrix is {X.shape} and the counts are "
                        f"{counts.shape}; they must be the same cells and genes")
    Xc, rows = _rows(X, cells)
    Cc = counts.tocsr() if (sparse is not None and sparse.issparse(counts)) \
        else np.asarray(counts)
    out = []
    for lo, hi in COUNT_STRATA:
        tot, n = 0.0, 0
        for i in rows:
            if sparse is not None and sparse.issparse(Cc):
                sl = slice(Cc.indptr[i], Cc.indptr[i + 1])
                cj, cv = Cc.indices[sl], Cc.data[sl]
            else:
                cj = np.flatnonzero(Cc[i])
                cv = Cc[i][cj]
            take = cj[(cv >= lo) & (cv <= hi)]
            if take.size == 0:
                continue
            if sparse is not None and sparse.issparse(Xc):
                xs = slice(Xc.indptr[i], Xc.indptr[i + 1])
                pos = {j: k for k, j in enumerate(Xc.indices[xs])}
                vals = np.array([Xc.data[xs][pos[j]] for j in take if j in pos],
                                dtype=np.float64)
            else:
                vals = np.asarray(Xc[i][take], dtype=np.float64)
            if vals.size:
                tot += float(vals.sum())
                n += vals.size
        label = f"{lo}" if lo == hi else (f"{lo}+" if hi > 10 ** 8 else f"{lo}-{hi}")
        out.append((label, n, tot / n if n else float("nan")))
    return out


def assert_counts_ladder(X, counts, source):
    """Refuse a matrix whose values do not rise with the counts behind them."""
    rungs = [r for r in count_ladder(X, counts)
             if r[1] >= MIN_STRATUM_ENTRIES and np.isfinite(r[2])]
    if len(rungs) < 3:
        raise BadMatrix(
            f"{source}: only {len(rungs)} count strata have at least "
            f"{MIN_STRATUM_ENTRIES} entries; too few to test the ladder")
    ratio = rungs[-1][2] / rungs[0][2]
    if ratio < COUNT_LADDER_MIN_RATIO:
        shape = "  ".join(f"{lab}:{m:.3f}" for lab, _, m in rungs)
        raise BadMatrix(
            f"{source}: the values barely rise with the counts behind them - "
            f"mean {rungs[0][2]:.3f} at {rungs[0][0]} count(s) and "
            f"{rungs[-1][2]:.3f} at {rungs[-1][0]}, a ratio of {ratio:.2f} "
            f"against a floor of {COUNT_LADDER_MIN_RATIO}. A correct log1p CP10K "
            f"matrix measures about 3.9 here. This is the signature of a second "
            f"normalisation, which preserves rank and destroys scale, so a "
            f"rank-based statistic will not show it and every mean will be "
            f"wrong. Ladder: {shape}. "
            f"Use a correctly normalised input.")
    return ratio


if __name__ == "__main__":
    import sys
    from matrix_check_selftest import _self_test
    sys.exit(_self_test() if "--self-test" in sys.argv else 0)
