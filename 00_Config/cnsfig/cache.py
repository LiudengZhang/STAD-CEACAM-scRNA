"""
Drawing reads a table; computing writes it. Nothing else touches both.

A panel script that computes and draws in one pass re-runs Milo, or an exact
permutation test, or a whole h5ad read, every time a margin moves by a
millimetre. This module is the seam that makes the drawing step separable.

    from cnsfig.cache import table

    pts = table(PANEL_DIR, "points", lambda: compute_points(adata_path))
    ax.scatter(pts.x, pts.y)

`table()` returns `PANEL_DIR/data/<name>.csv` if it exists, and otherwise
calls `compute()`, writes what it returns, reads it back, checks the round trip
value for value, and returns the read copy - so the drawing is made from the
file in every case, never from the in-memory result. Pass `--recompute` on the
command line (or `recompute=True`) to run `compute()` regardless.

Every write is recorded in `data/MANIFEST.csv` with its md5, so
`freeze_baseline.py` can hold the tables the way it holds the analysis outputs,
and `verify(PANEL_DIR)` says whether any table has been edited since it was
written. Rule 3 becomes checkable rather than asserted: a script whose drawing
half has no data path but `data/` cannot move a number.

What this is not: a refactor of the tree. A panel is moved onto the seam only
with compare_panel_content.py's word that not one value moved.

    python cache.py --self-test
"""

import csv
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = ["table", "wants_recompute", "verify", "manifest_path"]

# Floats are written as repr (shortest exact form) and read back with the
# round-trip parser, so the value drawn is the value computed to the last bit.
# pandas' default C parser is not correctly rounded.
FLOAT_FORMAT = None
READ_KW = {"float_precision": "round_trip"}


def wants_recompute(argv=None) -> bool:
    return "--recompute" in (sys.argv[1:] if argv is None else argv)


def _md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def manifest_path(panel_dir) -> Path:
    return Path(panel_dir) / "data" / "MANIFEST.csv"


def _record(panel_dir: Path, name: str, path: Path, df: pd.DataFrame, source: str):
    mp = manifest_path(panel_dir)
    rows = []
    if mp.exists():
        with mp.open(newline="") as fh:
            rows = [r for r in csv.DictReader(fh) if r["table"] != name]
    rows.append({"table": name, "file": path.name, "md5": _md5(path),
                 "rows": len(df), "columns": len(df.columns),
                 "source_script": source})
    rows.sort(key=lambda r: r["table"])
    with mp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["table", "file", "md5", "rows",
                                           "columns", "source_script"])
        w.writeheader()
        w.writerows(rows)


def _same(a: pd.DataFrame, b: pd.DataFrame) -> str:
    """'' if b is a's round trip through CSV, else what differs."""
    if list(a.columns) != list(b.columns):
        return f"columns {list(a.columns)} -> {list(b.columns)}"
    if len(a) != len(b):
        return f"rows {len(a)} -> {len(b)}"
    for c in a.columns:
        x, y = a[c].to_numpy(), b[c].to_numpy()
        if np.issubdtype(np.asarray(x).dtype, np.number) or \
                np.issubdtype(np.asarray(y).dtype, np.number):
            x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float)
            y = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(float)
            bad = ~(np.isclose(x, y, rtol=0, atol=0, equal_nan=True))
            if bad.any():
                i = int(np.argmax(bad))
                return f"column {c!r} row {i}: {x[i]!r} -> {y[i]!r}"
        else:
            xs = pd.Series(x).astype(str).where(pd.notna(pd.Series(x)), "")
            ys = pd.Series(y).astype(str).where(pd.notna(pd.Series(y)), "")
            bad = (xs != ys).to_numpy()
            if bad.any():
                i = int(np.argmax(bad))
                return f"column {c!r} row {i}: {x[i]!r} -> {y[i]!r}"
    return ""


def table(panel_dir, name, compute, *, recompute=None, source=None) -> pd.DataFrame:
    """The table `data/<name>.csv`, computed only if absent or asked for.

    `compute` returns a DataFrame with a default RangeIndex; anything that
    matters must be a column, because the index is not written.
    """
    panel_dir = Path(panel_dir)
    path = panel_dir / "data" / f"{name}.csv"
    if recompute is None:
        recompute = wants_recompute()
    if path.exists() and not recompute:
        return pd.read_csv(path, **READ_KW)

    df = compute()
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"{name}: compute() must return a DataFrame, "
                        f"got {type(df).__name__}")
    if not isinstance(df.index, pd.RangeIndex):
        raise ValueError(f"{name}: the index is not written; reset_index() "
                         f"so that it becomes a column")
    # float32 columns are widened so they are written as exact doubles.
    df = df.copy()
    for c in df.columns:
        if df[c].dtype == np.float32:
            df[c] = df[c].astype(np.float64)
    path.parent.mkdir(exist_ok=True)
    df.to_csv(path, index=False, float_format=FLOAT_FORMAT)
    back = pd.read_csv(path, **READ_KW)
    diff = _same(df, back)
    if diff:
        path.unlink()
        raise ValueError(f"{name}: does not survive the CSV round trip - {diff}")
    src = source or (Path(sys.argv[0]).name if sys.argv and sys.argv[0] else "")
    _record(panel_dir, name, path, back, src)
    print(f"  cache: wrote {path.relative_to(panel_dir)} "
          f"({len(back)} rows x {len(back.columns)} columns)")
    return back


def verify(panel_dir) -> list:
    """Tables whose md5 no longer matches the manifest: [(table, why)]."""
    mp = manifest_path(panel_dir)
    if not mp.exists():
        return [("MANIFEST.csv", "absent")]
    out = []
    with mp.open(newline="") as fh:
        for r in csv.DictReader(fh):
            p = mp.parent / r["file"]
            if not p.exists():
                out.append((r["table"], "file missing"))
            elif _md5(p) != r["md5"]:
                out.append((r["table"], "md5 differs from the manifest"))
    return out


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        from cache_selftest import self_test
        sys.exit(self_test())
    print(__doc__)
