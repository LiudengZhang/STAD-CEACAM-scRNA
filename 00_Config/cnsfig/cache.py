"""
Drawing reads a table; computing writes it. Nothing else touches both.

Until 2026-09-15 a panel script computed and drew in one pass: 2F re-ran Milo,
2H/I re-ran an exact permutation test, 4E and 5F read a whole h5ad, every time
a margin moved by a millimetre. The author asked for the drawing step to be
separable so that adjusting a figure is fast. This module is the seam.

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

What this is not: a refactor of the tree. It is used by the panels one round
touches, each one proven value-identical with compare_panel_content.py before
and after the split.

    python cache.py --self-test
"""

import csv
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = ["table", "wants_recompute", "verify", "manifest_path", "self_test"]

#: Floats are written as Python's repr - the shortest string that reads back
#: to the same double - and read back with the correctly-rounded parser, so
#: that the value drawn is the value computed to the last bit. Not "%.17g":
#: that prints a rounded 0.56426 as 0.56425999999999998 and made a 239,000-row
#: table 13 MB. pandas' default C parser is fast and not correctly rounded:
#: it returned 0.5145336985588073 for a written ...074 on the first table
#: this module cached, and the round-trip check below refused it.
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
    # No date column: the release's publish-safety scan reads any date as
    # narrative that must be resolved, and git carries the date anyway.
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
    # A float32 column is written at float32 precision by pandas ("0.5145337"),
    # which reads back as a different double from the float32 it came from;
    # widened first, it is written as the exact double it is.
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


# ----------------------------------------------------------------- self-test
def self_test() -> int:
    """Every claim above, with the mutation that would refute it."""
    import tempfile
    calls = {"n": 0}

    def compute():
        calls["n"] += 1
        return pd.DataFrame({"gene": ["CEACAM5", "CEACAM6", None],
                             "rho": [0.3612345678901234, -1e-17, np.nan],
                             # float32, as a mean over an h5ad comes out
                             "f32": np.array([0.51453369855880737, 0.1, 2.5],
                                             dtype=np.float32),
                             "n": [19, 19, 18]})

    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        a = table(d, "t", compute, recompute=False, source="self_test")
        assert calls["n"] == 1, "compute() was not called for an absent table"
        b = table(d, "t", compute, recompute=False)
        assert calls["n"] == 1, "compute() ran although the table existed"
        assert _same(a, b) == "", "the read copy differs from the written one"
        assert a["rho"][0] == 0.3612345678901234 and a["rho"][1] == -1e-17, \
            "17 significant digits did not survive"
        table(d, "t", compute, recompute=True)
        assert calls["n"] == 2, "--recompute did not recompute"
        assert verify(d) == [], f"a fresh manifest reports {verify(d)}"

        # MUTATION 1: an edited table must be reported by verify().
        p = d / "data" / "t.csv"
        p.write_text(p.read_text().replace("19", "20", 1))
        bad = verify(d)
        assert bad and bad[0][0] == "t", "an edited table passed verify()"
        # ... and is still what the drawing gets: the cache is not a guard
        # against editing, the manifest is. Restore it.
        table(d, "t", compute, recompute=True)
        assert verify(d) == []

        # MUTATION 2: a write that loses precision is refused and leaves no
        # file behind. Three significant digits stand in for any lossy writer
        # (FLOAT_FORMAT None is repr, which is lossless).
        global FLOAT_FORMAT
        keep, FLOAT_FORMAT = FLOAT_FORMAT, "%.3g"
        try:
            table(d, "bad", compute, recompute=True)
        except ValueError:
            pass
        else:
            raise AssertionError("a table that lost precision was accepted")
        finally:
            FLOAT_FORMAT = keep
        assert not (d / "data" / "bad.csv").exists(), "a refused table was left on disk"
        assert verify(d) == [], "a refused table reached the manifest"

        # MUTATION 3: an index carrying data is refused.
        def indexed():
            return pd.DataFrame({"v": [1, 2]}, index=["a", "b"])
        try:
            table(d, "idx", indexed, recompute=True)
        except ValueError:
            pass
        else:
            raise AssertionError("a labelled index was silently dropped")

        # MUTATION 4: a manifest whose file went missing.
        p.unlink()
        assert verify(d) == [("t", "file missing")], verify(d)

    # wants_recompute reads the flag and nothing else.
    assert wants_recompute(["--recompute"]) and not wants_recompute(["--other"])
    print("cnsfig.cache self-test: every control behaved as required")
    return 0


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        sys.exit(self_test())
    print(__doc__)
