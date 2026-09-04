"""
Rescore the existing MAST differential-expression tables against the pinned
Hallmark gene sets, so the whole gsea/ directory comes from one library.

Why this is a separate script. MAST needs rpy2 to open R, which fails on this
machine under the driver's inherited LD_LIBRARY_PATH, so recompute_deg.py
records mast_status and moves on. Its MAST *differential-expression* tables from
the 2026-08-28 run are still valid - nothing about them depends on the gene sets
- but their enrichment scores were computed with the library fetched from
Enrichr rather than from the pinned copy. The two are content-identical, so
this rescoring is bookkeeping: it makes every file in gsea/ trace to the same
declared input.

This reads deg/*_mast.csv and rewrites gsea/*_mast_hallmark.csv using
HALLMARK_GMT, with the same ranking metric, seed and size filters as
recompute_deg.py. It does not touch the t-test files and does not need R.

Run: python rescore_mast_gsea.py
"""

from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))

from recompute_deg import OUT, gsea  # noqa: E402


def main():
    tables = sorted((OUT / "deg").glob("*_mast.csv"))
    if not tables:
        sys.exit(f"no MAST tables in {OUT / 'deg'}")
    changed, unchanged, skipped = 0, 0, []
    for deg in tables:
        tag = deg.stem                      # e.g. MoMac_post_mast
        out = OUT / "gsea" / f"{tag}_hallmark.csv"
        df = pd.read_csv(deg)
        res = gsea(df, tag)
        if res is None:
            skipped.append(tag)
            continue
        before = None
        if out.exists():
            old = pd.read_csv(out).set_index("Term")
            term = [t for t in old.index if "TNF-alpha" in str(t)]
            before = float(old.loc[term[0], "NES"]) if term else None
        res.to_csv(out, index=False)
        after = res.set_index("Term")
        term = [t for t in after.index if "TNF-alpha" in str(t)]
        now = float(after.loc[term[0], "NES"]) if term else None
        if before is not None and now is not None and abs(before - now) > 1e-9:
            changed += 1
            print(f"  {tag:<28} NF-kB NES {before:+.4f} -> {now:+.4f}")
        else:
            unchanged += 1
    print(f"\n{len(tables)} MAST tables rescored against {Path(__file__).name}'s "
          f"pinned library: {changed} moved, {unchanged} unchanged")
    if skipped:
        print(f"  too few ranked genes, no enrichment written: {', '.join(skipped)}")


if __name__ == "__main__":
    main()
