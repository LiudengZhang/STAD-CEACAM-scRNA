"""
Pin the Hallmark gene sets to a file, so the analysis does not depend on Enrichr.

recompute_deg.py passed the library NAME "MSigDB_Hallmark_2020" to gseapy, which
fetches it from Enrichr at run time. That is an undeclared dependency on a remote
service: the analysis cannot be repeated offline, and nothing records which
release was used. This writes ../../../00_Reference/MSigDB_Hallmark_2020.gmt and
a note with the date and checksum; recompute_deg.py reads that file.

Pinning does not change any result. The pinned copy was compared against the
live library on the day it was written and is content-identical, set for set and
gene for gene. An earlier note in this file claimed the library had changed
between two runs and that this explained a shift in the reported scores; that
claim was wrong and is withdrawn. The shift is Monte Carlo error in the NES
normalisation - see 00_Data_Audit/FINDINGS.md, section 10.2.

Re-run this only to deliberately adopt a newer release of the library, and re-run
recompute_deg.py after it.

Run: python pin_gene_sets.py
"""

from pathlib import Path
import datetime
import hashlib
import sys

import gseapy as gp

LIBRARY = "MSigDB_Hallmark_2020"
REF = Path(__file__).resolve().parents[1] / "reference"


def main():
    REF.mkdir(parents=True, exist_ok=True)
    gmt = REF / f"{LIBRARY}.gmt"
    if gmt.exists():
        print(f"{gmt.name} already exists; delete it first to re-pin.")
        print(f"  md5 {hashlib.md5(gmt.read_bytes()).hexdigest()}")
        return

    lib = gp.get_library(LIBRARY)
    if not lib:
        sys.exit(f"Enrichr returned nothing for {LIBRARY}")
    with open(gmt, "w") as fh:
        for term in sorted(lib):
            fh.write("\t".join([term, ""] + sorted(lib[term])) + "\n")

    digest = hashlib.md5(gmt.read_bytes()).hexdigest()
    today = datetime.date.today().isoformat()
    n_genes = sum(len(v) for v in lib.values())
    (REF / "README.txt").write_text(
        f"{LIBRARY}.gmt\n"
        f"{'=' * (len(LIBRARY) + 4)}\n\n"
        f"Downloaded from Enrichr with gseapy {gp.__version__} "
        f"(gseapy.get_library) on {today}.\n"
        f"{len(lib)} gene sets, {n_genes} gene entries, md5 {digest}.\n\n"
        "Why this file exists\n"
        "--------------------\n"
        "recompute_deg.py used to pass the library name to gseapy, which fetches\n"
        "it from Enrichr at run time. That is an undeclared dependency on a\n"
        "remote service: the analysis could not be repeated offline, and nothing\n"
        "recorded which release had been used. Pinning removes that, and changes\n"
        "no result - the pinned copy was checked against the live library and is\n"
        "content-identical.\n\n"
        "Do not replace this file without re-running recompute_deg.py and\n"
        "re-checking every enrichment number the manuscript quotes with\n"
        "04_Manuscript_R1/verify_numbers.py.\n")
    print(f"terms {len(lib)}   gene entries {n_genes}")
    print(f"md5   {digest}")
    print(f"saved {gmt}")


if __name__ == "__main__":
    main()
