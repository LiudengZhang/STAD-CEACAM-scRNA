MSigDB_Hallmark_2020.gmt
========================

Downloaded from Enrichr with gseapy 1.1.4 (gseapy.get_library) on 2026-08-29.
50 gene sets, 7321 gene entries, md5 cbd26bdb22d1f84cd6452bd2f86eac48.

Why this file exists
--------------------
recompute_deg.py used to pass the library NAME to gseapy, which fetches it from
Enrichr at run time. That is an undeclared dependency on a remote service: the
analysis could not be repeated offline, and nothing recorded which release had
been used. Pinning removes that.

What this file is NOT the explanation for
-----------------------------------------
It was first pinned on the belief that Enrichr had updated the library between
two runs. That belief was wrong and is withdrawn. This file was checked against
the live library on the day it was written and the two are content-identical,
set for set and gene for gene.

The difference that prompted the investigation - the pre-treatment B cells
TNFa/NF-kB score reading NES 1.183727 on 2026-08-28 and 1.173631 on 2026-08-29
from the same differential-expression table - is Monte Carlo error in the NES
normalisation, not a change of gene sets. See 00_Data_Audit/FINDINGS.md,
section 10.2, for the evidence.

Do not replace this file without re-running
02_New_Analyses/12_R1.8_DEG_Recompute/scripts/recompute_deg.py and re-checking
every enrichment number the manuscript quotes with
04_Manuscript_R1/verify_numbers.py.
