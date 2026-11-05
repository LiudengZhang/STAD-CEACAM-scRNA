#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
"""Rewrite the original graphst_config.yaml for the rsrch8 reproduction.

Positional, identifier-free textual substitution, so the diff against the
original is readable and every changed line is recorded.  NOTHING but input
paths, the device and results_dir moves; every analysis parameter is untouched.
"""
import re, sys, hashlib, pathlib

ORIG, OUT, RESULTS_DIR = sys.argv[1], sys.argv[2], sys.argv[3]
ROOT = "/path/to/home/graphst_repro"
N_SAMPLES = 10

txt = pathlib.Path(ORIG).read_text()
changes = []

def sub1(pat, rep, why):
    global txt
    new, n = re.subn(pat, rep, txt, flags=re.M)
    if n != 1:
        raise SystemExit(f"REFUSING: expected 1 match, got {n} for {why}")
    changes.append(why)
    txt = new

sub1(r'^  path: ".*stomach_14types_reference\.h5ad"$',
     f'  path: "{ROOT}/data/reference/stomach_14types_reference.h5ad"',
     "reference.path -> staged copy")

for i in range(1, N_SAMPLES + 1):
    s = f"sample_{i:02d}"
    def nth(pat, rep, why, _i=i):
        global txt
        hits = list(re.finditer(pat, txt, flags=re.M))
        if len(hits) != N_SAMPLES - (_i - 1):
            raise SystemExit(f"REFUSING: {why}: expected {N_SAMPLES-(_i-1)} remaining, got {len(hits)}")
        m = hits[0]
        txt = txt[:m.start()] + rep + txt[m.end():]
        changes.append(why)
    nth(r'^    path: "(?!/path/to/machine/).*"$', f'    path: "{ROOT}/data/spatial/{s}"',
        f"{s}.path -> staged copy")
    nth(r'^    sample_id: "(?!sample_).*"$', f'    sample_id: "{s}"',
        f"{s}.sample_id -> de-identified (field is never read by the runner)")

# 4. device pinned explicitly to cpu (the original auto-detected, found no GPU,
#    and fell back to cpu -- see the original run log)
sub1(r'^  device: "auto".*$',
     '  device: "cpu"              # PINNED explicitly (original auto-detected -> cpu)',
     "device auto -> explicit cpu")

sub1(r'^  results_dir: ".*"$', f'  results_dir: "{RESULTS_DIR}"',
     "output.results_dir -> per-run absolute dir")

pathlib.Path(OUT).write_text(txt)
