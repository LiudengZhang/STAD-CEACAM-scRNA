#!/usr/bin/env python3
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
"""Rewrite the original graphst_config.yaml for the rsrch8 reproduction.

Positional, identifier-free textual substitution, so the diff against the
original is readable and every changed line is recorded.  NOTHING but input
paths, the device and results_dir moves; every analysis parameter is untouched.

usage: make_config.py <original_yaml> <out_yaml> <abs_results_dir>
"""
import re, sys, hashlib, pathlib

ORIG, OUT, RESULTS_DIR = sys.argv[1], sys.argv[2], sys.argv[3]
ROOT = "/path/to/home/graphst_repro"
N_SAMPLES = 10

txt = pathlib.Path(ORIG).read_text()
changes = []

def sub1(pat, rep, why):
    """Substitute exactly one occurrence; refuse if the count is not 1."""
    global txt
    new, n = re.subn(pat, rep, txt, flags=re.M)
    if n != 1:
        raise SystemExit(f"REFUSING: expected 1 match, got {n} for {why}")
    changes.append(why)
    txt = new

# 1. reference h5ad -> staged copy (matched by filename, not by path)
sub1(r'^  path: ".*stomach_14types_reference\.h5ad"$',
     f'  path: "{ROOT}/data/reference/stomach_14types_reference.h5ad"',
     "reference.path -> staged copy")

# 2/3. the ten Visium sample paths and their (unused) sample_id fields, taken
#      positionally in the order the original config lists them.
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

# 5. per-run absolute output dir (the runner resolves this against the config dir)
sub1(r'^  results_dir: ".*"$', f'  results_dir: "{RESULTS_DIR}"',
     "output.results_dir -> per-run absolute dir")

pathlib.Path(OUT).write_text(txt)
print(f"wrote {OUT}  md5={hashlib.md5(txt.encode()).hexdigest()}")
print(f"  {len(changes)} lines changed")
