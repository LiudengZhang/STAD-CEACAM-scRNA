"""
Positive control. Before any negative result from this module is trusted, the
tool must be shown able to reproduce a number that is already known and printed.

Three controls, in increasing strictness:

  C1  MoMac post, `sound13`, seed 42, published metric and settings, against
      13_.../outputs/nfkb_per_celltype_sound13.csv  (NES 2.228, q 0.000, rank 1)
  C2  the same for all 26 sound13 contrasts, NES / q / rank
  C3  MoMac post, `live`, against 07_.../outputs/nfkb_per_celltype.csv - the
      table verify_numbers.py checks the manuscript's 2.23 against.

A control that only ever returns "matches" proves nothing, so C4 is a
NEGATIVE control: the same machinery is pointed at a deliberately wrong input
(the ranked list reversed) and must NOT reproduce the published value.
"""
import sys
import numpy as np
import pandas as pd
from gsea_common import (load_rank, run_one, PUB_SOUND13, PUB_LIVE,
                         PUB_MANUSCRIPT, NFKB_TERM, CELLS, PHASES, OUT)

rows = []


def pub_row(path, cell, phase, method="ttest"):
    d = pd.read_csv(path)
    d = d[(d["method"] == method) & (d["phase"] == phase)
          & (d["cell_type"] == cell)]
    return None if d.empty else d.iloc[0]


def one(control, degset, cell, phase, pubpath, reverse=False):
    rnk = load_rank(cell, phase, degset)
    if reverse:
        rnk = (-rnk).sort_values(ascending=False)
    tab = run_one(rnk, 42)
    r = tab[tab["Term"] == NFKB_TERM].iloc[0]
    p = pub_row(pubpath, cell, phase)
    rows.append(dict(control=control, degset=degset, cell=cell, phase=phase,
                     reversed_input=reverse,
                     measured_nes=float(r["NES"]), published_nes=float(p["nes"]),
                     d_nes=abs(float(r["NES"]) - float(p["nes"])),
                     measured_q=float(r["FDR q-val"]), published_q=float(p["fdr_q"]),
                     d_q=abs(float(r["FDR q-val"]) - float(p["fdr_q"])),
                     measured_rank=int(r["rank"]), published_rank=int(p["rank"]),
                     measured_nsets=len(tab), published_nsets=int(p["n_sets"])))
    print(rows[-1], flush=True)


one("C1 MoMac post sound13", "sound13", "MoMac", "post", PUB_SOUND13)
for c in CELLS:
    for ph in PHASES:
        if (c, ph) == ("MoMac", "post"):
            continue
        one("C2 all sound13", "sound13", c, ph, PUB_SOUND13)
one("C3 MoMac post live vs the manuscript's own table", "live", "MoMac", "post",
    PUB_MANUSCRIPT)
one("C4 NEGATIVE: reversed ranked list", "sound13", "MoMac", "post",
    PUB_SOUND13, reverse=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "positive_control.csv", index=False)

c12 = df[df["control"].str.startswith(("C1", "C2"))]
print("\n--- summary ---")
print(f"C1+C2  max |dNES| = {c12['d_nes'].max():.3g}   "
      f"max |dq| = {c12['d_q'].max():.3g}   "
      f"ranks identical {int((c12['measured_rank']==c12['published_rank']).sum())}/{len(c12)}   "
      f"n_sets identical {int((c12['measured_nsets']==c12['published_nsets']).sum())}/{len(c12)}")
c3 = df[df["control"].str.startswith("C3")].iloc[0]
print(f"C3     dNES = {c3['d_nes']:.3g}  dq = {c3['d_q']:.3g}  "
      f"rank {c3['measured_rank']} vs {c3['published_rank']}")
c4 = df[df["control"].str.startswith("C4")].iloc[0]
print(f"C4     reversed-input dNES = {c4['d_nes']:.3g}  "
      f"(must be large; {'PASS' if c4['d_nes'] > 1 else 'FAIL'})")
sys.exit(0)
