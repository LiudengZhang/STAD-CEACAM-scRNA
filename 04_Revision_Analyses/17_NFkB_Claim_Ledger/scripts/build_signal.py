"""
signal.csv - one row per contrast per DE set: does this contrast carry any
gene-set signal at all, and how much of its TNFa/NF-kB number is permutation
noise?

Built from outputs/signal_runs.csv (364 prerank runs; see signal_scan.py).
The header block written at the top of signal.csv states the verdict criterion
in full, so the file is readable without this script.
"""
import sys
import numpy as np
import pandas as pd
from gsea_common import OUT, PUB_SEED, SEEDS

runs = pd.read_csv(OUT / "signal_runs.csv")

HEADER = f"""# signal.csv - the gate. 02_New_Analyses/17_NFkB_Claim_Ledger/
#
# One row per contrast (13 cell types x pre/post) per differential-expression
# set. Everything except the permutation seed is held at the published run's
# settings: ranking metric logFC x -log10(P), pinned
# 00_Reference/MSigDB_Hallmark_2020.gmt, permutation_num=1000, min_size=15,
# max_size=500, PYTHONHASHSEED=0. Seed {PUB_SEED} is the published seed; the noise
# columns are measured over seeds {SEEDS}.
#
# POSITIVE CONTROL (outputs/positive_control.csv). Before any negative result
# here is trusted, the same code was required to reproduce a published value.
# On all 26 `sound13` contrasts it reproduces
# 13_.../outputs/nfkb_per_celltype_sound13.csv to max |dNES| 3.4e-9, max |dq|
# 1.1e-16, with ranks identical 26/26 and set counts identical 26/26. A
# NEGATIVE control (the same machinery on a deliberately reversed ranked list)
# moves MoMac post from +2.228 to -2.218, i.e. the tool can fail when it should.
# On `live`, MoMac post measures 2.2249 against the published 2.2331 (dNES
# 0.008): the `live` tables were scored against a differently-ordered gene-set
# object, so their ES is reproducible and their permutation null is not.
# This is 16_GSEA_Metric_Sensitivity section 1's finding, reproduced here.
#
# VERDICT CRITERION - has_signal, stated once, applied mechanically:
#   yes  at the published seed, the contrast has at least one Hallmark set at
#        FDR q < 0.25 AND max |NES| over all sets >= 1.5, and both hold at 4 or
#        more of the 7 seeds.
#   weak exactly one of the two conditions holds at the published seed.
#   no   neither holds at the published seed.
# The threshold FDR<0.25 is GSEA's own conventional discovery threshold
# (Subramanian 2005); |NES| 1.5 is the level below which no set in this study
# is ever called significant. `no` means: this contrast does not distinguish
# any Hallmark set from the permutation null, so no NES read off it - the
# NF-kB set's included - measures anything.
#
# nfkb_quotable - a separate, stricter question about the NF-kB number itself:
#   yes  |NES| SD over the 7 seeds < 0.02, rank range <= 2, no degenerate
#        NES = 1.000 draw, and has_signal is yes.
#   no   otherwise. `no` does not mean the number is wrong; it means the
#        printed decimals are not reproducible from a different random seed.
#
# n_sets_same_es_sign is the size of the pool gseapy normalises against. When
# it is small the NES denominator is estimated from few permutation values and
# is unstable draw to draw; 16_ traced the MoMac-pre artefact to exactly this.
"""

rows = []
for (degset, cell, phase), g in runs.groupby(["degset", "cell", "phase"]):
    g = g.sort_values("seed")
    pub = g[g["seed"] == PUB_SEED].iloc[0]
    ok = (g["n_sets_q25"] >= 1) & (g["max_abs_nes"] >= 1.5)
    pub_ok_q = pub["n_sets_q25"] >= 1
    pub_ok_n = pub["max_abs_nes"] >= 1.5
    if pub_ok_q and pub_ok_n and int(ok.sum()) >= 4:
        verdict = "yes"
    elif pub_ok_q != pub_ok_n:
        verdict = "weak"
    else:
        verdict = "no"
    nes, q, p, rk = g["nfkb_nes"], g["nfkb_fdr_q"], g["nfkb_nom_p"], g["nfkb_rank"]
    quotable = ("yes" if (verdict == "yes" and nes.std(ddof=1) < 0.02
                          and (rk.max() - rk.min()) <= 2
                          and not g["nfkb_nes_degenerate"].any())
                else "no")
    rows.append(dict(
        degset=degset, cell_type=cell, contrast=phase,
        n_genes_ranked=int(pub["n_genes_ranked"]),
        n_hallmark_sets_tested=int(pub["n_sets_tested"]),
        tnfa_genes_in_ranked_list=int(pub["tnfa_overlap"]),
        max_abs_nes_any_set=round(float(pub["max_abs_nes"]), 4),
        n_sets_fdr25=int(pub["n_sets_q25"]), n_sets_fdr05=int(pub["n_sets_q05"]),
        nfkb_nes_seed42=round(float(pub["nfkb_nes"]), 4),
        nfkb_es_seed42=round(float(pub["nfkb_es"]), 4),
        nfkb_es_sign=pub["nfkb_es_sign"],
        n_sets_same_es_sign=int(pub["n_sets_same_es_sign"]),
        nfkb_fdr_q_seed42=round(float(pub["nfkb_fdr_q"]), 4),
        nfkb_rank_seed42=int(pub["nfkb_rank"]),
        nes_min=round(float(nes.min()), 4), nes_max=round(float(nes.max()), 4),
        nes_range=round(float(nes.max() - nes.min()), 4),
        nes_sd=round(float(nes.std(ddof=1)), 5),
        nom_p_min=round(float(p.min()), 4), nom_p_max=round(float(p.max()), 4),
        nom_p_sd=round(float(p.std(ddof=1)), 5),
        q_min=round(float(q.min()), 4), q_max=round(float(q.max()), 4),
        q_range=round(float(q.max() - q.min()), 4),
        q_sd=round(float(q.std(ddof=1)), 5),
        rank_min=int(rk.min()), rank_max=int(rk.max()),
        rank_range=int(rk.max() - rk.min()),
        rank_sd=round(float(rk.std(ddof=1)), 3),
        n_degenerate_seeds=int(g["nfkb_nes_degenerate"].sum()),
        n_seeds=len(g), has_signal=verdict, nfkb_quotable=quotable))

d = pd.DataFrame(rows).sort_values(["degset", "contrast", "cell_type"])
path = OUT / "signal.csv"
with open(path, "w") as fh:
    fh.write(HEADER)
    d.to_csv(fh, index=False)

for degset in ("live", "sound13"):
    s = d[d["degset"] == degset]
    print(f"\n=== {degset} ===")
    print(s[["cell_type", "contrast", "n_genes_ranked", "n_hallmark_sets_tested",
             "max_abs_nes_any_set", "n_sets_fdr25", "n_sets_fdr05",
             "nfkb_nes_seed42", "nfkb_es_seed42", "n_sets_same_es_sign",
             "nes_range", "q_range", "rank_min", "rank_max",
             "has_signal", "nfkb_quotable"]].to_string(index=False))
print("\nhas_signal counts:")
print(d.groupby(["degset", "contrast", "has_signal"]).size())
sys.exit(0)
