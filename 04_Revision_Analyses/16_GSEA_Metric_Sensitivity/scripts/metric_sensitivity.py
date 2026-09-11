"""
How much of the NF-kB result is the ranking metric, and how much is the data?

The published analysis (12_R1.8_DEG_Recompute/scripts/recompute_deg.py:343-351)
ranks genes for GSEA by

    logfoldchanges * -log10(clip(pvals, 1e-300))

That is a defensible choice, but it is not the field default and it weights the
p-value very heavily - a gene whose Welch p underflows to 0.0 is handed a
multiplier of 300, and 182 of MoMac's 4,693 post-treatment genes are in that
state. This module holds the differential expression completely fixed and
varies only the ranking metric and GSEA's own settings, so that any movement is
attributable to the scoring choice and not to the model or the cohort.

Nothing here is adopted. No manuscript file, figure, panel, provenance row or
verifier is touched; everything is written under this module's outputs/.

INPUTS - read only, never regenerated
    The thirteen-cell-type "sound" Welch t-test tables, which are the same
    tables `nfkb_per_celltype_sound13.csv` was built from:
      twelve types   07_Archive/2026-08-31_deg_recompute_on_sound_per_cell_type_inputs/
                     04_Revision_Analyses/12_R1.8_DEG_Recompute/outputs/deg/
      neutrophils    04_Revision_Analyses/13_R1.8_Neutrophil_Rebuilt_Recompute/
                     outputs/deg/
    (13_.../outputs/deg/ holds neutrophils only; merge_thirteen.sh shows the
    other twelve come from the archive, so the archive is the honest door to
    them. It is read, never written.)
    Gene sets: 00_Reference/MSigDB_Hallmark_2020.gmt, the file the published
    run used.

METRICS
    published   logFC * -log10(P)          the metric under test
    signed_p    sign(logFC) * -log10(P)    drops the effect size
    logfc       logFC                      drops the p-value
    tstat       the `scores` column        the Welch t statistic scanpy
                                           computed; a real t-like statistic,
                                           not a reconstruction

SETTINGS  permutation_num 1000 vs 10000; the min_size / max_size filters.
SEEDS     42 is the published seed. 1, 7, 13, 101 are added to measure the
          permutation noise the reproduction check is judged against.

Run: PYTHONHASHSEED=0 conda run -n stad_ceacam python metric_sensitivity.py
"""

from pathlib import Path
import argparse
import itertools
import os
import sys
import time

import numpy as np
import pandas as pd

os.environ.setdefault("PYTHONHASHSEED", "0")

HERE = Path(__file__).resolve().parent
MOD = HERE.parent
ROOT = MOD.parents[1]
OUT = MOD / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT / "00_Config"))
from paths import HALLMARK_GMT, SOUND_DEG_DIR                       # noqa: E402

HALLMARK = HALLMARK_GMT
# The twelve sound-input DEG tables, named through paths.py rather than as a
# literal into 07_Archive/. A literal there would be a live input read out of an
# archive - the fault 03_Final_Panels/verify_panel_provenance.py check 9
# exists to catch - and a path that resolves to nothing in the deposit. They are
# an intermediate the record carries as
# 02_Preparation_for_Panels/DEG_Sound_Recompute/deg.
SOUND12 = SOUND_DEG_DIR
NEUT = (ROOT / "04_Revision_Analyses" / "13_R1.8_Neutrophil_Rebuilt_Recompute"
        / "outputs" / "deg")
# The `live` tables - module 12's own outputs/, from full_dataset.h5ad. These
# are the ones verify_numbers.py checks the main text against, so the counted
# claims ("12 of 13", "q < 0.05 in four", "6 of 13") are counts over these and
# not over the sound thirteen. Both sets are swept, separately; within each the
# differential expression is held completely fixed.
LIVE = (ROOT / "04_Revision_Analyses" / "12_R1.8_DEG_Recompute" / "outputs" / "deg")
PUBLISHED = (ROOT / "04_Revision_Analyses" / "13_R1.8_Neutrophil_Rebuilt_Recompute"
             / "outputs" / "nfkb_per_celltype_sound13.csv")

CELLS = ["B_cells", "DC_cells", "Endothelial_cells", "Epithelial", "Fibroblast",
         "Mast_cells", "MoMac", "Neutrophils", "NK_cells", "Pericyte",
         "Plasma_cells", "TCD4_cells", "TCD8_cells"]
PHASES = ["pre", "post"]
NFKB_TERM = "TNF-alpha Signaling via NF-kB"
PUB_SEED = 42


# ----------------------------------------------------------------- metrics
def m_published(d):
    return d["logfoldchanges"] * -np.log10(d["pvals"].clip(lower=1e-300))


def m_signed_p(d):
    return np.sign(d["logfoldchanges"]) * -np.log10(d["pvals"].clip(lower=1e-300))


def m_logfc(d):
    return d["logfoldchanges"]


def m_tstat(d):
    # scanpy's rank_genes_groups(method="t-test_overestim_var") writes the Welch
    # t statistic into `scores`. It is present in every table, so the t-like
    # statistic is read, not reconstructed from logFC and P.
    return d["scores"]


METRICS = {"published": m_published, "signed_p": m_signed_p,
           "logfc": m_logfc, "tstat": m_tstat}


def deg_path(cell, phase, degset):
    if degset == "live":
        p = LIVE
    else:
        p = NEUT if cell == "Neutrophils" else SOUND12
    return p / f"{cell}_{phase}_ttest.csv"


def load_rank(cell, phase, metric, degset):
    d = pd.read_csv(deg_path(cell, phase, degset))
    d = d.dropna(subset=["logfoldchanges", "pvals", "scores"]).copy()
    d["metric"] = METRICS[metric](d)
    d = d[np.isfinite(d["metric"])]
    return (d.set_index("gene")["metric"].groupby(level=0).first()
            .sort_values(ascending=False))


def read_gmt():
    sets = {}
    for line in HALLMARK.read_text().splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            sets[f[0]] = {g for g in f[2:] if g}
    return sets


GMT = read_gmt()


def run_one(rnk, perm, min_size, max_size, seed, threads):
    import gseapy as gp
    res = gp.prerank(rnk=rnk, gene_sets=str(HALLMARK), permutation_num=perm,
                     min_size=min_size, max_size=max_size, seed=seed,
                     no_plot=True, outdir=None, threads=threads)
    out = res.res2d.copy()
    out["Term"] = out["Term"].astype(str).str.replace(r"^.*__", "", regex=True)
    out = out.sort_values("NES", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def summarise(tab, rnk, **meta):
    """The NF-kB row, plus the shape of the list it was scored against."""
    genes = set(rnk.index)
    overlap = len(GMT[NFKB_TERM] & genes)
    r = tab[tab["Term"] == NFKB_TERM]
    row = dict(meta)
    row.update(n_genes_ranked=len(rnk), n_sets=len(tab),
               tnfa_genes_in_gmt=len(GMT[NFKB_TERM]), tnfa_overlap=overlap,
               tnfa_overlap_pct=round(100 * overlap / len(GMT[NFKB_TERM]), 1),
               n_ties_at_max=int((rnk == rnk.max()).sum()),
               n_ties_at_min=int((rnk == rnk.min()).sum()),
               n_distinct_values=int(rnk.nunique()))
    if len(r):
        r = r.iloc[0]
        row.update(nes=float(r["NES"]), es=float(r["ES"]),
                   nom_p=float(r["NOM p-val"]), fdr_q=float(r["FDR q-val"]),
                   rank=int(r["rank"]), in_output=True)
    else:
        row.update(nes=np.nan, es=np.nan, nom_p=np.nan, fdr_q=np.nan,
                   rank=np.nan, in_output=False)
    return row


def main():
    ap = argparse.ArgumentParser()
    # No-argument default: every stage on both DE sets, which is exactly what
    # scripts/run_all.sh and scripts/run_live.sh between them ran and exactly
    # what outputs/ holds - nfkb_{metrics,settings,seeds}[_live].csv and the
    # matching all_terms_*. The figure driver launches every scripts/*.py with
    # no arguments, and `required=True` made that an argparse error rather than
    # a run. About 45 min at --threads 4.
    ap.add_argument("--stage", default="all",
                    choices=["all", "metrics", "settings", "seeds"])
    ap.add_argument("--degset", default="both",
                    choices=["both", "sound13", "live"])
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    stages = (["metrics", "settings", "seeds"] if args.stage == "all"
              else [args.stage])
    degsets = (["sound13", "live"] if args.degset == "both" else [args.degset])
    for stage in stages:
        for degset in degsets:
            one_sweep(stage, degset, args.threads)
    return 0


def one_sweep(stage, degset, threads):
    rows, dists = [], []
    t0 = time.time()

    if stage == "metrics":
        # A: four metrics, published GSEA settings, published seed.
        grid = [(m, PUB_SEED, 1000, 15, 500) for m in METRICS]
    elif stage == "settings":
        # B: published metric only; permutations, then the size filters.
        grid = [("published", PUB_SEED, 10000, 15, 500),
                ("published", PUB_SEED, 1000, 5, 500),
                ("published", PUB_SEED, 1000, 10, 500),
                ("published", PUB_SEED, 1000, 25, 500),
                ("published", PUB_SEED, 1000, 15, 200),
                ("published", PUB_SEED, 1000, 15, 1000)]
    else:
        # C: permutation noise, all four metrics, four extra seeds.
        grid = [(m, s, 1000, 15, 500)
                for m in METRICS for s in (1, 7, 13, 101)]

    for cell, phase in itertools.product(CELLS, PHASES):
        for metric, seed, perm, mn, mx in grid:
            rnk = load_rank(cell, phase, metric, degset)
            tag = f"{metric}_p{perm}_min{mn}_max{mx}_s{seed}"
            tab = run_one(rnk, perm, mn, mx, seed, threads)
            meta = dict(stage=stage, degset=degset, cell=cell,
                        phase=phase, metric=metric,
                        permutations=perm, min_size=mn, max_size=mx, seed=seed,
                        run_id=tag)
            rows.append(summarise(tab, rnk, **meta))
            d = tab[["Term", "NES", "NOM p-val", "FDR q-val", "rank"]].copy()
            for k, v in meta.items():
                d[k] = v
            dists.append(d)
            print(f"[{time.time()-t0:7.0f}s] {cell:<20}{phase:<5}{tag:<40}"
                  f"NES={rows[-1]['nes']}", flush=True)

    sfx = "" if degset == "sound13" else f"_{degset}"
    pd.DataFrame(rows).to_csv(OUT / f"nfkb_{stage}{sfx}.csv", index=False)
    pd.concat(dists, ignore_index=True).to_csv(
        OUT / f"all_terms_{stage}{sfx}.csv", index=False)
    print(f"{stage} {degset}: done in {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    sys.exit(main())
