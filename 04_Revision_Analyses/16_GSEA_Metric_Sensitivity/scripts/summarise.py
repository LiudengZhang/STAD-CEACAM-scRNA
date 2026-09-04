"""
Tables and verdicts from the sweep. Reads only this module's outputs plus the
published nfkb_per_celltype_*.csv; writes only into this module's outputs.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "outputs"
ROOT = HERE.parents[2]
PUB = {"sound13": ROOT / "02_New_Analyses" / "13_R1.8_Neutrophil_Rebuilt_Recompute"
       / "outputs" / "nfkb_per_celltype_sound13.csv",
       "live": ROOT / "02_New_Analyses" / "13_R1.8_Neutrophil_Rebuilt_Recompute"
       / "outputs" / "nfkb_per_celltype_live.csv"}
METRIC_ORDER = ["published", "signed_p", "logfc", "tstat"]


def load(stage, degset):
    sfx = "" if degset == "sound13" else f"_{degset}"
    f = OUT / f"nfkb_{stage}{sfx}.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f)
    d["degset"] = degset
    return d


def published(degset):
    p = pd.read_csv(PUB[degset])
    p = p[p["method"] == "ttest"].rename(columns={"cell_type": "cell"})
    return p[["cell", "phase", "nes", "nom_p", "fdr_q", "rank", "n_sets"]]


def repro(degset):
    m = load("metrics", degset)
    b = m[(m.metric == "published") & (m.seed == 42) & (m.permutations == 1000)]
    p = published(degset)
    j = b.merge(p, on=["cell", "phase"], suffixes=("_rerun", "_pub"))
    j["d_nes"] = j.nes_rerun - j.nes_pub
    j["d_q"] = j.fdr_q_rerun - j.fdr_q_pub
    j["rank_same"] = j.rank_rerun == j.rank_pub
    j["nsets_same"] = j.n_sets_rerun == j.n_sets_pub
    return j[["cell", "phase", "nes_pub", "nes_rerun", "d_nes", "fdr_q_pub",
              "fdr_q_rerun", "d_q", "rank_pub", "rank_rerun", "rank_same",
              "n_sets_pub", "n_sets_rerun", "nsets_same"]]


def wide(degset):
    m = load("metrics", degset)
    m = m[(m.seed == 42) & (m.permutations == 1000)]
    p = published(degset).rename(columns={c: f"pub_{c}" for c in
                                          ["nes", "nom_p", "fdr_q", "rank"]})
    piv = m.pivot_table(index=["cell", "phase"], columns="metric",
                        values=["nes", "nom_p", "fdr_q", "rank", "n_sets"])
    piv.columns = [f"{b}_{a}" for a, b in piv.columns]
    piv = piv.reset_index().merge(p, on=["cell", "phase"])
    return piv


def counted(degset):
    """The five counted claims, recomputed under every metric."""
    m = load("metrics", degset)
    m = m[(m.seed == 42) & (m.permutations == 1000)]
    rows = []
    for metric in METRIC_ORDER:
        s = m[m.metric == metric]
        r = dict(degset=degset, metric=metric)
        for ph in ("post", "pre"):
            x = s[s.phase == ph]
            pos = x[x.nes > 0]
            r[f"{ph}_n_positive"] = f"{len(pos)} of {len(x)}"
            r[f"{ph}_n_pos_q05"] = int((pos.fdr_q < 0.05).sum())
            r[f"{ph}_n_pos_q25"] = int((pos.fdr_q < 0.25).sum())
            r[f"{ph}_n_rank1"] = int((x["rank"] == 1).sum())
            r[f"{ph}_rank1_types"] = "; ".join(sorted(x.loc[x["rank"] == 1, "cell"]))
        mm = s[(s.cell == "MoMac") & (s.phase == "post")]
        r["MoMac_post_nes"] = round(float(mm.nes.iloc[0]), 3)
        r["MoMac_post_q"] = round(float(mm.fdr_q.iloc[0]), 4)
        r["MoMac_post_rank"] = int(mm["rank"].iloc[0])
        r["MoMac_post_nsets"] = int(mm["n_sets"].iloc[0])
        mp = s[(s.cell == "MoMac") & (s.phase == "pre")]
        r["MoMac_pre_nes"] = round(float(mp.nes.iloc[0]), 3)
        r["MoMac_pre_q"] = round(float(mp.fdr_q.iloc[0]), 4)
        r["MoMac_pre_rank"] = int(mp["rank"].iloc[0])
        r["MoMac_pre_nsets"] = int(mp["n_sets"].iloc[0])
        rows.append(r)
    return pd.DataFrame(rows)


def seed_spread(degset):
    """NES and rank across five seeds, same metric, same settings."""
    a = load("metrics", degset)
    b = load("seeds", degset)
    if b is None:
        return None
    a = a[(a.permutations == 1000) & (a.min_size == 15)]
    d = pd.concat([a, b], ignore_index=True)
    g = d.groupby(["cell", "phase", "metric"]).agg(
        n_seeds=("seed", "nunique"),
        nes_mean=("nes", "mean"), nes_sd=("nes", "std"),
        nes_min=("nes", "min"), nes_max=("nes", "max"),
        rank_min=("rank", "min"), rank_max=("rank", "max"),
        q_min=("fdr_q", "min"), q_max=("fdr_q", "max")).reset_index()
    g["nes_range"] = g.nes_max - g.nes_min
    g["rank_range"] = g.rank_max - g.rank_min
    return g


def settings(degset):
    s = load("settings", degset)
    if s is None:
        return None
    base = load("metrics", degset)
    base = base[(base.metric == "published") & (base.permutations == 1000)
                & (base.min_size == 15) & (base.max_size == 500)]
    d = pd.concat([base, s], ignore_index=True)
    d["setting"] = ("perm" + d.permutations.astype(str) + "_min"
                    + d.min_size.astype(str) + "_max" + d.max_size.astype(str))
    b = base.set_index(["cell", "phase"])
    d["d_nes_vs_base"] = d.apply(
        lambda r: r.nes - b.loc[(r.cell, r.phase), "nes"], axis=1)
    d["d_q_vs_base"] = d.apply(
        lambda r: r.fdr_q - b.loc[(r.cell, r.phase), "fdr_q"], axis=1)
    d["d_rank_vs_base"] = d.apply(
        lambda r: r["rank"] - b.loc[(r.cell, r.phase), "rank"], axis=1)
    return d


def momac_pre_distribution(degset):
    """Every Hallmark NES for MoMac pre, under every metric."""
    sfx = "" if degset == "sound13" else f"_{degset}"
    f = OUT / f"all_terms_metrics{sfx}.csv"
    d = pd.read_csv(f)
    d = d[(d.cell == "MoMac") & (d.phase == "pre") & (d.seed == 42)
          & (d.permutations == 1000)]
    return d.sort_values(["metric", "NES"], ascending=[True, False])


def main():
    for degset in ("sound13", "live"):
        if load("metrics", degset) is None:
            print(f"-- {degset}: not run yet")
            continue
        sfx = "" if degset == "sound13" else f"_{degset}"
        repro(degset).to_csv(OUT / f"reproduction_check{sfx}.csv", index=False)
        wide(degset).to_csv(OUT / f"nfkb_by_metric_wide{sfx}.csv", index=False)
        counted(degset).to_csv(OUT / f"counted_claims_by_metric{sfx}.csv",
                               index=False)
        momac_pre_distribution(degset).to_csv(
            OUT / f"momac_pre_nes_distribution{sfx}.csv", index=False)
        ss = settings(degset)
        if ss is not None:
            ss.to_csv(OUT / f"settings_effect{sfx}.csv", index=False)
        sp = seed_spread(degset)
        if sp is not None:
            sp.to_csv(OUT / f"seed_spread{sfx}.csv", index=False)
        print(f"-- wrote {degset} summaries")
    return 0


if __name__ == "__main__":
    sys.exit(main())
