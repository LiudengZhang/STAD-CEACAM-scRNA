"""
The counted NF-kB claims, from any one set of Hallmark tables.

Every number the main text prints about "every population exhibited elevated
NF-kB pathway activity" is a count over one directory of
12_R1.8_DEG_Recompute/outputs/gsea/*_hallmark.csv. This reads a directory,
counts them, and prints them; it computes nothing that
07_R1.8_NFkB_Specificity/scripts/nfkb_specificity.py does not, and it reuses
that script's own `load_gsea` rather than reimplementing the parsing, so a
disagreement here would be a disagreement with the shipped analysis.

Nothing is written outside the directory given with --out, and no panel is
drawn.

Run:
    python count_nfkb_claims.py --gsea <dir> [--gsea <dir> ...] [--out DIR]
"""

from pathlib import Path
import argparse
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
NFKB = HERE.parents[1] / "07_R1.8_NFkB_Specificity" / "scripts"
sys.path.insert(0, str(NFKB))

import nfkb_specificity as nf                                    # noqa: E402

PRIMARY = nf.PRIMARY


def claims(g, label):
    """The counted claims, under the primary test."""
    p = g[g["method"] == PRIMARY]
    out = dict(label=label)
    for phase in ("pre", "post"):
        s = p[p["phase"] == phase]
        pos = s[s["nes"] > 0]
        out[f"{phase}_n_types"] = len(s)
        out[f"{phase}_n_positive"] = len(pos)
        out[f"{phase}_n_positive_q05"] = int((pos["fdr_q"] < 0.05).sum())
        out[f"{phase}_types_q05"] = "; ".join(
            sorted(pos.loc[pos["fdr_q"] < 0.05, "cell_type"]))
        out[f"{phase}_n_rank1"] = int((s["rank"] == 1).sum())
        out[f"{phase}_types_rank1"] = "; ".join(
            sorted(s.loc[s["rank"] == 1, "cell_type"]))

    def one(cell, phase, field):
        r = p[(p["cell_type"] == cell) & (p["phase"] == phase)]
        return None if r.empty else r.iloc[0][field]

    out["MoMac_pre_nes"] = one("MoMac", "pre", "nes")
    out["MoMac_pre_q"] = one("MoMac", "pre", "fdr_q")
    out["B_cells_post_rank"] = one("B_cells", "post", "rank")
    out["B_cells_post_nsets"] = one("B_cells", "post", "n_sets")
    out["Pericyte_post_rank"] = one("Pericyte", "post", "rank")
    out["Pericyte_post_nsets"] = one("Pericyte", "post", "n_sets")
    out["Epithelial_post_nes"] = one("Epithelial", "post", "nes")
    out["Neutrophils_post_nes"] = one("Neutrophils", "post", "nes")
    out["Neutrophils_post_q"] = one("Neutrophils", "post", "fdr_q")
    out["Neutrophils_post_rank"] = one("Neutrophils", "post", "rank")
    out["Neutrophils_post_nsets"] = one("Neutrophils", "post", "n_sets")
    out["Neutrophils_pre_nes"] = one("Neutrophils", "pre", "nes")
    out["Neutrophils_pre_q"] = one("Neutrophils", "pre", "fdr_q")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gsea", action="append", required=True,
                    metavar="LABEL=DIR")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rows, tables = [], {}
    for spec in args.gsea:
        label, d = spec.split("=", 1)
        nf.RECOMPUTE_GSEA = Path(d)
        g = nf.load_gsea()
        tables[label] = g
        rows.append(claims(g, label))
        print(f"\n=== {label}  ({d}) ===")
        for phase in ("pre", "post"):
            s = (g[(g["method"] == PRIMARY) & (g["phase"] == phase)]
                 .sort_values("nes", ascending=False))
            print(f"  [{phase}]  {'cell type':<22}{'NES':>8}{'FDR q':>10}"
                  f"{'rank':>8}")
            for _, r in s.iterrows():
                print(f"        {r['cell_type']:<24}{r['nes']:>8.3f}"
                      f"{r['fdr_q']:>10.4f}{r['rank']:>6} / {r['n_sets']}")
    df = pd.DataFrame(rows)
    print("\n" + df.to_string(index=False))
    if args.out:
        o = Path(args.out)
        o.mkdir(parents=True, exist_ok=True)
        df.to_csv(o / "counted_claims.csv", index=False)
        for label, g in tables.items():
            g.to_csv(o / f"nfkb_per_celltype_{label}.csv", index=False)
        print(f"\nwrote {o}/counted_claims.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
