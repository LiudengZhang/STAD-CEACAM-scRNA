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

Run with no arguments it counts the three Hallmark directories the shipped
tables were counted over, and writes them where the shipped tables sit:

    live     12_R1.8_DEG_Recompute/outputs/gsea          full_dataset.h5ad
    sound12  the archived sound-input recompute          twelve populations
    sound13  outputs/gsea_13types                        + rebuilt neutrophils

which reproduces outputs/counted_claims.csv and the three
outputs/nfkb_per_celltype_<label>.csv. claims.csv rows C061-C096 are checked
against nfkb_per_celltype_sound13.csv, so this is a reproduction step and not a
session tool. --gsea overrides the set.

Run:
    python count_nfkb_claims.py [--gsea LABEL=DIR ...] [--out DIR]
"""

from pathlib import Path
import argparse
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
MOD = HERE.parent
NFKB = HERE.parents[1] / "07_R1.8_NFkB_Specificity" / "scripts"
sys.path.insert(0, str(NFKB))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))

import nfkb_specificity as nf                                    # noqa: E402
from paths import SOUND_GSEA_DIR                                 # noqa: E402

PRIMARY = nf.PRIMARY

# The three sets the shipped tables were counted over. sound12 is named through
# paths.py rather than written out here: the directory it points at sits inside
# 07_Archive/, and an archive path written as a literal is both unresolvable in
# the deposit and invisible to the input checks.
DEFAULT_GSEA = [
    f"live={HERE.parents[1] / '12_R1.8_DEG_Recompute' / 'outputs' / 'gsea'}",
    f"sound12={SOUND_GSEA_DIR}",
    f"sound13={MOD / 'outputs' / 'gsea_13types'}",
]
DEFAULT_OUT = MOD / "outputs"


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
    ap.add_argument("--gsea", action="append", default=None,
                    metavar="LABEL=DIR")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    specs = args.gsea or DEFAULT_GSEA
    out_dir = args.out if args.out is not None else (
        None if args.gsea else str(DEFAULT_OUT))

    rows, tables = [], {}
    for spec in specs:
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
    if out_dir:
        o = Path(out_dir)
        o.mkdir(parents=True, exist_ok=True)
        df.to_csv(o / "counted_claims.csv", index=False)
        for label, g in tables.items():
            g.to_csv(o / f"nfkb_per_celltype_{label}.csv", index=False)
        print(f"\nwrote {o}/counted_claims.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
