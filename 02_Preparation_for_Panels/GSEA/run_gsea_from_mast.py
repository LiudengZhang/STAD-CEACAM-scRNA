#!/usr/bin/env python3
"""
Rebuild full Hallmark GSEA from MAST DEGs for all 13 cell types × pre+post.

Parameters:
  - Ranking: logFC × -log10(pval)
  - Gene sets: MSigDB_Hallmark_2020
  - min_size=5, max_size=500, permutation_num=1000, seed=42
"""

import pandas as pd
import numpy as np
import gseapy as gp
from pathlib import Path
import sys
import shutil

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from paths import DEG_PRE_DIR, DEG_POST_DIR, GSEA_DIR, DEG_CELL_TYPES

GSEA_PRE_DIR = GSEA_DIR / "pre"
GSEA_POST_DIR = GSEA_DIR / "post"
ARCHIVED_DIR = GSEA_DIR / "_archived"

NFKB_TERM = "TNF-alpha Signaling via NF-kB"


def compute_ranking(deg_csv: Path) -> pd.Series:
    df = pd.read_csv(deg_csv)
    pvals_clipped = df['pvals'].clip(lower=1e-300)
    # Negate logFC: MAST DEGs are R-vs-NR; flip to NR-vs-R convention
    df['rank_metric'] = -df['logfoldchanges'] * (-np.log10(pvals_clipped))
    rnk = df.set_index('gene')['rank_metric'].dropna()
    rnk = rnk[~rnk.index.duplicated(keep='first')].sort_values(ascending=False)
    return rnk


def run_gsea_prerank(rnk: pd.Series, outdir: Path) -> pd.DataFrame:
    try:
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets='MSigDB_Hallmark_2020',
            outdir=str(outdir),
            min_size=5,
            max_size=500,
            permutation_num=1000,
            seed=42,
            verbose=False,
        )
        return pre_res.res2d
    except LookupError:
        print("    GSEA failed: no gene sets passed filtering (too few genes)")
        return None


def main():
    print("=" * 60)
    print("Rebuilding GSEA from MAST DEGs — full Hallmark set")
    print("=" * 60)

    # Archive old files (if any remain)
    ARCHIVED_DIR.mkdir(parents=True, exist_ok=True)
    old_files = list(GSEA_DIR.glob("*_mast_prerank_gsea.csv"))
    if old_files:
        for f in old_files:
            shutil.move(str(f), str(ARCHIVED_DIR / f.name))
            print(f"Archived: {f.name}")
    else:
        print("No old files to archive (already moved)")

    # Create output dirs
    GSEA_PRE_DIR.mkdir(parents=True, exist_ok=True)
    GSEA_POST_DIR.mkdir(parents=True, exist_ok=True)

    # Track NF-kB NES for rankings CSV
    nfkb_records = []

    for ct in DEG_CELL_TYPES:
        print(f"\n--- {ct} ---")

        for comparison, deg_dir, gsea_out_dir in [
            ("pre", DEG_PRE_DIR, GSEA_PRE_DIR),
            ("post", DEG_POST_DIR, GSEA_POST_DIR),
        ]:
            deg_csv = deg_dir / f"{ct}_mast_deg.csv"
            if not deg_csv.exists():
                print(f"  {comparison}: DEG CSV not found, skipping")
                continue

            rnk = compute_ranking(deg_csv)
            print(f"  {comparison}: {len(rnk)} genes ranked")

            if len(rnk) < 15:
                print(f"  {comparison}: too few genes, skipping")
                continue

            tmp_dir = GSEA_DIR / f"_tmp_{ct}_{comparison}"
            res2d = run_gsea_prerank(rnk, tmp_dir)
            if res2d is None:
                nfkb_records.append({'cell_type': ct, 'comparison': comparison, 'nes': np.nan})
                shutil.rmtree(tmp_dir, ignore_errors=True)
                continue

            # Save full results
            out_csv = gsea_out_dir / f"{ct}_gsea_hallmark.csv"
            res2d.to_csv(out_csv, index=False)
            print(f"  {comparison}: saved {len(res2d)} terms → {out_csv.name}")

            # Extract NF-kB NES
            nfkb_row = res2d[res2d['Term'] == NFKB_TERM]
            if len(nfkb_row) > 0:
                nes_val = nfkb_row['NES'].values[0]
                print(f"  {comparison}: NF-kB NES = {nes_val:.4f}")
            else:
                nes_val = np.nan
                print(f"  {comparison}: NF-kB term not found")

            nfkb_records.append({
                'cell_type': ct,
                'comparison': comparison,
                'nes': nes_val,
            })

            # Clean up gseapy temp dir
            shutil.rmtree(tmp_dir, ignore_errors=True)

    # Build nfkb_rankings CSV
    print("\n" + "=" * 60)
    print("Building nfkb_rankings_13types_mast.csv")
    print("=" * 60)

    LABELS = {
        "B_cells": "B cells", "DC_cells": "DC", "Endothelial_cells": "Endothelial",
        "Epithelial": "Epithelial", "Fibroblast": "Fibroblast", "Mast_cells": "Mast cells",
        "MoMac": "MoMac", "Neutrophils": "Neutrophils", "NK_cells": "NK cells",
        "Pericyte": "Pericyte", "Plasma_cells": "Plasma cells",
        "TCD4_cells": "CD4+ T cells", "TCD8_cells": "CD8+ T cells",
    }

    nfkb_df = pd.DataFrame(nfkb_records)
    rows = []
    for ct in DEG_CELL_TYPES:
        pre_row = nfkb_df[(nfkb_df['cell_type'] == ct) & (nfkb_df['comparison'] == 'pre')]
        post_row = nfkb_df[(nfkb_df['cell_type'] == ct) & (nfkb_df['comparison'] == 'post')]
        pre_nes = pre_row['nes'].values[0] if len(pre_row) > 0 else np.nan
        post_nes = post_row['nes'].values[0] if len(post_row) > 0 else np.nan

        # Rank by absolute NES (descending) — compute after collecting all
        rows.append({
            'cell_type': ct,
            'label': LABELS.get(ct, ct),
            'pre_nes': pre_nes,
            'post_nes': post_nes,
        })

    rankings = pd.DataFrame(rows)

    # Compute ranks (by descending absolute NES), handle NaN gracefully
    rankings['pre_rank'] = rankings['pre_nes'].abs().rank(ascending=False, method='min', na_option='bottom')
    rankings['post_rank'] = rankings['post_nes'].abs().rank(ascending=False, method='min', na_option='bottom')
    rankings['pre_inverted'] = (50 - rankings['pre_rank']).clip(lower=0)
    rankings['post_inverted'] = (50 - rankings['post_rank']).clip(lower=0)

    # Reorder columns to match old format
    rankings = rankings[['cell_type', 'label', 'pre_rank', 'pre_nes', 'pre_inverted',
                          'post_rank', 'post_nes', 'post_inverted']]

    out_path = GSEA_DIR / "nfkb_rankings_13types_mast.csv"
    rankings.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")
    print(rankings.to_string(index=False))

    print("\n\nDone!")


if __name__ == "__main__":
    main()
