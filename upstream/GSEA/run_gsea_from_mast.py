#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
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
    # Negate logFC: MAST DEGs are R-vs-NR; flip to NR-vs-R to match upstream-pipeline convention
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
        return None


def main():

    ARCHIVED_DIR.mkdir(parents=True, exist_ok=True)
    old_files = list(GSEA_DIR.glob("*_mast_prerank_gsea.csv"))
    if old_files:
        for f in old_files:
            shutil.move(str(f), str(ARCHIVED_DIR / f.name))
    else:
        print("No old files to archive (already moved)")

    GSEA_PRE_DIR.mkdir(parents=True, exist_ok=True)
    GSEA_POST_DIR.mkdir(parents=True, exist_ok=True)

    nfkb_records = []

    for ct in DEG_CELL_TYPES:

        for comparison, deg_dir, gsea_out_dir in [
            ("pre", DEG_PRE_DIR, GSEA_PRE_DIR),
            ("post", DEG_POST_DIR, GSEA_POST_DIR),
        ]:
            deg_csv = deg_dir / f"{ct}_mast_deg.csv"
            if not deg_csv.exists():
                continue

            rnk = compute_ranking(deg_csv)

            if len(rnk) < 15:
                continue

            tmp_dir = GSEA_DIR / f"_tmp_{ct}_{comparison}"
            res2d = run_gsea_prerank(rnk, tmp_dir)
            if res2d is None:
                nfkb_records.append({'cell_type': ct, 'comparison': comparison, 'nes': np.nan})
                shutil.rmtree(tmp_dir, ignore_errors=True)
                continue

            out_csv = gsea_out_dir / f"{ct}_gsea_hallmark.csv"
            res2d.to_csv(out_csv, index=False)

            nfkb_row = res2d[res2d['Term'] == NFKB_TERM]
            if len(nfkb_row) > 0:
                nes_val = nfkb_row['NES'].values[0]
            else:
                nes_val = np.nan

            nfkb_records.append({
                'cell_type': ct,
                'comparison': comparison,
                'nes': nes_val,
            })

            shutil.rmtree(tmp_dir, ignore_errors=True)

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

        rows.append({
            'cell_type': ct,
            'label': LABELS.get(ct, ct),
            'pre_nes': pre_nes,
            'post_nes': post_nes,
        })

    rankings = pd.DataFrame(rows)

    rankings['pre_rank'] = rankings['pre_nes'].abs().rank(ascending=False, method='min', na_option='bottom')
    rankings['post_rank'] = rankings['post_nes'].abs().rank(ascending=False, method='min', na_option='bottom')
    rankings['pre_inverted'] = (50 - rankings['pre_rank']).clip(lower=0)
    rankings['post_inverted'] = (50 - rankings['post_rank']).clip(lower=0)

    rankings = rankings[['cell_type', 'label', 'pre_rank', 'pre_nes', 'pre_inverted',
                          'post_rank', 'post_nes', 'post_inverted']]

    out_path = GSEA_DIR / "nfkb_rankings_13types_mast.csv"
    rankings.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()
