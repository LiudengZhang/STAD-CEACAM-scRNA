#!/usr/bin/env python3
"""
Generate per-sample mean CEACAM5/6 expression for liver and stomach epithelial cells.
Output: liver_stomach_ceacam_per_sample.csv
Columns: site, sample, response, CEACAM5, CEACAM6
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import EPITHELIAL_H5AD
from shared.sample_ids import sample_id_map, to_study_ids

import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR = Path(__file__).parent


def main():
    print("Generating liver + stomach per-sample CEACAM5/6 expression...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)
    # This file's .X is not the output of the annotation pipeline: the
    # scale(max_value=10) clip never fired, and a third of it is NaN, in whole
    # cell rows concentrated on the deepest-sequenced cells. Read the
    # log1p CP10K matrix in .raw, which reproduces the counts deposit to 1e-2.
    if adata.raw is None:
        raise SystemExit("Epithelial.h5ad has no .raw - refusing to fall back to .X")
    adata = adata.raw.to_adata()[adata.obs_names]   # .raw.to_adata() carries obs
    # Specimens are named the way Supplementary Table 1 names them, never by the
    # identifier .obs['sample'] carries. Resolved here, from the whole object,
    # so the crosswalk is checked against every specimen in it.
    ids = sample_id_map(adata.obs)
    print(f"  Total epithelial cells: {adata.n_obs}")

    rows = []

    # --- Liver ---
    liver = adata[adata.obs['Sample site'] == 'Liver'].copy()
    liver.obs['response'] = liver.obs['liver_pre_grouping'].map({
        'Responsed': 'R', 'No-response': 'NR'
    })
    liver = liver[liver.obs['response'].isin(['R', 'NR'])].copy()
    print(f"  Liver R/NR cells: {liver.n_obs}")

    for sample in liver.obs['sample'].unique():
        mask = liver.obs['sample'] == sample
        sub = liver[mask]
        resp = sub.obs['response'].iloc[0]
        vals5 = np.asarray(sub[:, 'CEACAM5'].X.todense()).flatten()
        vals6 = np.asarray(sub[:, 'CEACAM6'].X.todense()).flatten()
        c5 = vals5.mean()
        c6 = vals6.mean()
        rows.append({'site': 'Liver', 'sample': to_study_ids([sample], ids)[0],
                      'response': resp, 'CEACAM5': c5, 'CEACAM6': c6})

    # --- Stomach ---
    stomach = adata[adata.obs['Sample site'] == 'Stomach'].copy()
    stomach.obs['response'] = stomach.obs['stomach_pre_grouping'].map({
        'Responsed': 'R', 'No-response': 'NR'
    })
    stomach = stomach[stomach.obs['response'].isin(['R', 'NR'])].copy()
    print(f"  Stomach R/NR cells: {stomach.n_obs}")

    for sample in stomach.obs['sample'].unique():
        mask = stomach.obs['sample'] == sample
        sub = stomach[mask]
        resp = sub.obs['response'].iloc[0]
        vals5 = np.asarray(sub[:, 'CEACAM5'].X.todense()).flatten()
        vals6 = np.asarray(sub[:, 'CEACAM6'].X.todense()).flatten()
        c5 = vals5.mean()
        c6 = vals6.mean()
        rows.append({'site': 'Stomach', 'sample': to_study_ids([sample], ids)[0],
                      'response': resp, 'CEACAM5': c5, 'CEACAM6': c6})

    df = pd.DataFrame(rows)
    out_path = OUTPUT_DIR / 'liver_stomach_ceacam_per_sample.csv'
    df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")
    print(f"  Liver:   {(df['site']=='Liver').sum()} samples (R={((df['site']=='Liver')&(df['response']=='R')).sum()}, NR={((df['site']=='Liver')&(df['response']=='NR')).sum()})")
    print(f"  Stomach: {(df['site']=='Stomach').sum()} samples (R={((df['site']=='Stomach')&(df['response']=='R')).sum()}, NR={((df['site']=='Stomach')&(df['response']=='NR')).sum()})")
    print(df.to_string(index=False))


if __name__ == '__main__':
    main()
