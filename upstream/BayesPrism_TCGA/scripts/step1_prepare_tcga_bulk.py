#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Step 1: Prepare TCGA-STAD bulk data for BayesPrism deconvolution.

- Converts TCGA FPKM matrix from Ensembl IDs to gene symbols
- Matches gene space to the existing scRNA reference (5,020 genes)
- Copies scRNA reference files from the TIGER pipeline (already prepared)
- Outputs: tcga_bulk_counts.tsv (samples × genes, GEP format)
"""
import pandas as pd
from pathlib import Path
import shutil

TCGA_BASE = Path("/path/to/upstream-pipeline/04_Final_Panels/00_Set_Ups/00_Data/98_External/Bulk/02_TCGA_STAD")
FPKM_FILE = TCGA_BASE / "02_Raw_Data/Gene_Expression/FPKM/TCGA_STAD_FPKM_matrix.tsv"
GENE_INFO = TCGA_BASE / "02_Raw_Data/Gene_Expression/FPKM/TCGA_STAD_gene_info.tsv"
CLINICAL = TCGA_BASE / "02_Raw_Data/Clinical/TCGA_STAD_clinical_data.tsv"

TIGER_BP = Path("/path/to/submission-tree/02_Preparation_for_Panels/BayesPrism")

OUTPUT_DIR = Path(__file__).parent


def main():
    for fname in ["sc_counts.tsv", "sc_cell_types.tsv"]:
        src = TIGER_BP / fname
        dst = OUTPUT_DIR / fname
        if not dst.exists():
            shutil.copy2(src, dst)
        else:
            print(f"  Already exists: {fname}")

    sc_genes = pd.read_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t', index_col=0, nrows=0).columns.tolist()
    gene_info = pd.read_csv(GENE_INFO, sep='\t')
    gene_info = gene_info.dropna(subset=['gene_name'])
    gene_info = gene_info[gene_info['gene_id'].str.startswith('ENSG')]
    ens_to_symbol = dict(zip(gene_info['gene_id'], gene_info['gene_name']))
    fpkm = pd.read_csv(FPKM_FILE, sep='\t', index_col=0)
    fpkm = fpkm[fpkm.index.str.startswith('ENSG')]

    fpkm['gene_symbol'] = fpkm.index.map(ens_to_symbol)
    fpkm = fpkm.dropna(subset=['gene_symbol'])
    fpkm['mean_expr'] = fpkm.drop(columns=['gene_symbol']).mean(axis=1)
    fpkm = fpkm.sort_values('mean_expr', ascending=False).drop_duplicates(subset='gene_symbol', keep='first')
    fpkm = fpkm.drop(columns=['gene_symbol', 'mean_expr'])
    fpkm.index = fpkm.index.map(ens_to_symbol)

    common_genes = sorted(set(sc_genes) & set(fpkm.index))

    for g in ['IL1B', 'IL6', 'TNF', 'NFKB1', 'CEACAM5', 'CEACAM6', 'CD274', 'OSM']:
        status = "found" if g in common_genes else "MISSING"

    tcga_bulk = fpkm.loc[common_genes].T

    tcga_bulk.to_csv(OUTPUT_DIR / "tcga_bulk_counts.tsv", sep='\t')
    sc_counts = pd.read_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t', index_col=0)
    sc_reduced = sc_counts[common_genes]
    sc_reduced.to_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t')
    clinical = pd.read_csv(CLINICAL, sep='\t')
    cols = ['submitter_id', 'vital_status', 'days_to_death', 'days_to_last_follow_up']
    cols = [c for c in cols if c in clinical.columns]
    clinical_out = clinical[cols].copy()
    # Match to FPKM samples
    tcga_samples = tcga_bulk.index.tolist()
    clinical_out = clinical_out[clinical_out['submitter_id'].isin(tcga_samples)]
    clinical_out.to_csv(OUTPUT_DIR / "tcga_clinical.tsv", sep='\t', index=False)


if __name__ == '__main__':
    main()
