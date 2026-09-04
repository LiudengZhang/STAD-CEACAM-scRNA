#!/usr/bin/env python3
# Paths below refer to the upstream Round_4 processing pipeline, which is
# not part of this release. This script is included as a record of how the
# input was produced; it is not called by _run_all_panels.sh.
"""
Step 1: Prepare TCGA-STAD bulk data for BayesPrism deconvolution.

- Converts TCGA FPKM matrix from Ensembl IDs to gene symbols
- Matches gene space to the existing scRNA reference (5,020 genes)
- Copies scRNA reference files from the TIGER pipeline (already prepared)
- Outputs: tcga_bulk_counts.tsv (samples × genes, GEP format)

Reuses the SAME scRNA reference as the TIGER deconvolution.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import shutil

# One change from step1_prepare_tcga_bulk.py. The original reaches the TCGA
# data through Round_4/.../98_External/Bulk/02_TCGA_STAD, which is a symlink
# to Round_5/01_Raw_Inputs/04_External/Bulk/TCGA_STAD - a directory that no
# longer exists, because 04_External was renamed 02_External. The files
# themselves are all present; only the route to them is broken. This points
# at where they actually are.
TCGA_BASE = Path("/path/to/Project_4_05232025/Round_5/01_Raw_Inputs/02_External/Bulk/TCGA_STAD")
FPKM_FILE = TCGA_BASE / "02_Raw_Data/Gene_Expression/FPKM/TCGA_STAD_FPKM_matrix.tsv"
GENE_INFO = TCGA_BASE / "02_Raw_Data/Gene_Expression/FPKM/TCGA_STAD_gene_info.tsv"
CLINICAL = TCGA_BASE / "02_Raw_Data/Clinical/TCGA_STAD_clinical_data.tsv"

# Existing scRNA reference from TIGER pipeline
TIGER_BP = Path("/path/to/Project_4_05232025/Round_5/02_Preparation_for_Panels/BayesPrism")

OUTPUT_DIR = Path(__file__).parent


def main():
    print("=" * 60)
    print("Step 1: Prepare TCGA-STAD bulk for BayesPrism")
    print("=" * 60)

    # --- Copy scRNA reference from TIGER pipeline ---
    print("\nCopying scRNA reference from TIGER pipeline...")
    for fname in ["sc_counts.tsv", "sc_cell_types.tsv"]:
        src = TIGER_BP / fname
        dst = OUTPUT_DIR / fname
        if not dst.exists():
            shutil.copy2(src, dst)
            print(f"  Copied: {fname}")
        else:
            print(f"  Already exists: {fname}")

    # Get reference gene list
    sc_genes = pd.read_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t', index_col=0, nrows=0).columns.tolist()
    print(f"scRNA reference genes: {len(sc_genes)}")

    # --- Load gene mapping ---
    print("\nLoading Ensembl → symbol mapping...")
    gene_info = pd.read_csv(GENE_INFO, sep='\t')
    gene_info = gene_info.dropna(subset=['gene_name'])
    # Remove non-gene rows (N_unmapped etc.)
    gene_info = gene_info[gene_info['gene_id'].str.startswith('ENSG')]
    # Keep protein-coding genes (primary)
    ens_to_symbol = dict(zip(gene_info['gene_id'], gene_info['gene_name']))
    print(f"Mapping entries: {len(ens_to_symbol)}")

    # --- Load TCGA FPKM ---
    print("\nLoading TCGA-STAD FPKM matrix...")
    fpkm = pd.read_csv(FPKM_FILE, sep='\t', index_col=0)
    # Remove non-gene rows
    fpkm = fpkm[fpkm.index.str.startswith('ENSG')]
    print(f"FPKM: {fpkm.shape[0]:,} genes × {fpkm.shape[1]} samples")

    # Map to symbols
    fpkm['gene_symbol'] = fpkm.index.map(ens_to_symbol)
    fpkm = fpkm.dropna(subset=['gene_symbol'])
    # Handle duplicates: keep the one with highest mean expression
    fpkm['mean_expr'] = fpkm.drop(columns=['gene_symbol']).mean(axis=1)
    fpkm = fpkm.sort_values('mean_expr', ascending=False).drop_duplicates(subset='gene_symbol', keep='first')
    fpkm = fpkm.drop(columns=['gene_symbol', 'mean_expr'])
    fpkm.index = fpkm.index.map(ens_to_symbol)
    print(f"After symbol mapping: {fpkm.shape[0]:,} unique genes")

    # --- Match to scRNA reference gene space ---
    common_genes = sorted(set(sc_genes) & set(fpkm.index))
    print(f"Genes common with scRNA reference: {len(common_genes)} / {len(sc_genes)}")

    # Key genes check
    for g in ['IL1B', 'IL6', 'TNF', 'NFKB1', 'CEACAM5', 'CEACAM6', 'CD274', 'OSM']:
        status = "found" if g in common_genes else "MISSING"
        print(f"  {g}: {status}")

    # Subset and transpose to samples × genes
    tcga_bulk = fpkm.loc[common_genes].T
    print(f"\nTCGA bulk matrix: {tcga_bulk.shape[0]} samples × {tcga_bulk.shape[1]} genes")

    # Save
    tcga_bulk.to_csv(OUTPUT_DIR / "tcga_bulk_counts.tsv", sep='\t')
    print(f"Saved: tcga_bulk_counts.tsv")

    # --- Also reduce scRNA reference to common genes only ---
    print("\nReducing scRNA reference to common genes...")
    sc_counts = pd.read_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t', index_col=0)
    sc_reduced = sc_counts[common_genes]
    sc_reduced.to_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t')
    print(f"sc_counts.tsv reduced: {sc_reduced.shape}")

    # --- Prepare clinical metadata ---
    print("\nPreparing clinical metadata...")
    clinical = pd.read_csv(CLINICAL, sep='\t')
    # Keep columns relevant for survival
    cols = ['submitter_id', 'vital_status', 'days_to_death', 'days_to_last_follow_up']
    cols = [c for c in cols if c in clinical.columns]
    clinical_out = clinical[cols].copy()
    # Match to FPKM samples
    tcga_samples = tcga_bulk.index.tolist()
    clinical_out = clinical_out[clinical_out['submitter_id'].isin(tcga_samples)]
    clinical_out.to_csv(OUTPUT_DIR / "tcga_clinical.tsv", sep='\t', index=False)
    print(f"Clinical data: {len(clinical_out)} samples (of {len(tcga_samples)} in FPKM)")

    print("\nStep 1 complete!")
    print(f"Files ready for BayesPrism in: {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
