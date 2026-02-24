#!/usr/bin/env python3
"""
Step 1b: Reduce gene set from 30K to ~5K for faster BayesPrism.
Selects top variable genes from bulk + ensures CEACAM genes are included.
Overwrites sc_counts.tsv and bulk_counts.tsv with reduced gene sets.
"""

import pandas as pd
import numpy as np
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent
N_GENES = 5000


def main():
    print("Loading bulk data to select variable genes...")
    bulk = pd.read_csv(OUTPUT_DIR / "bulk_counts.tsv", sep='\t', index_col=0)
    print(f"Original bulk: {bulk.shape}")

    # Variance across samples (in bulk)
    gene_var = bulk.var(axis=0).sort_values(ascending=False)
    top_genes = set(gene_var.head(N_GENES).index)

    # Force-include CEACAM genes and CD274
    force = [g for g in bulk.columns if g.startswith('CEACAM') or g == 'CD274']
    top_genes.update(force)
    print(f"Forced genes: {force}")

    top_genes = sorted(top_genes)
    print(f"Selected genes: {len(top_genes)}")

    # Subset bulk
    bulk_sub = bulk[top_genes]
    bulk_sub.to_csv(OUTPUT_DIR / "bulk_counts.tsv", sep='\t')
    print(f"Bulk reduced: {bulk_sub.shape}")

    # Subset scRNA-seq
    print("\nLoading scRNA-seq counts...")
    sc = pd.read_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t', index_col=0)
    print(f"Original sc: {sc.shape}")
    common = [g for g in top_genes if g in sc.columns]
    sc_sub = sc[common]
    sc_sub.to_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t')
    print(f"sc reduced: {sc_sub.shape}")

    print(f"\nCEACAM genes in final set: {[g for g in common if 'CEACAM' in g]}")
    print("Done!")


if __name__ == '__main__':
    main()
