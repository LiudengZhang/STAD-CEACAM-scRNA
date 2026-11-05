#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Step 1b: Reduce gene set from 30K to ~5K for faster BayesPrism.
Selects top variable genes from bulk + ensures CEACAM genes are included.
Overwrites sc_counts.tsv and bulk_counts.tsv with reduced gene sets.
"""

import pandas as pd
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent
N_GENES = 5000


def main():
    bulk = pd.read_csv(OUTPUT_DIR / "bulk_counts.tsv", sep='\t', index_col=0)

    gene_var = bulk.var(axis=0).sort_values(ascending=False)
    top_genes = set(gene_var.head(N_GENES).index)

    force = [g for g in bulk.columns if g.startswith('CEACAM') or g == 'CD274']
    top_genes.update(force)

    top_genes = sorted(top_genes)

    bulk_sub = bulk[top_genes]
    bulk_sub.to_csv(OUTPUT_DIR / "bulk_counts.tsv", sep='\t')
    sc = pd.read_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t', index_col=0)
    common = [g for g in top_genes if g in sc.columns]
    sc_sub = sc[common]
    sc_sub.to_csv(OUTPUT_DIR / "sc_counts.tsv", sep='\t')


if __name__ == '__main__':
    main()
