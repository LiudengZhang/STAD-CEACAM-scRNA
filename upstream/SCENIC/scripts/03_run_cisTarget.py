#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Step 3: cisTarget Motif Enrichment Analysis - Full MoMac (12 Post-Stomach Samples)
==================================================================================
Purpose: Prune co-expression networks using motif databases to identify
         direct TF-target regulatory relationships

Input:
  - adjacencies.tsv (from Step 2 on compute node)
  - SCENIC databases (hg38 motif rankings + annotations)

Output:
  - regulons.csv (validated TF-target regulons)
  - regulons.pkl (pickle format for AUCell)
  - motif_enrichment.csv (motif enrichment results)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import pickle
import loompy as lp

def main():

    base_dir = Path('/path/to/upstream-pipeline')
    project_dir = base_dir / '02_Playground/02_Cell-Cell_Interaction/01_Analysis/09_SCENIC/17_Full_MoMac_SCENIC'
    db_dir = base_dir / '02_Playground/02_Cell-Cell_Interaction/01_Analysis/09_SCENIC/00_Database'

    loom_file = project_dir / 'Results/MoMac_12sample_scenic.loom'
    adjacencies_file = project_dir / 'Results/adjacencies.tsv'

    db_10kb = db_dir / 'hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
    motif_annotations = db_dir / 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'

    output_regulons = project_dir / 'Results/regulons.csv'
    output_motifs = project_dir / 'Results/motif_enrichment.csv'
    output_regulons_pkl = project_dir / 'Results/regulons.pkl'

    if not loom_file.exists():
        sys.exit(1)

    if not adjacencies_file.exists():
        sys.exit(1)

    db_10kb_resolved = db_10kb.resolve()
    motif_annotations_resolved = motif_annotations.resolve()

    if not db_10kb_resolved.exists():
        sys.exit(1)

    if not motif_annotations_resolved.exists():
        sys.exit(1)

    with lp.connect(str(loom_file), mode='r', validate=False) as ds:
        ex_matrix = pd.DataFrame(
            data=ds[:, :],
            index=ds.ra.Gene,
            columns=ds.ca.CellID
        ).T  # Transpose so cells are rows

    adjacencies = pd.read_csv(adjacencies_file, sep='\t')

    from pyscenic.utils import modules_from_adjacencies

    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

    module_sizes = [len(module) for module in modules]

    from ctxcore.rnkdb import FeatherRankingDatabase

    dbs = [FeatherRankingDatabase(fname=str(db_10kb_resolved), name=db_10kb.stem)]

    from pyscenic.prune import prune2df, df2regulons

    df_motifs = prune2df(
        rnkdbs=dbs,
        modules=modules,
        motif_annotations_fname=str(motif_annotations_resolved),
        num_workers=8,
        client_or_address='custom_multiprocessing'
    )

    regulons = df2regulons(df_motifs)

    regulon_sizes = [len(regulon) for regulon in regulons]

    regulon_info = [(reg.name, len(reg)) for reg in regulons]
    regulon_info_sorted = sorted(regulon_info, key=lambda x: x[1], reverse=True)

    for tf_name, size in regulon_info_sorted[:20]:
        print(f"     {tf_name}: {size} genes")

    regulon_data = []
    for regulon in regulons:
        for gene in regulon.gene2weight:
            regulon_data.append({
                'TF': regulon.name,
                'target': gene,
                'importance': regulon.gene2weight[gene]
            })

    regulon_df = pd.DataFrame(regulon_data)
    regulon_df.to_csv(output_regulons, index=False)

    with open(output_regulons_pkl, 'wb') as f:
        pickle.dump(regulons, f)

    df_motifs.to_csv(output_motifs, index=False)

    key_targets = ['IL1B', 'TNF', 'IL6', 'CCL2', 'CCL3', 'CCL4', 'CXCL8', 'S100A8', 'S100A9']
    for target in key_targets:
        regulating_tfs = regulon_df[regulon_df['target'] == target]['TF'].unique()
        if len(regulating_tfs) > 0:
            print(f"     {target} regulated by: {', '.join(regulating_tfs[:5])}")

if __name__ == "__main__":
    main()
