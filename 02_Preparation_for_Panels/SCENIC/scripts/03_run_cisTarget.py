#!/usr/bin/env python3
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

Runtime: ~1-2 hours
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import pickle
import loompy as lp

def main():
    print("=" * 80)
    print("SCENIC Step 3: cisTarget Motif Enrichment - Full MoMac (12 samples)")
    print("=" * 80)
    print()

    # Central config
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
    from paths import SCENIC_RESULTS_DIR, SCENIC_DB_DIR

    # Input files
    loom_file = SCENIC_RESULTS_DIR / 'MoMac_12sample_scenic.loom'
    adjacencies_file = SCENIC_RESULTS_DIR / 'adjacencies.tsv'

    # Database files (download from https://resources.aertslab.org/cistarget/)
    db_10kb = SCENIC_DB_DIR / 'hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
    motif_annotations = SCENIC_DB_DIR / 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'

    # Output files
    output_regulons = SCENIC_RESULTS_DIR / 'regulons.csv'
    output_motifs = SCENIC_RESULTS_DIR / 'motif_enrichment.csv'
    output_regulons_pkl = SCENIC_RESULTS_DIR / 'regulons.pkl'

    # Check input files
    if not loom_file.exists():
        print(f"ERROR: Loom file not found: {loom_file}")
        sys.exit(1)

    if not adjacencies_file.exists():
        print(f"ERROR: Adjacencies file not found: {adjacencies_file}")
        sys.exit(1)

    # Check database files (resolve symlinks)
    db_10kb_resolved = db_10kb.resolve()
    motif_annotations_resolved = motif_annotations.resolve()

    if not db_10kb_resolved.exists():
        print(f"ERROR: Database file not found: {db_10kb}")
        sys.exit(1)

    if not motif_annotations_resolved.exists():
        print(f"ERROR: Motif annotations not found: {motif_annotations}")
        sys.exit(1)

    # Step 1: Load expression matrix and adjacencies
    print("[1/5] Loading expression matrix and adjacency matrix...")
    print(f"   Loom file: {loom_file}")
    print(f"   Adjacencies file: {adjacencies_file}")

    # Load expression matrix from loom
    with lp.connect(str(loom_file), mode='r', validate=False) as ds:
        ex_matrix = pd.DataFrame(
            data=ds[:, :],
            index=ds.ra.Gene,
            columns=ds.ca.CellID
        ).T  # Transpose so cells are rows
        print(f"   Loaded expression matrix: {ex_matrix.shape[0]:,} cells x {ex_matrix.shape[1]:,} genes")

    # Load adjacencies
    adjacencies = pd.read_csv(adjacencies_file, sep='\t')
    print(f"   Loaded {len(adjacencies):,} TF-gene interactions")
    print(f"   Unique TFs: {adjacencies['TF'].nunique()}")
    print(f"   Unique targets: {adjacencies['target'].nunique()}")
    print()

    # Step 2: Create modules from adjacencies
    print("[2/5] Creating co-expression modules...")

    from pyscenic.utils import modules_from_adjacencies

    # Convert adjacencies to modules
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

    print(f"   Created {len(modules)} modules")
    print()

    # Print module statistics
    module_sizes = [len(module) for module in modules]
    print("   Module size statistics:")
    print(f"     Mean: {np.mean(module_sizes):.1f} genes")
    print(f"     Median: {np.median(module_sizes):.1f} genes")
    print(f"     Min: {np.min(module_sizes)} genes")
    print(f"     Max: {np.max(module_sizes)} genes")
    print()

    # Step 3: Load ranking databases
    print("[3/5] Loading SCENIC ranking databases...")
    print(f"   Database: {db_10kb.name}")

    from ctxcore.rnkdb import FeatherRankingDatabase

    dbs = [FeatherRankingDatabase(fname=str(db_10kb_resolved), name=db_10kb.stem)]

    print(f"   Loaded {len(dbs)} ranking database(s)")
    print()

    # Step 4: Run cisTarget motif enrichment
    print("[4/5] Running cisTarget motif enrichment...")
    print("   This step prunes co-expression modules using motif information")
    print("   to identify direct TF-target relationships")
    print("   Runtime: ~1-2 hours")
    print()

    from pyscenic.prune import prune2df, df2regulons

    # Prune modules using cisTarget
    df_motifs = prune2df(
        rnkdbs=dbs,
        modules=modules,
        motif_annotations_fname=str(motif_annotations_resolved),
        num_workers=8,
        client_or_address='custom_multiprocessing'
    )

    print("   cisTarget enrichment complete!")
    print(f"   Enriched {len(df_motifs)} TF-motif-target associations")
    print()

    # Step 5: Convert to regulons and save
    print("[5/5] Creating regulons and saving results...")

    # Convert to regulons
    regulons = df2regulons(df_motifs)

    print(f"   Created {len(regulons)} regulons")
    print()

    # Print regulon statistics
    regulon_sizes = [len(regulon) for regulon in regulons]
    print("   Regulon size statistics:")
    print(f"     Mean: {np.mean(regulon_sizes):.1f} genes")
    print(f"     Median: {np.median(regulon_sizes):.1f} genes")
    print(f"     Min: {np.min(regulon_sizes)} genes")
    print(f"     Max: {np.max(regulon_sizes)} genes")
    print()

    # Top regulons by size
    regulon_info = [(reg.name, len(reg)) for reg in regulons]
    regulon_info_sorted = sorted(regulon_info, key=lambda x: x[1], reverse=True)

    print("   Top 20 regulons by size:")
    for tf_name, size in regulon_info_sorted[:20]:
        print(f"     {tf_name}: {size} genes")
    print()

    # Save regulons as CSV
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
    print(f"   Saved regulons: {output_regulons}")
    print(f"   File size: {output_regulons.stat().st_size / 1e6:.2f} MB")

    # Save regulons as pickle for later use
    with open(output_regulons_pkl, 'wb') as f:
        pickle.dump(regulons, f)
    print(f"   Saved regulons (pickle): {output_regulons_pkl}")

    # Save motif enrichment results
    df_motifs.to_csv(output_motifs, index=False)
    print(f"   Saved motif enrichment: {output_motifs}")
    print()

    # Check for key TFs regulating inflammatory genes
    key_targets = ['IL1B', 'TNF', 'IL6', 'CCL2', 'CCL3', 'CCL4', 'CXCL8', 'S100A8', 'S100A9']
    print("   Checking regulons for key inflammatory targets:")
    for target in key_targets:
        regulating_tfs = regulon_df[regulon_df['target'] == target]['TF'].unique()
        if len(regulating_tfs) > 0:
            print(f"     {target} regulated by: {', '.join(regulating_tfs[:5])}")
    print()

    # Summary
    print("=" * 80)
    print("cisTarget Motif Enrichment Complete!")
    print("=" * 80)
    print(f"Results:")
    print(f"  - Regulons: {len(regulons)}")
    print(f"  - Total TF-target pairs: {len(regulon_df):,}")
    print(f"  - Output files:")
    print(f"      {output_regulons}")
    print(f"      {output_regulons_pkl}")
    print(f"      {output_motifs}")
    print()
    print("Next step: Run AUCell scoring with 04_run_aucell.py")
    print("=" * 80)

if __name__ == "__main__":
    main()
