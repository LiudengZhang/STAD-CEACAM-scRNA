#!/usr/bin/env python3
"""
CellRank Analysis - MoMac Fate Mapping
======================================
Uses PseudotimeKernel with Diffusion Pseudotime (no velocity required)
"""

import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def main():
    print("=" * 80)
    print("CellRank Analysis - MoMac Fate Mapping")
    print("=" * 80)
    print()

    # Central config
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
    from paths import SCENIC_RESULTS_DIR

    input_file = SCENIC_RESULTS_DIR / 'MoMac_12sample_post_stomach.h5ad'
    aucell_file = SCENIC_RESULTS_DIR / 'aucell_matrix.csv'
    output_dir = SCENIC_RESULTS_DIR / 'CellRank'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / 'figures'
    fig_dir.mkdir(exist_ok=True)

    # Step 1: Load data
    print("[1/8] Loading MoMac data...")
    adata = sc.read_h5ad(input_file)
    print(f"   Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    print(f"   Annotations: {adata.obs['MoMac_annotation'].value_counts().to_dict()}")
    print()

    # Step 2: Load AUCell scores
    print("[2/8] Loading SCENIC AUCell scores...")
    aucell = pd.read_csv(aucell_file, index_col=0)
    common_cells = adata.obs_names.intersection(aucell.index)
    print(f"   Matched {len(common_cells):,} cells with AUCell scores")

    for col in aucell.columns:
        adata.obs[f'AUCell_{col}'] = aucell.loc[adata.obs_names, col].values
    print(f"   Added {len(aucell.columns)} regulon scores to adata.obs")
    print()

    # Step 3: Compute diffusion map and pseudotime
    print("[3/8] Computing Diffusion Pseudotime...")

    # Need to recompute neighbors on the processed data
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, use_rep='X_pca')

    # Compute diffusion map
    sc.tl.diffmap(adata, n_comps=15)
    print("   Diffusion map computed")

    # Find root cell (most undifferentiated - likely a monocyte)
    # Use gene count as proxy for stemness (more genes = less differentiated)
    adata.obs['n_genes_temp'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, 'A1') else (adata.X > 0).sum(axis=1)

    # Find cells in clusters likely to be monocytes (less differentiated)
    # Based on annotation, look for cluster patterns
    root_candidates = adata.obs['n_genes_temp'].values
    root_idx = np.argmax(root_candidates)
    adata.uns['iroot'] = root_idx
    print(f"   Root cell index: {root_idx}")
    print(f"   Root cell annotation: {adata.obs['MoMac_annotation'].iloc[root_idx]}")

    # Compute DPT
    sc.tl.dpt(adata, n_dcs=10)
    print(f"   DPT range: {adata.obs['dpt_pseudotime'].min():.3f} - {adata.obs['dpt_pseudotime'].max():.3f}")
    print()

    # Step 4: Initialize PseudotimeKernel
    print("[4/8] Initializing PseudotimeKernel...")
    pk = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
    pk.compute_transition_matrix(threshold_scheme="soft", frac_to_keep=0.3)
    print("   Transition matrix computed")
    print()

    # Step 5: Identify terminal states
    print("[5/8] Identifying terminal states...")
    g = cr.estimators.GPCCA(pk)
    g.compute_schur(n_components=15)
    g.compute_macrostates(n_states=5, cluster_key="MoMac_annotation")
    print(f"   Macrostates: {list(g.macrostates.cat.categories)}")

    g.set_terminal_states()
    terminal_states = list(g.terminal_states.cat.categories)
    print(f"   Terminal states: {terminal_states}")
    print()

    # Step 6: Compute fate probabilities
    print("[6/8] Computing fate probabilities...")
    g.compute_fate_probabilities()
    print("   Fate probabilities computed")

    fate_probs = g.fate_probabilities
    fate_df = pd.DataFrame(
        fate_probs.X,
        index=adata.obs_names,
        columns=fate_probs.names
    )
    print()

    # Step 7: Correlate TF activity with fates
    print("[7/8] Correlating TF regulon activity with fate probabilities...")

    aucell_cols = [c for c in adata.obs.columns if c.startswith('AUCell_')]
    correlations = []

    for fate in fate_probs.names:
        fate_prob = fate_df[fate].values
        for tf_col in aucell_cols:
            tf_activity = adata.obs[tf_col].values
            corr = np.corrcoef(fate_prob, tf_activity)[0, 1]
            correlations.append({
                'fate': fate,
                'TF': tf_col.replace('AUCell_', ''),
                'correlation': corr
            })

    corr_df = pd.DataFrame(correlations)

    # Top TFs per fate
    print("   Top TFs correlated with each fate:")
    for fate in fate_probs.names:
        fate_corrs = corr_df[corr_df['fate'] == fate].sort_values('correlation', ascending=False)
        top5 = fate_corrs.head(5)
        print(f"\n   {fate}:")
        for _, row in top5.iterrows():
            print(f"      {row['TF']}: r={row['correlation']:.3f}")
    print()

    # Step 8: Save results and visualizations
    print("[8/8] Saving results and visualizations...")

    # Save data
    adata.obs[['MoMac_annotation', 'dpt_pseudotime']].to_csv(output_dir / 'pseudotime_scores.csv')
    print(f"   Saved: pseudotime_scores.csv")

    fate_df.to_csv(output_dir / 'fate_probabilities.csv')
    print(f"   Saved: fate_probabilities.csv")

    corr_df.to_csv(output_dir / 'tf_fate_correlations.csv', index=False)
    print(f"   Saved: tf_fate_correlations.csv")

    # Visualizations
    # Plot 1: DPT on UMAP
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.umap(adata, color='dpt_pseudotime', ax=axes[0], show=False,
               title='Diffusion Pseudotime\n(0=root/undifferentiated)')
    sc.pl.umap(adata, color='MoMac_annotation', ax=axes[1], show=False,
               title='MoMac Annotation')
    plt.tight_layout()
    plt.savefig(fig_dir / 'pseudotime_umap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: pseudotime_umap.png")

    # Plot 2: DPT by annotation
    fig, ax = plt.subplots(figsize=(10, 6))
    order = adata.obs.groupby('MoMac_annotation')['dpt_pseudotime'].median().sort_values().index
    adata.obs['MoMac_annotation_cat'] = pd.Categorical(adata.obs['MoMac_annotation'], categories=order)
    sc.pl.violin(adata, keys='dpt_pseudotime', groupby='MoMac_annotation_cat', ax=ax, show=False,
                 rotation=45)
    ax.set_title('Pseudotime by MoMac Annotation')
    plt.tight_layout()
    plt.savefig(fig_dir / 'pseudotime_by_annotation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: pseudotime_by_annotation.png")

    # Plot 3: Terminal states
    fig, ax = plt.subplots(figsize=(8, 6))
    g.plot_macrostates(which="terminal", basis="umap", show=False)
    plt.savefig(fig_dir / 'terminal_states_umap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: terminal_states_umap.png")

    # Plot 4: Fate probabilities
    try:
        g.plot_fate_probabilities(basis="umap", same_plot=False)
        plt.savefig(fig_dir / 'fate_probabilities_umap.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   Saved: fate_probabilities_umap.png")
    except:
        print("   (fate probability plot skipped)")

    # Plot 5: TF-Fate correlation heatmap (top 10 TFs per fate)
    top_tfs = set()
    for fate in fate_probs.names:
        fate_corrs = corr_df[corr_df['fate'] == fate].nlargest(10, 'correlation')
        top_tfs.update(fate_corrs['TF'].values)
        fate_corrs_neg = corr_df[corr_df['fate'] == fate].nsmallest(5, 'correlation')
        top_tfs.update(fate_corrs_neg['TF'].values)

    pivot_df = corr_df[corr_df['TF'].isin(top_tfs)].pivot(index='TF', columns='fate', values='correlation')

    fig, ax = plt.subplots(figsize=(10, 12))
    im = ax.imshow(pivot_df.values, cmap='RdBu_r', aspect='auto', vmin=-0.5, vmax=0.5)
    ax.set_xticks(range(len(pivot_df.columns)))
    ax.set_xticklabels(pivot_df.columns, rotation=45, ha='right')
    ax.set_yticks(range(len(pivot_df.index)))
    ax.set_yticklabels(pivot_df.index)
    ax.set_title('TF Regulon Activity vs Fate Probability Correlation')
    plt.colorbar(im, ax=ax, label='Pearson r')
    plt.tight_layout()
    plt.savefig(fig_dir / 'tf_fate_correlation_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: tf_fate_correlation_heatmap.png")

    # Save annotated adata
    adata.write_h5ad(output_dir / 'MoMac_cellrank.h5ad')
    print(f"   Saved: MoMac_cellrank.h5ad")

    # Summary
    print()
    print("=" * 80)
    print("CellRank Analysis Complete!")
    print("=" * 80)
    print(f"Results: {output_dir}")
    print()
    print("Key findings:")
    print(f"  - Terminal states identified: {terminal_states}")
    print(f"  - TF-fate correlations computed for {len(aucell_cols)} regulons")
    print()

if __name__ == "__main__":
    main()
