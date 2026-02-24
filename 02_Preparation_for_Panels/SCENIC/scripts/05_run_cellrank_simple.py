#!/usr/bin/env python3
"""
CellRank Simple Analysis - MoMac TF Dynamics
=============================================
Uses Diffusion Pseudotime + TF activity correlation (no heavy GPCCA)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def main():
    print("=" * 80)
    print("CellRank Simple Analysis - MoMac TF Dynamics")
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
    print("[1/6] Loading MoMac data...")
    adata = sc.read_h5ad(input_file)
    print(f"   Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    print()

    # Step 2: Load AUCell scores
    print("[2/6] Loading SCENIC AUCell scores...")
    aucell = pd.read_csv(aucell_file, index_col=0)
    for col in aucell.columns:
        adata.obs[f'AUCell_{col}'] = aucell.loc[adata.obs_names, col].values
    print(f"   Added {len(aucell.columns)} regulon scores")
    print()

    # Step 3: Compute Diffusion Pseudotime
    print("[3/6] Computing Diffusion Pseudotime...")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, use_rep='X_pca')
    sc.tl.diffmap(adata, n_comps=15)

    # Find root: cell with highest gene count (proxy for undifferentiated)
    n_genes = (adata.X > 0).sum(axis=1)
    if hasattr(n_genes, 'A1'):
        n_genes = n_genes.A1
    adata.obs['n_genes_expr'] = n_genes
    root_idx = np.argmax(n_genes)
    adata.uns['iroot'] = root_idx

    sc.tl.dpt(adata, n_dcs=10)
    print(f"   Root annotation: {adata.obs['MoMac_annotation'].iloc[root_idx]}")
    print(f"   Pseudotime range: {adata.obs['dpt_pseudotime'].min():.3f} - {adata.obs['dpt_pseudotime'].max():.3f}")
    print()

    # Step 4: Correlate TF activity with pseudotime
    print("[4/6] Correlating TF regulon activity with pseudotime...")
    pseudotime = adata.obs['dpt_pseudotime'].values
    aucell_cols = [c for c in adata.obs.columns if c.startswith('AUCell_')]

    correlations = []
    for tf_col in aucell_cols:
        tf_activity = adata.obs[tf_col].values
        r, p = stats.spearmanr(pseudotime, tf_activity)
        correlations.append({
            'TF': tf_col.replace('AUCell_', ''),
            'spearman_r': r,
            'p_value': p,
            'direction': 'increases' if r > 0 else 'decreases'
        })

    corr_df = pd.DataFrame(correlations)
    corr_df = corr_df.sort_values('spearman_r', ascending=False)

    # Top increasing TFs (late/differentiated)
    print("\n   TFs INCREASING with pseudotime (late/differentiated):")
    for _, row in corr_df.head(10).iterrows():
        sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else ''
        print(f"      {row['TF']}: r={row['spearman_r']:.3f} {sig}")

    # Top decreasing TFs (early/undifferentiated)
    print("\n   TFs DECREASING with pseudotime (early/undifferentiated):")
    for _, row in corr_df.tail(10).iloc[::-1].iterrows():
        sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else ''
        print(f"      {row['TF']}: r={row['spearman_r']:.3f} {sig}")
    print()

    # Step 5: Analyze by MoMac annotation
    print("[5/6] Analyzing pseudotime by MoMac subtype...")
    pt_by_annot = adata.obs.groupby('MoMac_annotation')['dpt_pseudotime'].agg(['mean', 'median', 'std'])
    pt_by_annot = pt_by_annot.sort_values('median')
    print("\n   Pseudotime by annotation (sorted):")
    for annot, row in pt_by_annot.iterrows():
        n_cells = (adata.obs['MoMac_annotation'] == annot).sum()
        print(f"      {annot}: median={row['median']:.3f}, n={n_cells}")
    print()

    # Step 6: Save results and visualizations
    print("[6/6] Saving results and visualizations...")

    # Save correlations
    corr_df.to_csv(output_dir / 'tf_pseudotime_correlations.csv', index=False)
    print(f"   Saved: tf_pseudotime_correlations.csv")

    # Save pseudotime data
    adata.obs[['MoMac_annotation', 'dpt_pseudotime', 'n_genes_expr']].to_csv(
        output_dir / 'pseudotime_scores.csv')
    print(f"   Saved: pseudotime_scores.csv")

    # Plot 1: Pseudotime UMAP
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.umap(adata, color='dpt_pseudotime', ax=axes[0], show=False,
               title='Diffusion Pseudotime\n(0=undifferentiated)', cmap='viridis')
    sc.pl.umap(adata, color='MoMac_annotation', ax=axes[1], show=False,
               title='MoMac Annotation')
    plt.tight_layout()
    plt.savefig(fig_dir / 'pseudotime_umap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: pseudotime_umap.png")

    # Plot 2: Pseudotime by annotation violin
    fig, ax = plt.subplots(figsize=(12, 6))
    order = pt_by_annot.index.tolist()
    sc.pl.violin(adata, keys='dpt_pseudotime', groupby='MoMac_annotation',
                 order=order, ax=ax, show=False, rotation=45)
    ax.set_title('Pseudotime by MoMac Annotation (ordered by median)')
    ax.set_ylabel('Diffusion Pseudotime')
    plt.tight_layout()
    plt.savefig(fig_dir / 'pseudotime_by_annotation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: pseudotime_by_annotation.png")

    # Plot 3: Top TF correlations bar plot
    fig, ax = plt.subplots(figsize=(10, 8))
    top_tfs = pd.concat([corr_df.head(15), corr_df.tail(15)])
    colors = ['#d62728' if r > 0 else '#1f77b4' for r in top_tfs['spearman_r']]
    ax.barh(range(len(top_tfs)), top_tfs['spearman_r'].values, color=colors)
    ax.set_yticks(range(len(top_tfs)))
    ax.set_yticklabels(top_tfs['TF'].values)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('Spearman correlation with pseudotime')
    ax.set_title('TF Regulon Activity vs Pseudotime\n(red=late/differentiated, blue=early/undifferentiated)')
    plt.tight_layout()
    plt.savefig(fig_dir / 'tf_pseudotime_correlation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: tf_pseudotime_correlation.png")

    # Plot 4: Key inflammatory TF dynamics
    key_tfs = ['AUCell_FOSB(+)', 'AUCell_FOS(+)', 'AUCell_JUN(+)', 'AUCell_ATF3(+)',
               'AUCell_NFKB1(+)', 'AUCell_NFKB2(+)', 'AUCell_ETS2(+)', 'AUCell_MAFB(+)']
    available_tfs = [tf for tf in key_tfs if tf in adata.obs.columns]

    if available_tfs:
        fig, axes = plt.subplots(2, 4, figsize=(16, 8))
        axes = axes.flatten()
        for i, tf in enumerate(available_tfs[:8]):
            ax = axes[i]
            x = adata.obs['dpt_pseudotime'].values
            y = adata.obs[tf].values
            ax.scatter(x, y, alpha=0.1, s=1)
            # Add trend line
            z = np.polyfit(x, y, 2)
            p = np.poly1d(z)
            x_line = np.linspace(0, 1, 100)
            ax.plot(x_line, p(x_line), 'r-', linewidth=2)
            ax.set_xlabel('Pseudotime')
            ax.set_ylabel('AUCell Score')
            ax.set_title(tf.replace('AUCell_', ''))
        plt.tight_layout()
        plt.savefig(fig_dir / 'key_tf_dynamics.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   Saved: key_tf_dynamics.png")

    # Save adata
    adata.write_h5ad(output_dir / 'MoMac_cellrank.h5ad')
    print(f"   Saved: MoMac_cellrank.h5ad")

    # Summary
    print()
    print("=" * 80)
    print("Analysis Complete!")
    print("=" * 80)
    print()
    print("Key findings:")
    print(f"  - Early TFs (undifferentiated): {', '.join(corr_df.tail(5)['TF'].values)}")
    print(f"  - Late TFs (differentiated): {', '.join(corr_df.head(5)['TF'].values)}")
    print()
    print(f"Results: {output_dir}")

if __name__ == "__main__":
    main()
