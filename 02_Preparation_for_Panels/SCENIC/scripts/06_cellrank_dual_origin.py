#!/usr/bin/env python3
"""
CellRank Analysis with Dual Initial States
==========================================
Two independent origins converging to IL1B+ macrophages:
  - Origin 1: C5_Mac_Prolif_MKI67 (proliferating progenitors)
  - Origin 2: C4_Mono_Alternative_CD16 (CD16+ monocytes)
  - Terminal: C3_Mac_Inflam_IL1B (IL1B+ inflammatory macrophages)

Integrates with SCENIC AUCell scores to identify TF dynamics per branch.
"""

import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def main():
    print("=" * 80)
    print("CellRank Dual-Origin Analysis - MoMac Fate Mapping")
    print("=" * 80)
    print()

    # Central config
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
    from paths import SCENIC_RESULTS_DIR

    input_file = SCENIC_RESULTS_DIR / 'MoMac_12sample_post_stomach.h5ad'
    aucell_file = SCENIC_RESULTS_DIR / 'aucell_matrix.csv'
    output_dir = SCENIC_RESULTS_DIR / 'CellRank_DualOrigin'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / 'figures'
    fig_dir.mkdir(exist_ok=True)

    # Define biological states
    INITIAL_STATES = ['C5_Mac_Prolif_MKI67', 'C4_Mono_Alternative_CD16']
    TERMINAL_STATE = 'C3_Mac_Inflam_IL1B'
    STATE_COL = 'minor_cell_state'

    # Step 1: Load data
    print("[1/8] Loading MoMac data...")
    adata = sc.read_h5ad(input_file)
    print(f"   Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    print(f"   Cell states: {adata.obs[STATE_COL].value_counts().to_dict()}")
    print()

    # Step 2: Load AUCell scores
    print("[2/8] Loading SCENIC AUCell scores...")
    aucell = pd.read_csv(aucell_file, index_col=0)
    for col in aucell.columns:
        adata.obs[f'AUCell_{col}'] = aucell.loc[adata.obs_names, col].values
    print(f"   Added {len(aucell.columns)} regulon scores")
    print()

    # Step 3: Compute connectivity kernel
    print("[3/8] Computing connectivity kernel...")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, use_rep='X_pca')

    # Use ConnectivityKernel (no velocity/pseudotime needed)
    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()
    print("   Transition matrix computed")
    print()

    # Step 4: Initialize GPCCA estimator
    print("[4/8] Initializing GPCCA estimator...")
    g = cr.estimators.GPCCA(ck)
    g.compute_schur(n_components=20)
    print("   Schur decomposition complete")
    print()

    # Step 5: Compute macrostates
    print("[5/8] Computing macrostates...")
    g.compute_macrostates(n_states=7, cluster_key=STATE_COL)
    print(f"   Macrostates identified: {list(g.macrostates.cat.categories)}")
    print()

    # Step 6: Manually set initial and terminal states
    print("[6/8] Setting initial and terminal states...")
    print(f"   Initial states: {INITIAL_STATES}")
    print(f"   Terminal state: {TERMINAL_STATE}")

    # Find macrostates that correspond to our biological states
    macrostate_names = list(g.macrostates.cat.categories)

    # Set terminal states
    terminal_macrostates = [m for m in macrostate_names if TERMINAL_STATE in m]
    if not terminal_macrostates:
        # Try exact match
        terminal_macrostates = [m for m in macrostate_names if m == TERMINAL_STATE]

    if terminal_macrostates:
        g.set_terminal_states(terminal_macrostates)
        print(f"   Terminal macrostates set: {terminal_macrostates}")
    else:
        print(f"   WARNING: Could not find terminal state {TERMINAL_STATE} in macrostates")
        print(f"   Available macrostates: {macrostate_names}")
        # Fall back to automatic
        g.predict_terminal_states()
        print(f"   Using automatic terminal states: {list(g.terminal_states.cat.categories)}")
    print()

    # Step 7: Compute fate probabilities
    print("[7/8] Computing fate probabilities...")
    g.compute_fate_probabilities()
    print("   Fate probabilities computed")

    fate_probs = g.fate_probabilities
    fate_df = pd.DataFrame(
        fate_probs.X,
        index=adata.obs_names,
        columns=fate_probs.names
    )

    # Add fate probabilities to adata
    for col in fate_df.columns:
        adata.obs[f'fate_prob_{col}'] = fate_df[col].values
    print()

    # Step 8: Correlate TF activity with fate and trajectories
    print("[8/8] Analyzing TF dynamics per trajectory...")

    aucell_cols = [c for c in adata.obs.columns if c.startswith('AUCell_')]

    # Analyze each initial state's trajectory to terminal
    results = []

    for origin in INITIAL_STATES:
        print(f"\n   Analyzing trajectory: {origin} → {TERMINAL_STATE}")

        # Get cells from this origin
        origin_cells = adata.obs[adata.obs[STATE_COL] == origin].index
        terminal_cells = adata.obs[adata.obs[STATE_COL] == TERMINAL_STATE].index

        print(f"      Origin cells: {len(origin_cells)}")
        print(f"      Terminal cells: {len(terminal_cells)}")

        # For each TF, compute:
        # 1. Mean activity in origin vs terminal
        # 2. Fold change
        # 3. Statistical test

        for tf_col in aucell_cols:
            tf_name = tf_col.replace('AUCell_', '')

            origin_activity = adata.obs.loc[origin_cells, tf_col].values
            terminal_activity = adata.obs.loc[terminal_cells, tf_col].values

            origin_mean = np.mean(origin_activity)
            terminal_mean = np.mean(terminal_activity)

            # Log2 fold change (terminal vs origin)
            if origin_mean > 0:
                log2fc = np.log2((terminal_mean + 0.001) / (origin_mean + 0.001))
            else:
                log2fc = np.nan

            # Mann-Whitney U test
            stat, pval = stats.mannwhitneyu(origin_activity, terminal_activity, alternative='two-sided')

            results.append({
                'origin': origin,
                'terminal': TERMINAL_STATE,
                'TF': tf_name,
                'origin_mean': origin_mean,
                'terminal_mean': terminal_mean,
                'log2FC': log2fc,
                'pvalue': pval,
                'direction': 'up' if log2fc > 0 else 'down'
            })

    results_df = pd.DataFrame(results)

    # Multiple testing correction
    from scipy.stats import false_discovery_control
    results_df['padj'] = false_discovery_control(results_df['pvalue'].values, method='bh')

    # Print top TFs per trajectory
    print("\n   Top TFs INCREASING toward IL1B+ Mac:")
    for origin in INITIAL_STATES:
        branch_df = results_df[(results_df['origin'] == origin) & (results_df['direction'] == 'up')]
        branch_df = branch_df.sort_values('log2FC', ascending=False)
        print(f"\n      {origin} → IL1B+ Mac:")
        for _, row in branch_df.head(5).iterrows():
            sig = '***' if row['padj'] < 0.001 else '**' if row['padj'] < 0.01 else '*' if row['padj'] < 0.05 else ''
            print(f"         {row['TF']}: log2FC={row['log2FC']:.2f} {sig}")

    print("\n   Top TFs DECREASING toward IL1B+ Mac:")
    for origin in INITIAL_STATES:
        branch_df = results_df[(results_df['origin'] == origin) & (results_df['direction'] == 'down')]
        branch_df = branch_df.sort_values('log2FC', ascending=True)
        print(f"\n      {origin} → IL1B+ Mac:")
        for _, row in branch_df.head(5).iterrows():
            sig = '***' if row['padj'] < 0.001 else '**' if row['padj'] < 0.01 else '*' if row['padj'] < 0.05 else ''
            print(f"         {row['TF']}: log2FC={row['log2FC']:.2f} {sig}")

    # Save results
    print("\n" + "=" * 80)
    print("Saving results...")

    # Save TF dynamics
    results_df.to_csv(output_dir / 'tf_trajectory_dynamics.csv', index=False)
    print(f"   Saved: tf_trajectory_dynamics.csv")

    # Save fate probabilities
    fate_df.to_csv(output_dir / 'fate_probabilities.csv')
    print(f"   Saved: fate_probabilities.csv")

    # Visualizations
    print("\nGenerating visualizations...")

    # Plot 1: UMAP with cell states
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    sc.pl.umap(adata, color=STATE_COL, ax=axes[0], show=False,
               title='Cell States (minor_cell_state)')

    # Highlight origins and terminal
    adata.obs['trajectory_role'] = 'Other'
    adata.obs.loc[adata.obs[STATE_COL].isin(INITIAL_STATES), 'trajectory_role'] = 'Origin'
    adata.obs.loc[adata.obs[STATE_COL] == TERMINAL_STATE, 'trajectory_role'] = 'Terminal (IL1B+)'

    colors = {'Origin': '#E69F00', 'Terminal (IL1B+)': '#D55E00', 'Other': '#CCCCCC'}
    sc.pl.umap(adata, color='trajectory_role', ax=axes[1], show=False,
               title='Trajectory: 2 Origins → IL1B+ Mac',
               palette=colors)

    # Fate probability to terminal state
    if fate_df.columns[0] in adata.obs.columns or f'fate_prob_{fate_df.columns[0]}' in adata.obs.columns:
        fate_col = f'fate_prob_{fate_df.columns[0]}' if f'fate_prob_{fate_df.columns[0]}' in adata.obs.columns else fate_df.columns[0]
        sc.pl.umap(adata, color=fate_col, ax=axes[2], show=False,
                   title=f'Fate Probability → {fate_df.columns[0]}', cmap='Reds')

    plt.tight_layout()
    plt.savefig(fig_dir / 'trajectory_overview.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: trajectory_overview.png")

    # Plot 2: TF dynamics heatmap (comparing two trajectories)
    fig, ax = plt.subplots(figsize=(12, 10))

    # Pivot for heatmap: TF x Origin, values = log2FC
    pivot_df = results_df.pivot(index='TF', columns='origin', values='log2FC')

    # Select top variable TFs
    pivot_df['variance'] = pivot_df.var(axis=1)
    top_tfs = pivot_df.nlargest(30, 'variance').drop('variance', axis=1)

    sns.heatmap(top_tfs, cmap='RdBu_r', center=0, ax=ax,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'log2FC (Terminal/Origin)'})
    ax.set_title('TF Dynamics: Two Trajectories to IL1B+ Mac\n(log2 Fold Change)')
    ax.set_xlabel('Origin Population')
    ax.set_ylabel('Transcription Factor Regulon')
    plt.tight_layout()
    plt.savefig(fig_dir / 'tf_dynamics_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: tf_dynamics_heatmap.png")

    # Plot 3: Shared vs unique TFs
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Significantly changing TFs (padj < 0.05, |log2FC| > 0.5)
    sig_up = {}
    sig_down = {}
    for origin in INITIAL_STATES:
        branch_df = results_df[results_df['origin'] == origin]
        sig_up[origin] = set(branch_df[(branch_df['padj'] < 0.05) & (branch_df['log2FC'] > 0.5)]['TF'])
        sig_down[origin] = set(branch_df[(branch_df['padj'] < 0.05) & (branch_df['log2FC'] < -0.5)]['TF'])

    # Shared upregulated
    shared_up = sig_up[INITIAL_STATES[0]] & sig_up[INITIAL_STATES[1]]
    unique_up_0 = sig_up[INITIAL_STATES[0]] - sig_up[INITIAL_STATES[1]]
    unique_up_1 = sig_up[INITIAL_STATES[1]] - sig_up[INITIAL_STATES[0]]

    # Bar plot for upregulated
    categories = ['Shared', f'Only {INITIAL_STATES[0][:15]}', f'Only {INITIAL_STATES[1][:15]}']
    counts_up = [len(shared_up), len(unique_up_0), len(unique_up_1)]
    axes[0].bar(categories, counts_up, color=['#4DAF4A', '#E69F00', '#377EB8'])
    axes[0].set_title('TFs Upregulated in Both vs One Trajectory')
    axes[0].set_ylabel('Number of TFs')
    axes[0].tick_params(axis='x', rotation=15)

    # Shared downregulated
    shared_down = sig_down[INITIAL_STATES[0]] & sig_down[INITIAL_STATES[1]]
    unique_down_0 = sig_down[INITIAL_STATES[0]] - sig_down[INITIAL_STATES[1]]
    unique_down_1 = sig_down[INITIAL_STATES[1]] - sig_down[INITIAL_STATES[0]]

    counts_down = [len(shared_down), len(unique_down_0), len(unique_down_1)]
    axes[1].bar(categories, counts_down, color=['#4DAF4A', '#E69F00', '#377EB8'])
    axes[1].set_title('TFs Downregulated in Both vs One Trajectory')
    axes[1].set_ylabel('Number of TFs')
    axes[1].tick_params(axis='x', rotation=15)

    plt.tight_layout()
    plt.savefig(fig_dir / 'shared_vs_unique_tfs.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: shared_vs_unique_tfs.png")

    # Plot 4: Key inflammatory TF comparison
    key_tfs = ['FOS(+)', 'FOSB(+)', 'JUN(+)', 'JUNB(+)', 'ATF3(+)',
               'NFKB1(+)', 'NFKB2(+)', 'ETS2(+)', 'MAFB(+)', 'MAF(+)']

    key_results = results_df[results_df['TF'].isin(key_tfs)]

    if len(key_results) > 0:
        fig, ax = plt.subplots(figsize=(10, 6))

        pivot_key = key_results.pivot(index='TF', columns='origin', values='log2FC')

        x = np.arange(len(pivot_key))
        width = 0.35

        bars1 = ax.bar(x - width/2, pivot_key[INITIAL_STATES[0]], width,
                       label=f'{INITIAL_STATES[0][:20]}...', color='#E69F00')
        bars2 = ax.bar(x + width/2, pivot_key[INITIAL_STATES[1]], width,
                       label=f'{INITIAL_STATES[1][:20]}...', color='#377EB8')

        ax.axhline(0, color='black', linewidth=0.5)
        ax.set_xlabel('Transcription Factor')
        ax.set_ylabel('log2 Fold Change (Terminal/Origin)')
        ax.set_title('Key Inflammatory TF Changes Along Two Trajectories')
        ax.set_xticks(x)
        ax.set_xticklabels(pivot_key.index, rotation=45, ha='right')
        ax.legend()

        plt.tight_layout()
        plt.savefig(fig_dir / 'key_inflammatory_tfs.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   Saved: key_inflammatory_tfs.png")

    # Save adata
    adata.write_h5ad(output_dir / 'MoMac_cellrank_dual_origin.h5ad')
    print(f"   Saved: MoMac_cellrank_dual_origin.h5ad")

    # Summary
    print()
    print("=" * 80)
    print("CellRank Dual-Origin Analysis Complete!")
    print("=" * 80)
    print(f"\nResults: {output_dir}")
    print()
    print("Key findings:")
    print(f"  - Two independent origins identified:")
    print(f"      1. {INITIAL_STATES[0]} ({(adata.obs[STATE_COL] == INITIAL_STATES[0]).sum()} cells)")
    print(f"      2. {INITIAL_STATES[1]} ({(adata.obs[STATE_COL] == INITIAL_STATES[1]).sum()} cells)")
    print(f"  - Terminal state: {TERMINAL_STATE} ({(adata.obs[STATE_COL] == TERMINAL_STATE).sum()} cells)")
    print(f"  - TF dynamics analyzed for {len(aucell_cols)} regulons")
    if len(shared_up) > 0:
        print(f"  - Shared upregulated TFs: {', '.join(list(shared_up)[:5])}...")
    print()

if __name__ == "__main__":
    main()
