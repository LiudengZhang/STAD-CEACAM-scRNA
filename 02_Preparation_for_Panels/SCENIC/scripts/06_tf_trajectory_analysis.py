#!/usr/bin/env python3
"""
TF Trajectory Analysis - Dual Origin to IL1B+ Mac
==================================================
Direct comparison of TF regulon activity between:
  - Origin 1: C5_Mac_Prolif_MKI67 (proliferating progenitors)
  - Origin 2: C4_Mono_Alternative_CD16 (CD16+ monocytes)
  - Terminal: C3_Mac_Inflam_IL1B (IL1B+ inflammatory macrophages)

Uses SCENIC AUCell scores to identify TF dynamics per trajectory branch.
No heavy GPCCA computation needed - direct population comparison.
"""

import scanpy as sc
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
    print("TF Trajectory Analysis - Dual Origin to IL1B+ Mac")
    print("=" * 80)
    print()

    # Central config
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
    from paths import SCENIC_RESULTS_DIR

    input_file = SCENIC_RESULTS_DIR / 'MoMac_12sample_post_stomach.h5ad'
    aucell_file = SCENIC_RESULTS_DIR / 'aucell_matrix.csv'
    output_dir = SCENIC_RESULTS_DIR / 'TF_Trajectory_Analysis'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / 'figures'
    fig_dir.mkdir(exist_ok=True)

    # Define biological states
    INITIAL_STATES = ['C5_Mac_Prolif_MKI67', 'C4_Mono_Alternative_CD16']
    TERMINAL_STATE = 'C3_Mac_Inflam_IL1B'
    STATE_COL = 'minor_cell_state'

    # Step 1: Load data
    print("[1/5] Loading MoMac data...")
    adata = sc.read_h5ad(input_file)
    print(f"   Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    print()

    # Step 2: Load AUCell scores
    print("[2/5] Loading SCENIC AUCell scores...")
    aucell = pd.read_csv(aucell_file, index_col=0)
    for col in aucell.columns:
        adata.obs[f'AUCell_{col}'] = aucell.loc[adata.obs_names, col].values
    print(f"   Added {len(aucell.columns)} regulon scores")

    # Print cell counts
    print(f"\n   Population sizes:")
    for state in INITIAL_STATES + [TERMINAL_STATE]:
        n = (adata.obs[STATE_COL] == state).sum()
        print(f"      {state}: {n} cells")
    print()

    # Step 3: Compute TF dynamics for each trajectory
    print("[3/5] Computing TF dynamics per trajectory...")

    aucell_cols = [c for c in adata.obs.columns if c.startswith('AUCell_')]
    results = []

    for origin in INITIAL_STATES:
        print(f"\n   Trajectory: {origin} → {TERMINAL_STATE}")

        origin_cells = adata.obs[adata.obs[STATE_COL] == origin].index
        terminal_cells = adata.obs[adata.obs[STATE_COL] == TERMINAL_STATE].index

        for tf_col in aucell_cols:
            tf_name = tf_col.replace('AUCell_', '')

            origin_activity = adata.obs.loc[origin_cells, tf_col].values
            terminal_activity = adata.obs.loc[terminal_cells, tf_col].values

            origin_mean = np.mean(origin_activity)
            terminal_mean = np.mean(terminal_activity)
            origin_std = np.std(origin_activity)
            terminal_std = np.std(terminal_activity)

            # Log2 fold change
            log2fc = np.log2((terminal_mean + 0.001) / (origin_mean + 0.001))

            # Mann-Whitney U test
            stat, pval = stats.mannwhitneyu(origin_activity, terminal_activity, alternative='two-sided')

            results.append({
                'origin': origin,
                'terminal': TERMINAL_STATE,
                'TF': tf_name,
                'origin_mean': origin_mean,
                'origin_std': origin_std,
                'terminal_mean': terminal_mean,
                'terminal_std': terminal_std,
                'log2FC': log2fc,
                'pvalue': pval
            })

    results_df = pd.DataFrame(results)

    # FDR correction
    results_df['padj'] = stats.false_discovery_control(results_df['pvalue'].values, method='bh')
    results_df['significant'] = results_df['padj'] < 0.05
    results_df['direction'] = results_df['log2FC'].apply(lambda x: 'up' if x > 0 else 'down')

    # Print top results
    print("\n   Top TFs INCREASING (upregulated in IL1B+ Mac):")
    for origin in INITIAL_STATES:
        branch_df = results_df[(results_df['origin'] == origin) & (results_df['log2FC'] > 0)]
        branch_df = branch_df.sort_values('log2FC', ascending=False)
        print(f"\n      From {origin[:25]}:")
        for _, row in branch_df.head(5).iterrows():
            sig = '***' if row['padj'] < 0.001 else '**' if row['padj'] < 0.01 else '*' if row['padj'] < 0.05 else ''
            print(f"         {row['TF']}: log2FC={row['log2FC']:.2f} {sig}")

    print("\n   Top TFs DECREASING (downregulated in IL1B+ Mac):")
    for origin in INITIAL_STATES:
        branch_df = results_df[(results_df['origin'] == origin) & (results_df['log2FC'] < 0)]
        branch_df = branch_df.sort_values('log2FC', ascending=True)
        print(f"\n      From {origin[:25]}:")
        for _, row in branch_df.head(5).iterrows():
            sig = '***' if row['padj'] < 0.001 else '**' if row['padj'] < 0.01 else '*' if row['padj'] < 0.05 else ''
            print(f"         {row['TF']}: log2FC={row['log2FC']:.2f} {sig}")

    # Step 4: Identify shared vs unique TF changes
    print("\n[4/5] Identifying shared vs trajectory-specific TF changes...")

    sig_up = {}
    sig_down = {}
    for origin in INITIAL_STATES:
        branch_df = results_df[results_df['origin'] == origin]
        sig_up[origin] = set(branch_df[(branch_df['padj'] < 0.05) & (branch_df['log2FC'] > 0.5)]['TF'])
        sig_down[origin] = set(branch_df[(branch_df['padj'] < 0.05) & (branch_df['log2FC'] < -0.5)]['TF'])

    # Shared
    shared_up = sig_up[INITIAL_STATES[0]] & sig_up[INITIAL_STATES[1]]
    shared_down = sig_down[INITIAL_STATES[0]] & sig_down[INITIAL_STATES[1]]

    # Unique
    unique_up_prolif = sig_up[INITIAL_STATES[0]] - sig_up[INITIAL_STATES[1]]
    unique_up_mono = sig_up[INITIAL_STATES[1]] - sig_up[INITIAL_STATES[0]]
    unique_down_prolif = sig_down[INITIAL_STATES[0]] - sig_down[INITIAL_STATES[1]]
    unique_down_mono = sig_down[INITIAL_STATES[1]] - sig_down[INITIAL_STATES[0]]

    print(f"\n   SHARED TFs (change in BOTH trajectories):")
    print(f"      Upregulated: {len(shared_up)} TFs")
    if shared_up:
        print(f"         {', '.join(sorted(shared_up))}")
    print(f"      Downregulated: {len(shared_down)} TFs")
    if shared_down:
        print(f"         {', '.join(sorted(shared_down))}")

    print(f"\n   TRAJECTORY-SPECIFIC TFs:")
    print(f"      Only in Prolif→IL1B+: {len(unique_up_prolif)} up, {len(unique_down_prolif)} down")
    print(f"      Only in CD16 Mono→IL1B+: {len(unique_up_mono)} up, {len(unique_down_mono)} down")

    # Step 5: Generate visualizations
    print("\n[5/5] Generating visualizations...")

    # Plot 1: UMAP overview
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    sc.pl.umap(adata, color=STATE_COL, ax=axes[0], show=False, title='Cell States')

    # Highlight trajectory
    adata.obs['trajectory_role'] = 'Other'
    adata.obs.loc[adata.obs[STATE_COL] == INITIAL_STATES[0], 'trajectory_role'] = 'Origin: Prolif'
    adata.obs.loc[adata.obs[STATE_COL] == INITIAL_STATES[1], 'trajectory_role'] = 'Origin: CD16 Mono'
    adata.obs.loc[adata.obs[STATE_COL] == TERMINAL_STATE, 'trajectory_role'] = 'Terminal: IL1B+ Mac'

    colors = {'Origin: Prolif': '#E69F00', 'Origin: CD16 Mono': '#377EB8',
              'Terminal: IL1B+ Mac': '#D55E00', 'Other': '#CCCCCC'}
    sc.pl.umap(adata, color='trajectory_role', ax=axes[1], show=False,
               title='Two Trajectories → IL1B+ Mac', palette=colors)

    # Key TF example
    if 'AUCell_NFKB1(+)' in adata.obs.columns:
        sc.pl.umap(adata, color='AUCell_NFKB1(+)', ax=axes[2], show=False,
                   title='NFKB1 Regulon Activity', cmap='Reds')

    plt.tight_layout()
    plt.savefig(fig_dir / 'trajectory_overview.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: trajectory_overview.png")

    # Plot 2: TF dynamics heatmap
    fig, ax = plt.subplots(figsize=(10, 12))

    pivot_df = results_df.pivot(index='TF', columns='origin', values='log2FC')
    pivot_df.columns = ['Prolif→IL1B+', 'CD16 Mono→IL1B+']

    # Sort by mean absolute change
    pivot_df['mean_abs'] = pivot_df.abs().mean(axis=1)
    pivot_df = pivot_df.sort_values('mean_abs', ascending=False).drop('mean_abs', axis=1)

    # Top 40 TFs
    sns.heatmap(pivot_df.head(40), cmap='RdBu_r', center=0, ax=ax,
                annot=True, fmt='.1f', annot_kws={'size': 8},
                cbar_kws={'label': 'log2FC (IL1B+ Mac / Origin)'})
    ax.set_title('TF Regulon Changes Along Two Trajectories to IL1B+ Mac')
    ax.set_ylabel('Transcription Factor Regulon')
    plt.tight_layout()
    plt.savefig(fig_dir / 'tf_dynamics_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: tf_dynamics_heatmap.png")

    # Plot 3: Scatter plot comparing two trajectories
    fig, ax = plt.subplots(figsize=(8, 8))

    compare_df = results_df.pivot(index='TF', columns='origin', values='log2FC')
    compare_df.columns = ['x', 'y']

    # Color by category
    colors = []
    for tf in compare_df.index:
        if tf in shared_up or tf in shared_down:
            colors.append('#4DAF4A')  # Green - shared
        elif tf in unique_up_prolif or tf in unique_down_prolif:
            colors.append('#E69F00')  # Orange - prolif specific
        elif tf in unique_up_mono or tf in unique_down_mono:
            colors.append('#377EB8')  # Blue - mono specific
        else:
            colors.append('#CCCCCC')  # Gray - not significant

    ax.scatter(compare_df['x'], compare_df['y'], c=colors, alpha=0.7, s=50)

    # Add diagonal
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', alpha=0.3)
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axvline(0, color='gray', linewidth=0.5)

    # Label key TFs
    key_tfs_to_label = ['FOS(+)', 'FOSB(+)', 'NFKB1(+)', 'MAFB(+)', 'E2F8(+)', 'MAF(+)']
    for tf in key_tfs_to_label:
        if tf in compare_df.index:
            x, y = compare_df.loc[tf, 'x'], compare_df.loc[tf, 'y']
            ax.annotate(tf, (x, y), fontsize=8, alpha=0.8)

    ax.set_xlabel('log2FC: Prolif → IL1B+ Mac')
    ax.set_ylabel('log2FC: CD16 Mono → IL1B+ Mac')
    ax.set_title('TF Changes: Comparing Two Trajectories\n(Green=Shared, Orange=Prolif-specific, Blue=Mono-specific)')

    plt.tight_layout()
    plt.savefig(fig_dir / 'trajectory_comparison_scatter.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: trajectory_comparison_scatter.png")

    # Plot 4: Key inflammatory TFs boxplot
    key_tfs = ['FOS(+)', 'FOSB(+)', 'JUN(+)', 'JUNB(+)', 'ATF3(+)',
               'NFKB1(+)', 'ETS2(+)', 'MAFB(+)', 'MAF(+)', 'MITF(+)']

    available_key_tfs = [tf for tf in key_tfs if f'AUCell_{tf}' in adata.obs.columns]

    if available_key_tfs:
        fig, axes = plt.subplots(2, 5, figsize=(20, 8))
        axes = axes.flatten()

        # Subset to trajectory cells
        trajectory_cells = adata.obs[STATE_COL].isin(INITIAL_STATES + [TERMINAL_STATE])
        adata_traj = adata[trajectory_cells].copy()

        # Rename for plotting
        state_order = [INITIAL_STATES[0], INITIAL_STATES[1], TERMINAL_STATE]
        short_names = {'C5_Mac_Prolif_MKI67': 'Prolif',
                       'C4_Mono_Alternative_CD16': 'CD16 Mono',
                       'C3_Mac_Inflam_IL1B': 'IL1B+ Mac'}
        adata_traj.obs['state_short'] = adata_traj.obs[STATE_COL].map(short_names)

        for i, tf in enumerate(available_key_tfs[:10]):
            ax = axes[i]
            tf_col = f'AUCell_{tf}'

            data = []
            for state in state_order:
                vals = adata_traj.obs.loc[adata_traj.obs[STATE_COL] == state, tf_col].values
                data.append(vals)

            parts = ax.violinplot(data, positions=[1, 2, 3], showmeans=True, showmedians=False)
            for pc, color in zip(parts['bodies'], ['#E69F00', '#377EB8', '#D55E00']):
                pc.set_facecolor(color)
                pc.set_alpha(0.7)

            ax.set_xticks([1, 2, 3])
            ax.set_xticklabels(['Prolif', 'CD16\nMono', 'IL1B+\nMac'], fontsize=9)
            ax.set_title(tf, fontsize=10)
            ax.set_ylabel('AUCell Score')

        plt.suptitle('Key TF Regulon Activity Across Trajectory States', fontsize=12, y=1.02)
        plt.tight_layout()
        plt.savefig(fig_dir / 'key_tf_violins.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   Saved: key_tf_violins.png")

    # Save results
    results_df.to_csv(output_dir / 'tf_trajectory_dynamics.csv', index=False)
    print(f"   Saved: tf_trajectory_dynamics.csv")

    # Save summary
    summary = {
        'shared_upregulated': list(shared_up),
        'shared_downregulated': list(shared_down),
        'prolif_specific_up': list(unique_up_prolif),
        'prolif_specific_down': list(unique_down_prolif),
        'mono_specific_up': list(unique_up_mono),
        'mono_specific_down': list(unique_down_mono)
    }
    pd.DataFrame(dict([(k, pd.Series(v)) for k, v in summary.items()])).to_csv(
        output_dir / 'shared_vs_specific_tfs.csv', index=False)
    print(f"   Saved: shared_vs_specific_tfs.csv")

    # Save adata
    adata.write_h5ad(output_dir / 'MoMac_trajectory_annotated.h5ad')
    print(f"   Saved: MoMac_trajectory_annotated.h5ad")

    # Summary
    print()
    print("=" * 80)
    print("TF Trajectory Analysis Complete!")
    print("=" * 80)
    print(f"\nResults: {output_dir}")
    print()
    print("KEY FINDINGS:")
    print("-" * 40)
    print(f"Two independent origins → IL1B+ Macrophages:")
    print(f"  1. C5_Mac_Prolif_MKI67 (proliferating progenitors)")
    print(f"  2. C4_Mono_Alternative_CD16 (CD16+ monocytes)")
    print()
    print(f"SHARED TF changes (both trajectories):")
    print(f"  Upregulated in IL1B+ Mac: {', '.join(sorted(shared_up)) if shared_up else 'None'}")
    print(f"  Downregulated in IL1B+ Mac: {', '.join(sorted(shared_down)) if shared_down else 'None'}")
    print()
    print(f"TRAJECTORY-SPECIFIC TF changes:")
    print(f"  Prolif→IL1B+: {len(unique_up_prolif)} up, {len(unique_down_prolif)} down")
    print(f"  CD16 Mono→IL1B+: {len(unique_up_mono)} up, {len(unique_down_mono)} down")
    print()

if __name__ == "__main__":
    main()
