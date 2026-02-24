#!/usr/bin/env python3
"""
CellRank Standard Analysis - Dual Origin to IL1B+ Mac
======================================================
Standard CellRank workflow with:
  - 2 Initial states: C5_Mac_Prolif_MKI67, C4_Mono_Alternative_CD16
  - 1 Terminal state: C3_Mac_Inflam_IL1B

Uses ConnectivityKernel (no RNA velocity needed).
Note: Will be slow without petsc4py but will work.
"""

# =========================================================================
# SCIPY COMPATIBILITY FIX - Must be done BEFORE importing cellrank
# scipy.sparse.linalg.gmres renamed 'tol' to 'rtol' and changed 'atol' handling
# =========================================================================
import scipy.sparse.linalg as _ssl
_original_gmres = _ssl.gmres

def _patched_gmres(*args, **kwargs):
    """Wrapper to handle scipy compatibility issues in gmres"""
    # Handle tol -> rtol rename
    if 'tol' in kwargs and 'rtol' not in kwargs:
        kwargs['rtol'] = kwargs.pop('tol')
    # Handle atol='legacy' which is no longer valid
    if kwargs.get('atol') == 'legacy' or kwargs.get('atol') is None:
        kwargs['atol'] = 0.0  # Set to 0 as a reasonable default
    return _original_gmres(*args, **kwargs)

_ssl.gmres = _patched_gmres
print("[PATCH] Applied scipy.sparse.linalg.gmres compatibility fix (tol+atol)")
# =========================================================================

import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def main():
    print("=" * 80)
    print("CellRank Standard Analysis - Dual Origin to IL1B+ Mac")
    print("=" * 80)
    print()
    print("Note: This may take 1-2 hours without petsc4py. Please be patient.")
    print()

    # Central config
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
    from paths import SCENIC_RESULTS_DIR

    input_file = SCENIC_RESULTS_DIR / 'MoMac_12sample_post_stomach.h5ad'
    aucell_file = SCENIC_RESULTS_DIR / 'aucell_matrix.csv'
    output_dir = SCENIC_RESULTS_DIR / 'CellRank_DualOrigin'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / 'figures'
    fig_dir.mkdir(exist_ok=True)

    # Define biological states
    STATE_COL = 'minor_cell_state'
    INITIAL_STATES = ['C5_Mac_Prolif_MKI67', 'C4_Mono_Alternative_CD16']
    TERMINAL_STATE = 'C3_Mac_Inflam_IL1B'

    # =========================================================================
    # Step 1: Load data
    # =========================================================================
    print("[1/7] Loading MoMac data...")
    adata = sc.read_h5ad(input_file)
    print(f"   Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # Check populations
    print(f"\n   Population sizes:")
    for state in INITIAL_STATES + [TERMINAL_STATE]:
        n = (adata.obs[STATE_COL] == state).sum()
        print(f"      {state}: {n} cells")
    print()

    # =========================================================================
    # Step 2: Load SCENIC AUCell scores
    # =========================================================================
    print("[2/7] Loading SCENIC AUCell scores...")
    aucell = pd.read_csv(aucell_file, index_col=0)

    # Add to adata.obs
    for col in aucell.columns:
        adata.obs[f'AUCell_{col}'] = aucell.loc[adata.obs_names, col].values
    print(f"   Added {len(aucell.columns)} regulon scores")
    print()

    # =========================================================================
    # Step 3: Ensure proper preprocessing for CellRank
    # =========================================================================
    print("[3/7] Preparing data for CellRank...")

    # Recompute neighbors if needed
    if 'neighbors' not in adata.uns:
        print("   Computing neighbors...")
        sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

    # Ensure UMAP exists
    if 'X_umap' not in adata.obsm:
        print("   Computing UMAP...")
        sc.tl.umap(adata)

    print("   Data ready")
    print()

    # =========================================================================
    # Step 4: Build CellRank kernel (PseudotimeKernel - no velocity needed)
    # =========================================================================
    print("[4/7] Building CellRank PseudotimeKernel...")
    print("   (Using diffusion pseudotime, no RNA velocity required)")

    # Compute diffusion pseudotime first
    # Recompute neighbors with higher n_neighbors for better connectivity
    print("   Recomputing neighbors for better connectivity...")
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=30)

    print("   Computing diffusion components...")
    # Use fewer components and random_state for reproducibility
    sc.tl.diffmap(adata, n_comps=10, random_state=42)

    # Find root cell - use cell from one of the initial states with highest connectivity
    print("   Finding root cell from initial states...")
    initial_cells = adata.obs[adata.obs[STATE_COL].isin(INITIAL_STATES)].index

    # Use the cell with highest number of expressed genes as root (proxy for undifferentiated)
    n_genes = (adata[initial_cells].X > 0).sum(axis=1)
    if hasattr(n_genes, 'A1'):
        n_genes = n_genes.A1
    root_cell_idx = np.argmax(n_genes)
    root_cell = initial_cells[root_cell_idx]
    root_idx_global = list(adata.obs_names).index(root_cell)
    adata.uns['iroot'] = root_idx_global

    print(f"   Root cell: {root_cell} (state: {adata.obs.loc[root_cell, STATE_COL]})")

    # Compute DPT
    print("   Computing diffusion pseudotime...")
    sc.tl.dpt(adata, n_dcs=5)
    print(f"   DPT range: {adata.obs['dpt_pseudotime'].min():.3f} - {adata.obs['dpt_pseudotime'].max():.3f}")

    # Use PseudotimeKernel
    from cellrank.kernels import PseudotimeKernel

    pk = PseudotimeKernel(adata, time_key='dpt_pseudotime')
    pk.compute_transition_matrix(threshold_scheme='soft', nu=0.5)
    print("   Transition matrix computed")
    print()

    # =========================================================================
    # Step 5: GPCCA for macrostate identification
    # =========================================================================
    print("[5/7] Running GPCCA estimator...")
    print("   WARNING: This step is slow without petsc4py. May take 30-60 minutes.")
    print("   Computing Schur decomposition...")

    from cellrank.estimators import GPCCA

    g = GPCCA(pk)

    # Compute Schur decomposition (this is the slow part)
    g.compute_schur(n_components=20)
    print("   Schur decomposition complete")

    # Compute macrostates
    print("   Computing macrostates...")
    g.compute_macrostates(n_states=6, cluster_key=STATE_COL)
    print(f"   Macrostates identified: {list(g.macrostates.cat.categories)}")

    # Set terminal states (our known biological terminal)
    print(f"\n   Setting terminal state: {TERMINAL_STATE}")
    # Find which macrostate corresponds to our terminal
    terminal_macrostates = []
    for ms in g.macrostates.cat.categories:
        if TERMINAL_STATE in ms or 'IL1B' in ms:
            terminal_macrostates.append(ms)

    if not terminal_macrostates:
        # Use the macrostate with most IL1B+ cells
        print("   Auto-detecting terminal macrostate based on IL1B+ cell enrichment...")
        il1b_cells = adata.obs[adata.obs[STATE_COL] == TERMINAL_STATE].index
        best_ms = None
        best_overlap = 0
        for ms in g.macrostates.cat.categories:
            ms_cells = adata.obs[g.macrostates == ms].index
            overlap = len(set(il1b_cells) & set(ms_cells))
            if overlap > best_overlap:
                best_overlap = overlap
                best_ms = ms
        terminal_macrostates = [best_ms]
        print(f"   Selected: {best_ms} (contains {best_overlap} IL1B+ cells)")

    g.set_terminal_states(terminal_macrostates)
    print(f"   Terminal states set: {terminal_macrostates}")
    print()

    # =========================================================================
    # Step 6: Compute fate probabilities
    # =========================================================================
    print("[6/7] Computing fate probabilities toward terminal state...")

    # Compute fate probabilities (scipy gmres is patched at top of script)
    g.compute_fate_probabilities()
    print("   Fate probabilities computed")

    # Extract fate probabilities (Lineage object uses .names for column names)
    fate_probs = g.fate_probabilities
    fate_names = fate_probs.names if hasattr(fate_probs, 'names') else list(range(fate_probs.shape[1]))
    for i, name in enumerate(fate_names):
        adata.obs[f'fate_prob_{name}'] = fate_probs[:, i]
    print(f"   Fate probability columns: {list(fate_names)}")

    # Compute absorption times (with timeout handling)
    print("   Computing absorption times...")
    try:
        g.compute_absorption_times()
        if hasattr(g, 'absorption_times') and g.absorption_times is not None:
            abs_times = g.absorption_times
            if hasattr(abs_times, 'values'):
                adata.obs['absorption_time'] = abs_times.values.flatten()
            else:
                adata.obs['absorption_time'] = np.array(abs_times).flatten()
            print("   Absorption times computed")
    except Exception as e:
        print(f"   WARNING: Could not compute absorption times: {e}")
        # Use pseudotime as fallback for absorption time
        adata.obs['absorption_time'] = adata.obs['dpt_pseudotime']
        print("   Using pseudotime as absorption time fallback")
    print()

    # =========================================================================
    # Step 7: Correlate TF activity with pseudotime (trajectory progression)
    # =========================================================================
    print("[7/7] Correlating TF regulon activity with trajectory (pseudotime)...")
    print("   Note: Using pseudotime because single terminal state gives fate_prob=1 for all cells")

    from scipy import stats

    aucell_cols = [c for c in adata.obs.columns if c.startswith('AUCell_')]

    # Use pseudotime as the trajectory measure (more informative than fate probability when single terminal)
    pseudotime_vals = adata.obs['dpt_pseudotime'].values

    correlations = []
    for tf_col in aucell_cols:
        tf_vals = adata.obs[tf_col].values
        tf_name = tf_col.replace('AUCell_', '')

        # Spearman correlation with pseudotime
        r, p = stats.spearmanr(pseudotime_vals, tf_vals)

        correlations.append({
            'fate': 'IL1B_trajectory',  # trajectory toward IL1B+ Mac terminal state
            'TF': tf_name,
            'spearman_r': r,
            'pvalue': p,
            'direction': 'increases' if r > 0 else 'decreases'  # increases/decreases along trajectory
        })

    corr_df = pd.DataFrame(correlations)
    corr_df = corr_df.sort_values('spearman_r', key=abs, ascending=False)

    # Print top TFs
    print("\n   Top TFs correlated with trajectory toward IL1B+ Mac:")
    print("\n   EARLY TFs (decrease along trajectory - high in initial states):")
    early_tfs = corr_df[corr_df['spearman_r'] < 0].head(5)
    for _, row in early_tfs.iterrows():
        sig = '***' if row['pvalue'] < 0.001 else '**' if row['pvalue'] < 0.01 else '*'
        print(f"      {row['TF']}: r={row['spearman_r']:.3f} {sig}")

    print("\n   LATE TFs (increase along trajectory - high in terminal IL1B+ Mac):")
    late_tfs = corr_df[corr_df['spearman_r'] > 0].head(5)
    for _, row in late_tfs.iterrows():
        sig = '***' if row['pvalue'] < 0.001 else '**' if row['pvalue'] < 0.01 else '*'
        print(f"      {row['TF']}: r={row['spearman_r']:.3f} {sig}")
    print()

    # =========================================================================
    # Generate visualizations
    # =========================================================================
    print("Generating visualizations...")

    fate_cols = [c for c in adata.obs.columns if c.startswith('fate_prob_')]

    # Plot 1: UMAP with trajectory annotation
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
               title='Dual Origin → IL1B+ Mac', palette=colors)

    # Pseudotime (more informative than fate probability for single terminal)
    sc.pl.umap(adata, color='dpt_pseudotime', ax=axes[2], show=False,
               title='Pseudotime (0=Initial → 1=Terminal)', cmap='viridis')

    plt.tight_layout()
    plt.savefig(fig_dir / 'trajectory_overview.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: trajectory_overview.png")

    # Plot 2: Pseudotime by origin population
    fig, ax = plt.subplots(figsize=(10, 6))

    plot_data = []
    labels = []
    for state in INITIAL_STATES + [TERMINAL_STATE]:
        vals = adata.obs.loc[adata.obs[STATE_COL] == state, 'dpt_pseudotime'].values
        plot_data.append(vals)
        labels.append(state.replace('C5_Mac_Prolif_', '').replace('C4_Mono_Alternative_', '').replace('C3_Mac_Inflam_', ''))

    parts = ax.violinplot(plot_data, positions=range(len(labels)), showmeans=True)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(['#E69F00', '#377EB8', '#D55E00'][i])
        pc.set_alpha(0.7)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylabel('Pseudotime')
    ax.set_title(f'Pseudotime by Population (0=Initial → 1={TERMINAL_STATE})')

    plt.tight_layout()
    plt.savefig(fig_dir / 'pseudotime_by_origin.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: pseudotime_by_origin.png")

    # Plot 3: TF-Trajectory correlation heatmap
    # Sort by correlation value
    corr_sorted = corr_df.sort_values('spearman_r')

    fig, ax = plt.subplots(figsize=(8, 12))
    y_pos = range(len(corr_sorted))
    colors_bar = ['#D55E00' if r > 0 else '#0072B2' for r in corr_sorted['spearman_r']]
    ax.barh(y_pos, corr_sorted['spearman_r'], color=colors_bar)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(corr_sorted['TF'], fontsize=8)
    ax.set_xlabel('Spearman Correlation with Pseudotime')
    ax.set_title('TF Regulon Activity Along Trajectory\n(Blue=Early TFs, Red=Late TFs)')
    ax.axvline(x=0, color='black', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(fig_dir / 'tf_trajectory_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: tf_trajectory_heatmap.png")

    # Plot 4: Key inflammatory TFs vs pseudotime
    key_tfs = ['FOS(+)', 'FOSB(+)', 'JUN(+)', 'NFKB1(+)', 'NFKB2(+)',
               'ETS2(+)', 'MAFB(+)', 'MAF(+)']
    available_tfs = [tf for tf in key_tfs if f'AUCell_{tf}' in adata.obs.columns]

    if available_tfs:
        fig, axes = plt.subplots(2, 4, figsize=(16, 8))
        axes = axes.flatten()

        for i, tf in enumerate(available_tfs[:8]):
            ax = axes[i]
            tf_col = f'AUCell_{tf}'

            x = adata.obs['dpt_pseudotime'].values
            y = adata.obs[tf_col].values

            ax.scatter(x, y, alpha=0.1, s=1, c='steelblue')

            # Trend line
            z = np.polyfit(x, y, 2)
            p = np.poly1d(z)
            x_line = np.linspace(x.min(), x.max(), 100)
            ax.plot(x_line, p(x_line), 'r-', linewidth=2)

            # Correlation
            r, pval = stats.spearmanr(x, y)
            ax.set_xlabel('Pseudotime')
            ax.set_ylabel('AUCell Score')
            ax.set_title(f'{tf}\nr={r:.3f}')

        plt.suptitle('Key TF Activity vs Pseudotime (Trajectory)', fontsize=12, y=1.02)
        plt.tight_layout()
        plt.savefig(fig_dir / 'key_tf_vs_pseudotime.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   Saved: key_tf_vs_pseudotime.png")

    # =========================================================================
    # Save results
    # =========================================================================
    print("\nSaving results...")

    # Save TF-trajectory correlations
    corr_df.to_csv(output_dir / 'tf_trajectory_correlations.csv', index=False)
    print(f"   Saved: tf_trajectory_correlations.csv")

    # Save trajectory data (pseudotime + fate probabilities)
    traj_cols = [STATE_COL, 'trajectory_role', 'dpt_pseudotime'] + fate_cols
    if 'absorption_time' in adata.obs.columns:
        traj_cols.append('absorption_time')
    traj_df = adata.obs[traj_cols].copy()
    traj_df.to_csv(output_dir / 'trajectory_data.csv')
    print(f"   Saved: trajectory_data.csv")

    # Save adata
    adata.write_h5ad(output_dir / 'MoMac_cellrank_dual_origin.h5ad')
    print(f"   Saved: MoMac_cellrank_dual_origin.h5ad")

    # =========================================================================
    # Summary
    # =========================================================================
    print()
    print("=" * 80)
    print("CellRank Dual Origin Analysis Complete!")
    print("=" * 80)
    print(f"\nResults: {output_dir}")
    print()
    print("Summary:")
    print(f"  - Initial states: {INITIAL_STATES}")
    print(f"  - Terminal state: {TERMINAL_STATE}")
    print(f"  - Macrostates identified: {list(g.macrostates.cat.categories)}")
    print(f"  - Terminal macrostates: {terminal_macrostates}")
    print()

if __name__ == "__main__":
    main()
