#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
CellRank Standard Analysis - Dual Origin to IL1B+ Mac
======================================================
Standard CellRank workflow with:
  - 2 Initial states: C5_Mac_Prolif_MKI67, C4_Mono_Alternative_CD16
  - 1 Terminal state: C3_Mac_Inflam_IL1B
"""

# =========================================================================
# SCIPY COMPATIBILITY FIX - Must be done BEFORE importing cellrank
# scipy.sparse.linalg.gmres renamed 'tol' to 'rtol' and changed 'atol' handling
# =========================================================================
import scipy.sparse.linalg as _ssl
_original_gmres = _ssl.gmres

def _patched_gmres(*args, **kwargs):
    """Wrapper to handle scipy compatibility issues in gmres"""
    if 'tol' in kwargs and 'rtol' not in kwargs:
        kwargs['rtol'] = kwargs.pop('tol')
    # Handle atol='legacy' which is no longer valid
    if kwargs.get('atol') == 'legacy' or kwargs.get('atol') is None:
        kwargs['atol'] = 0.0  # Set to 0 as a reasonable default
    return _original_gmres(*args, **kwargs)

_ssl.gmres = _patched_gmres

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def main():

    base_dir = Path('/path/to/upstream-pipeline')
    project_dir = base_dir / '02_Playground/02_Cell-Cell_Interaction/01_Analysis/09_SCENIC/17_Full_MoMac_SCENIC'

    input_file = project_dir / 'Results/MoMac_12sample_post_stomach.h5ad'
    aucell_file = project_dir / 'Results/aucell_matrix.csv'
    output_dir = project_dir / 'Results/CellRank_DualOrigin'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / 'figures'
    fig_dir.mkdir(exist_ok=True)

    STATE_COL = 'minor_cell_state'
    INITIAL_STATES = ['C5_Mac_Prolif_MKI67', 'C4_Mono_Alternative_CD16']
    TERMINAL_STATE = 'C3_Mac_Inflam_IL1B'

    adata = sc.read_h5ad(input_file)
    for state in INITIAL_STATES + [TERMINAL_STATE]:
        n = (adata.obs[STATE_COL] == state).sum()
    aucell = pd.read_csv(aucell_file, index_col=0)

    for col in aucell.columns:
        adata.obs[f'AUCell_{col}'] = aucell.loc[adata.obs_names, col].values

    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

    if 'X_umap' not in adata.obsm:
        sc.tl.umap(adata)

    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=30)

    sc.tl.diffmap(adata, n_comps=10, random_state=42)

    initial_cells = adata.obs[adata.obs[STATE_COL].isin(INITIAL_STATES)].index

    n_genes = (adata[initial_cells].X > 0).sum(axis=1)
    if hasattr(n_genes, 'A1'):
        n_genes = n_genes.A1
    root_cell_idx = np.argmax(n_genes)
    root_cell = initial_cells[root_cell_idx]
    root_idx_global = list(adata.obs_names).index(root_cell)
    adata.uns['iroot'] = root_idx_global

    sc.tl.dpt(adata, n_dcs=5)

    from cellrank.kernels import PseudotimeKernel

    pk = PseudotimeKernel(adata, time_key='dpt_pseudotime')
    pk.compute_transition_matrix(threshold_scheme='soft', nu=0.5)

    from cellrank.estimators import GPCCA

    g = GPCCA(pk)

    g.compute_schur(n_components=20)
    g.compute_macrostates(n_states=6, cluster_key=STATE_COL)
    terminal_macrostates = []
    for ms in g.macrostates.cat.categories:
        if TERMINAL_STATE in ms or 'IL1B' in ms:
            terminal_macrostates.append(ms)

    if not terminal_macrostates:
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

    g.set_terminal_states(terminal_macrostates)

    g.compute_fate_probabilities()

    fate_probs = g.fate_probabilities
    fate_names = fate_probs.names if hasattr(fate_probs, 'names') else list(range(fate_probs.shape[1]))
    for i, name in enumerate(fate_names):
        adata.obs[f'fate_prob_{name}'] = fate_probs[:, i]
    try:
        g.compute_absorption_times()
        if hasattr(g, 'absorption_times') and g.absorption_times is not None:
            abs_times = g.absorption_times
            if hasattr(abs_times, 'values'):
                adata.obs['absorption_time'] = abs_times.values.flatten()
            else:
                adata.obs['absorption_time'] = np.array(abs_times).flatten()
    except Exception as e:
        adata.obs['absorption_time'] = adata.obs['dpt_pseudotime']

    from scipy import stats

    aucell_cols = [c for c in adata.obs.columns if c.startswith('AUCell_')]

    pseudotime_vals = adata.obs['dpt_pseudotime'].values

    correlations = []
    for tf_col in aucell_cols:
        tf_vals = adata.obs[tf_col].values
        tf_name = tf_col.replace('AUCell_', '')

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

    early_tfs = corr_df[corr_df['spearman_r'] < 0].head(5)
    for _, row in early_tfs.iterrows():
        sig = '***' if row['pvalue'] < 0.001 else '**' if row['pvalue'] < 0.01 else '*'

    late_tfs = corr_df[corr_df['spearman_r'] > 0].head(5)
    for _, row in late_tfs.iterrows():
        sig = '***' if row['pvalue'] < 0.001 else '**' if row['pvalue'] < 0.01 else '*'

    fate_cols = [c for c in adata.obs.columns if c.startswith('fate_prob_')]

    adata.obs['trajectory_role'] = 'Other'
    adata.obs.loc[adata.obs[STATE_COL] == INITIAL_STATES[0], 'trajectory_role'] = 'Origin: Prolif'
    adata.obs.loc[adata.obs[STATE_COL] == INITIAL_STATES[1], 'trajectory_role'] = 'Origin: CD16 Mono'
    adata.obs.loc[adata.obs[STATE_COL] == TERMINAL_STATE, 'trajectory_role'] = 'Terminal: IL1B+ Mac'


    corr_df.to_csv(output_dir / 'tf_trajectory_correlations.csv', index=False)

    traj_cols = [STATE_COL, 'trajectory_role', 'dpt_pseudotime'] + fate_cols
    if 'absorption_time' in adata.obs.columns:
        traj_cols.append('absorption_time')
    traj_df = adata.obs[traj_cols].copy()
    traj_df.to_csv(output_dir / 'trajectory_data.csv')

    adata.write_h5ad(output_dir / 'MoMac_cellrank_dual_origin.h5ad')

if __name__ == "__main__":
    main()
