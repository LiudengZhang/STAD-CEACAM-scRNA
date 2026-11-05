#!/usr/bin/env python3
"""Figure 5B: four Hallmark pathway scores (sc.tl.score_genes) on the MoMac
UMAP, in a 2 x 2 grid, colour limits at the 2nd and 98th percentiles.
Each pathway is scored over its highly variable genes only.

Input:  MOMAC_H5AD, HALLMARK_GMT (00_Config/paths.py)
Output: momac_4pathway_umap.{svg,pdf,png} in this directory
Run:    python create_4pathway_umap_matched.py
"""
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *
from cnsfig import Panel, style
from cnsfig.layout import umap_grid_mm

SCALE = 4                       # marker-size scaling constants
SMALL_PT = 5.0

PANEL = Panel(5, "B", __file__, small_pt=SMALL_PT)
#: Map size and grid offsets in mm, shared with panel C so the maps match.
MAP_MM, GRID_LEFT_MM, GRID_TOP_MM = 12.8, 4.3, 1.0

PATHWAYS = [
    ('TNF-alpha Signaling via NF-kB', 'TNFα/NF-κB'),
    ('Inflammatory Response', 'Inflammation'),
    ('Interferon Gamma Response', 'IFN-γ'),
    ('IL-6/JAK/STAT3 Signaling', 'IL-6/STAT3'),
]


def _read_gmt(path):
    """Read a GMT file into {set name: [genes]}."""
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def main():
    PANEL.apply()
    AREA = (style.tick_pt() / (SMALL_PT * SCALE)) ** 2

    adata = sc.read_h5ad(MOMAC_H5AD)
    gene_sets = _read_gmt(HALLMARK_GMT)

    # Score only the highly variable genes of each pathway.
    if "highly_variable" not in adata.var:
        raise SystemExit(
            f"{MOMAC_H5AD} has no var['highly_variable'] - cannot reproduce the "
            f"gene selection; rebuild it with 06_Clean_Data/build_clean_h5ad.py")
    scored_genes = set(adata.var_names[adata.var["highly_variable"].to_numpy(dtype=bool)])

    for pathway_name, display_name in PATHWAYS:
        matching_key = None
        for key in gene_sets.keys():
            if pathway_name.lower().replace('-', ' ').replace('/', ' ') in key.lower().replace('-', ' ').replace('/', ' '):
                matching_key = key
                break
        if matching_key is None:
            continue
        genes = gene_sets[matching_key]
        genes_in_data = [g for g in genes if g in scored_genes]
        if len(genes_in_data) < 5:
            continue
        score_name = pathway_name.replace(' ', '_').replace('-', '_').replace('/', '_')
        sc.tl.score_genes(adata, genes_in_data, score_name=score_name)

    fig = style.figure_mm(PANEL.w_mm, PANEL.h_mm)
    cells = umap_grid_mm(fig, n_rows=2, n_cols=2, map_mm=MAP_MM,
                         left_mm=GRID_LEFT_MM, top_mm=GRID_TOP_MM)
    axes = [ax for ax, _ in cells]
    caxes = [cax for _, cax in cells]

    for idx, (pathway_name, display_name) in enumerate(PATHWAYS):
        ax = axes[idx]
        score_name = pathway_name.replace(' ', '_').replace('-', '_').replace('/', '_')

        if score_name not in adata.obs.columns:
            ax.text(0.5, 0.5, f'{display_name}\n(N/A)', ha='center', va='center', transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])
            continue

        umap = adata.obsm['X_umap']
        scores = adata.obs[score_name].values
        scatter = ax.scatter(umap[:, 0], umap[:, 1], c=scores, cmap='Purples',
                             s=1 * AREA, alpha=0.8, rasterized=True,
                             vmin=np.percentile(scores, 2), vmax=np.percentile(scores, 98))
        ax.set_title(display_name, pad=2.0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(''); ax.set_ylabel('')

        cbar = plt.colorbar(scatter, cax=caxes[idx])
        cbar.ax.tick_params(labelsize=style.tick_pt(), width=style.RULE_PT, length=1.5, pad=1.0)
        cbar.outline.set_linewidth(style.RULE_PT)

    PANEL.finish(fig, 'momac_4pathway_umap')


if __name__ == "__main__":
    main()
