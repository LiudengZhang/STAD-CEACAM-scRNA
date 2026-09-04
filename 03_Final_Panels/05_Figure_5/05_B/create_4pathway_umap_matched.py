#!/usr/bin/env python3
"""
Panel B: 2×2 Pathway Score UMAPs for MoMac cells.
4× scaling method for Nature Cancer.
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
from shared.figure_config import use_panel_style

BASE_DIR = Path(__file__).parent

# 4× scaling
DPI = 300
SCALE = 4
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 4.5 * SCALE
PANEL_HEIGHT_CM = 4.2 * SCALE

PATHWAYS = [
    ('TNF-alpha Signaling via NF-kB', 'TNFα/NF-κB'),
    ('Inflammatory Response', 'Inflammation'),
    ('Interferon Gamma Response', 'IFN-γ'),
    ('IL-6/JAK/STAT3 Signaling', 'IL-6/STAT3'),
]


def _read_gmt(path):
    """The pinned Hallmark sets, in the shape gseapy.get_library returns."""
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def main():
    use_panel_style(font_pt=7)

    print("Loading MoMac data...")
    adata = sc.read_h5ad(MOMAC_H5AD)
    print(f"  Loaded {adata.n_obs} cells")

    print("Fetching Hallmark gene sets...")
    gene_sets = _read_gmt(HALLMARK_GMT)

    # The published panel scored only the highly-variable genes of each pathway
    # - 33 of the 87 in IL-6/JAK/STAT3, 106 of the 200 in TNF-alpha via NF-kB -
    # because the working file's .X carries just that subset and the filter
    # below tested membership of .X.var_names. The clean deposit promotes
    # .raw.X to .X, so there var_names is every gene and the same filter would
    # silently widen the pathway and redraw the panel. var['highly_variable']
    # names the subset in both files, so the selection no longer depends on
    # which one is loaded.
    if "highly_variable" not in adata.var:
        raise SystemExit(
            f"{MOMAC_H5AD} has no var['highly_variable'] - cannot reproduce the "
            f"published gene selection; rebuild it with 06_Clean_Data/build_clean_h5ad.py")
    scored_genes = set(adata.var_names[adata.var["highly_variable"].to_numpy(dtype=bool)])

    print("Computing pathway scores...")
    for pathway_name, display_name in PATHWAYS:
        matching_key = None
        for key in gene_sets.keys():
            if pathway_name.lower().replace('-', ' ').replace('/', ' ') in key.lower().replace('-', ' ').replace('/', ' '):
                matching_key = key
                break
        if matching_key is None:
            print(f"  Warning: Could not find {pathway_name}")
            continue
        genes = gene_sets[matching_key]
        genes_in_data = [g for g in genes if g in scored_genes]
        if len(genes_in_data) < 5:
            continue
        score_name = pathway_name.replace(' ', '_').replace('-', '_').replace('/', '_')
        sc.tl.score_genes(adata, genes_in_data, score_name=score_name)
        print(f"  {pathway_name}: {len(genes_in_data)} genes")

    fig, axes = plt.subplots(2, 2, figsize=(PANEL_WIDTH_CM * CM_TO_INCH, PANEL_HEIGHT_CM * CM_TO_INCH))
    axes = axes.flatten()

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
                             s=1, alpha=0.8, rasterized=True,
                             vmin=np.percentile(scores, 2), vmax=np.percentile(scores, 98))
        ax.set_title(display_name, fontsize=7 * SCALE)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(''); ax.set_ylabel('')

        cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.02)
        cbar.ax.tick_params(labelsize=5 * SCALE, width=1.0, length=3)
        cbar.outline.set_linewidth(1.0)

    plt.tight_layout()

    output = BASE_DIR / 'momac_4pathway_umap.png'
    plt.savefig(output, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.svg'), dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.pdf'), format='pdf', dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
