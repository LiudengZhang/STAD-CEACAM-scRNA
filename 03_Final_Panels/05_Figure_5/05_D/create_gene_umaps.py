#!/usr/bin/env python3
"""
Figure 5 panel C - TNF, IL1B, IL6 and IL1A expression UMAPs for MoMac cells.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed.

  printed panel  Figure 5 C       (PROVENANCE.csv - the directory is "05_D";
                                   do NOT read the directory as the letter)

The four genes, the read from `.raw` when it is present, and the 2nd/98th
percentile colour limits are unchanged.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the colorbar tick labels - at 5 * SCALE, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The only non-type size the script sets is the scatter marker area, `s=1`,
    which becomes `1 * AREA`. Colorbar tick widths and lengths and the colorbar
    outline width are not rescaled - those are style, and cnsplots sets them.

Layout: as panel B, `tight_layout` is kept rather than `style.fit_margins`,
because the four colorbars are made by `plt.colorbar(ax=...)` and live outside
the figure's gridspec, where `subplots_adjust` does not reach them. It is a
layout algorithm rather than a canvas rescale, so the 1:1 relationship is
untouched, and `style.overflow_mm` and `style.letter_clear` still have the last
word.

The gene name stays italic: that is a nomenclature convention, not type
styling, so `style='italic'` is kept while the explicit point size goes.
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
import panel_style_cns as style
import slots

BASE_DIR = Path(__file__).parent

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier smallest body type, before * SCALE

PANEL_LETTER = "C"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(5, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(5, PANEL_LETTER)

# Set in main() once the style is applied; used by create_gene_umap().
AREA = None

GENES = ['TNF', 'IL1B', 'IL6', 'IL1A']


def create_gene_umap(adata, gene, ax):
    if gene not in adata.var_names and (adata.raw is None or gene not in adata.raw.var_names):
        ax.text(0.5, 0.5, f'{gene}\n(not found)', ha='center', va='center', transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        return

    umap = adata.obsm['X_umap']
    if adata.raw is not None and gene in adata.raw.var_names:
        gene_idx = list(adata.raw.var_names).index(gene)
        expression = adata.raw.X[:, gene_idx]
    else:
        gene_idx = list(adata.var_names).index(gene)
        expression = adata.X[:, gene_idx]

    if hasattr(expression, 'toarray'):
        expression = expression.toarray().flatten()
    else:
        expression = np.array(expression).flatten()

    scatter = ax.scatter(umap[:, 0], umap[:, 1], c=expression, cmap='Reds',
                         s=1 * AREA, alpha=0.8,
                         vmin=np.percentile(expression, 2),
                         vmax=np.percentile(expression, 98))
    # One glyph per cell in the SVG makes the panel tens of megabytes, which
    # pushes the assembled page past the size at which the assembler stops
    # compositing vector and starts rasterising whole panels - taking the text
    # with it. Rasterising the point cloud alone keeps every string editable.
    scatter.set_rasterized(True)

    ax.set_title(gene, style='italic')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(''); ax.set_ylabel('')

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.02)


def main():
    global AREA
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("Loading MoMac data...")
    adata = sc.read_h5ad(MOMAC_H5AD)
    print(f"  Loaded {adata.n_obs} cells")

    if 'X_umap' not in adata.obsm:
        print("  Computing UMAP...")
        sc.pp.neighbors(adata, use_rep='X_pca')
        sc.tl.umap(adata)

    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, 2, 2)
    axes = axes.flatten()

    for idx, gene in enumerate(GENES):
        create_gene_umap(adata, gene, axes[idx])

    # The panel letter is drawn over the panel's top-left corner by the
    # assembler, so that corner is kept free. Reserving a left band of the
    # letter cell's width costs less paper here than a top band of its height,
    # which is the choice `style.fit_margins` makes for a panel it can move.
    plt.tight_layout(rect=(LETTER_CELL[0] / PANEL_W_MM, 0.0, 1.0, 1.0))

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    style.save_panel(fig, BASE_DIR / 'tnf_il1b_il6_il1a_cytokines')
    print(f"Saved: tnf_il1b_il6_il1a_cytokines.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
