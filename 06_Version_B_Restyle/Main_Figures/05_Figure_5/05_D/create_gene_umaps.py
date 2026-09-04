#!/usr/bin/env python3
"""
Figure 5 panel C, RESTYLED (Version B) - TNF, IL1B, IL6 and IL1A expression
UMAPs for MoMac cells.

Version A is
`03_Revised_Panels/Main_Figures/05_Figure_5/05_D/create_gene_umaps.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by AREA,
and the save is `style.save_panel`.

The four genes, the read from `.raw` when it is present, and the 2nd/98th
percentile colour limits are Version A's, unchanged. Expression is read from
`adata.raw`, exactly as Version A reads it; that is not touched.

  printed panel  Figure 5 C        (PROVENANCE.csv - the directory is "05_D";
                 do NOT read the directory as the letter)
  printed rect   43.7 x 38.4 mm    (panel_rects.csv)
  Version B box  84.0 x 72.0 mm

    The same 2x2-with-colorbars layout as panel B, and sized to match it: four
    maps, each with a colorbar at 7 pt and an italic gene name at 8 pt over it.

MARK
    Version A drew 20 x 16.8 cm at SCALE = 4 and set its smallest body type -
    the colorbar tick labels - at `5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2                                     = 0.1225

    The only non-type size the script sets is the scatter marker area, `s=1`,
    which becomes `1 * AREA`. Colorbar tick widths/lengths and the colorbar
    outline width are NOT rescaled - those are style, and cnsplots sets them.

Layout note: as panel B, `plt.tight_layout()` is kept rather than
`style.margins_mm`, because the four colorbars are outside the gridspec.
`style.overflow_mm` still has the last word.

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
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *
import panel_style_cns as style

BASE_DIR = Path(__file__).parent

SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 5.0                  # Version A's smallest body type, before * SCALE

PRINTED_MM = (43.7, 38.4)       # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 84.0, 72.0

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

    ax.set_title(gene, style='italic')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(''); ax.set_ylabel('')

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.02)


def main():
    global AREA
    family = style.apply()
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

    plt.tight_layout()

    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, BASE_DIR / 'tnf_il1b_il6_il1a_cytokines')
    print(f"Saved: tnf_il1b_il6_il1a_cytokines.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
