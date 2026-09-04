#!/usr/bin/env python3
"""
Figure 5 panel B, RESTYLED (Version B) - 2x2 pathway-score UMAPs for MoMac
cells.

Version A is
`03_Revised_Panels/Main_Figures/05_Figure_5/05_B/create_4pathway_umap_matched.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK
(areas by AREA), and the save is `style.save_panel`.

The four Hallmark pathways, the pinned GMT, the `highly_variable` gene
selection that reproduces the published scoring, `sc.tl.score_genes`, and the
2nd/98th percentile colour limits are Version A's, unchanged. The comment
explaining the `highly_variable` selection is Version A's too and is kept
verbatim: it records a real trap, not a preference.

  printed panel  Figure 5 B        (PROVENANCE.csv; "05_B" happens to agree,
                 but it was looked up, not inferred)
  printed rect   48.0 x 40.6 mm    (panel_rects.csv)
  Version B box  84.0 x 72.0 mm

    Four UMAPs, each with its own colorbar and a title, in a 2x2 grid: at 7 pt
    tick labels on the colorbars and an 8 pt title over each map, a quadrant
    needs about 40 x 34 mm.

MARK
    Version A drew 18 x 16.8 cm at SCALE = 4 and set its smallest body type -
    the colorbar tick labels - at `5 * SCALE`, so SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35
        AREA = MARK ** 2                                     = 0.1225

    The only non-type size the script sets is the scatter marker area, `s=1`,
    which becomes `1 * AREA`: a 0.35 pt dot, the same fraction of the type size
    it was before. Colorbar tick widths/lengths and the colorbar outline width
    are NOT rescaled - those are style, and cnsplots sets them
    (axes.linewidth 0.5, ticks size 2 width 0.6).

Layout note: this panel keeps Version A's `plt.tight_layout()` instead of
`style.margins_mm`. The four colorbars are made by `plt.colorbar(ax=...)`,
which creates axes outside the figure's gridspec; `subplots_adjust` moves the
maps and leaves the colorbars behind. `tight_layout` is a layout algorithm, not
a canvas rescale, so the 1:1 relationship is untouched, and
`style.overflow_mm` still has the last word.
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

PRINTED_MM = (48.0, 40.6)       # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 84.0, 72.0

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
    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

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

    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, 2, 2)
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
                             s=1 * AREA, alpha=0.8, rasterized=True,
                             vmin=np.percentile(scores, 2), vmax=np.percentile(scores, 98))
        ax.set_title(display_name)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(''); ax.set_ylabel('')

        cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.02)

    plt.tight_layout()

    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, BASE_DIR / 'momac_4pathway_umap')
    print(f"Saved: momac_4pathway_umap.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
