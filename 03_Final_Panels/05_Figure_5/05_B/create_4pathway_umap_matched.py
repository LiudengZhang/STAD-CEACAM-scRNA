#!/usr/bin/env python3
"""
Figure 5 panel B - 2x2 pathway-score UMAPs for MoMac cells.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed.

  printed panel  Figure 5 B       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

The four Hallmark pathways, the pinned GMT, the `highly_variable` gene
selection that reproduces the published scoring, `sc.tl.score_genes` and the
2nd/98th percentile colour limits are unchanged, and so is the comment
explaining the `highly_variable` selection: it records a real trap, not a
preference.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the colorbar tick labels - at 5 * SCALE, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The only non-type size the script sets is the scatter marker area, `s=1`,
    which becomes `1 * AREA`: the same fraction of the type size it was before.
    Colorbar tick widths and lengths and the colorbar outline width are not
    rescaled - those are style, and cnsplots sets them.

The point layer is rasterised and the labels are not. One glyph per cell would
make the panel tens of megabytes and push the assembled page past the size at
which the assembler stops compositing vector and flattens whole panels, taking
their text with them.

Layout: the four colorbars are made by `plt.colorbar(ax=...)`, which creates
axes outside the figure's gridspec, so `subplots_adjust` - and with it
`style.fit_margins` - moves the maps and leaves the colorbars behind.
`tight_layout` is a layout algorithm rather than a canvas rescale, so the 1:1
relationship is untouched, and `style.overflow_mm` and `style.letter_clear`
still have the last word.
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
from cnsfig.layout import umap_grid_mm
import slots

BASE_DIR = Path(__file__).parent

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier smallest body type, before * SCALE

PANEL_LETTER = "B"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(5, PANEL_LETTER)
#: The grid B and C share: a 12.8 mm map, the letter cell and 0.6 mm to its
#: left, 1.0 mm above. The same three numbers in both scripts.
MAP_MM, GRID_LEFT_MM, GRID_TOP_MM = 12.8, 4.3, 1.0
LETTER_CELL = slots.letter_cell_mm(5, PANEL_LETTER)

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
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
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

    # ONE GRID FOR B AND C  (2026-09-14)
    #   Panels B and C are the same drawing twice, and tight_layout gave them
    #   maps of different sizes (12.88 x 14.32 vs 13.23 x 13.31 mm). Both now
    #   call cnsfig.layout.umap_grid_mm with the same numbers, so the maps are
    #   the same square on the page. The slot is C's width for both.
    fig = style.figure_mm(PANEL_W_MM, PANEL_H_MM)
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

    # No tight_layout: every axes is at its millimetres already.

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    style.save_panel(fig, BASE_DIR / 'momac_4pathway_umap')
    print(f"Saved: momac_4pathway_umap.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
