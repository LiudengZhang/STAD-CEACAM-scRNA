#!/usr/bin/env python3
"""
Figure 2 panel F - Milo differential abundance over the pre-treatment
epithelial neighbourhood graph, responders versus non-responders.

RUN THIS IN THE pertpy ENVIRONMENT, NOT THE ONE THE OTHER PANELS USE.

pertpy pulls in scvi, which imports jax; where that import is broken the script
dies before it reads any data. Nothing is wrong with the panel; the two
environments simply disagree about jax. run_all.sh knows this and runs the
script in the environment that carries pertpy.

The type system is read from `00_Config/panel_style_cns.py`. cnsplots is not
installed beside pertpy, so that module falls back to the settings snapshot in
`00_Config/panel_style_rc.json`, which records the same type spec and refuses
to apply itself if a different one is asked for.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. The margins are millimetres of paper and are
set before the neighbourhood graph is drawn, because that call makes the colour
bar and matplotlib sizes a colour bar from the axes box it finds.

  printed panel  Figure 2 F       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the tick labels - at 5 * SCALE. MARK carries the
    non-type point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The spine and tick widths the earlier drawing set are not carried over:
    those are axes furniture, cnsplots has its own settings for them, and
    following the library rather than rescaling the old numbers is the
    standard-methods rule. Nothing they control is a plotted value.

Every cell, every filter, every Milo parameter and every string is the earlier
drawing's. The drawing code is the same code.
"""

import scanpy as sc
import pertpy as pt
import mudata as mu
import os
import warnings
warnings.filterwarnings('ignore')

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402
import slots                                              # noqa: E402

# Milo parameters
K_NEIGHBORS = 15
NHOOD_PROP = 0.1
SOLVER = "edger"
ALPHA = 0.1
MAX_CELLS = 8000

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = BASE_DIR

GROUP_COL = "stomach_pre_grouping"
GROUP_A = "Responsed"
GROUP_B = "No-response"
SAMPLE_COL = "sample"

PANEL_LETTER = "F"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(2, PANEL_LETTER)
MARGIN = dict(left=6.0, right=10.0, top=4.5, bottom=4.5)
LABEL_PAD_MM = 1.0                  # paper between the axes and its x label

SCALE = 4                           # the earlier canvas multiplier, for MARK
SMALL_PT = 5.0                      # the earlier smallest body type


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    print("=" * 60)
    print("Milo Differential Abundance Analysis")
    print("=" * 60)

    print("\n[1/7] Loading epithelial data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)                 # noqa: F405
    print(f"  Loaded {adata.n_obs} cells")

    print("\n[2/7] Filtering for Pre-treatment cells...")
    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs[GROUP_COL].isin([GROUP_A, GROUP_B])
    adata = adata[pre_mask & valid_mask].copy()
    print(f"  Pre-treatment cells: {adata.n_obs}")

    if adata.n_obs > MAX_CELLS:
        print(f"\n[3/7] Subsampling to {MAX_CELLS} cells...")
        sc.pp.subsample(adata, n_obs=MAX_CELLS, random_state=42)
    else:
        print("\n[3/7] No subsampling needed")

    if 'X_pca' not in adata.obsm:
        print("\n[4/7] Computing PCA...")
        sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    else:
        print("\n[4/7] Using existing PCA")

    print(f"\n[5/7] Computing neighbors (k={K_NEIGHBORS})...")
    sc.pp.neighbors(adata, n_neighbors=K_NEIGHBORS, n_pcs=50, key_added='neighbors')

    if 'X_umap' not in adata.obsm:
        print("  Computing UMAP...")
        sc.tl.umap(adata)

    print(f"\n[6/7] Running Milo analysis...")
    milo = pt.tl.Milo()
    mdata = mu.MuData({"rna": adata})
    milo.make_nhoods(mdata, neighbors_key="neighbors", prop=NHOOD_PROP)
    milo.count_nhoods(mdata, sample_col=SAMPLE_COL)
    milo.da_nhoods(mdata, design=f"~ {GROUP_COL}", solver=SOLVER)
    milo.build_nhood_graph(mdata)

    print("\n[7/7] Creating visualization...")
    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    # Before the neighbourhood graph: it makes the colour bar, and matplotlib
    # sizes a colour bar from the axes box it finds.
    style.margins_mm(fig, **MARGIN)

    milo.plot_nhood_graph(mdata, alpha=1.0, min_logFC=0, ax=ax, title="")

    ax.set_title('Pre-R vs Pre-NR')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    # scanpy anchors the x axis label to the figure rather than to its own
    # axes, so it stays at the foot of the canvas however the axes is placed:
    # on a canvas four times the printed size that is a wide empty margin, and
    # at print size it falls off the bottom edge. It is put back under its own
    # axes, LABEL_PAD_MM of paper below it. The string is unchanged.
    box = ax.get_position()
    ax.xaxis.set_label_coords(0.5, -LABEL_PAD_MM / (PANEL_H_MM * box.height))
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "epithelial_milo_pre_rvsnr")
    print(f"  Saved: epithelial_milo_pre_rvsnr.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")

    print("\nDone!")


if __name__ == '__main__':
    main()
