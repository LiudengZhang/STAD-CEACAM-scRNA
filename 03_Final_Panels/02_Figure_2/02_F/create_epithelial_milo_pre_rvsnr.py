#!/usr/bin/env python3
"""
Figure 2 panel F - Milo differential abundance over the pre-treatment
epithelial neighbourhood graph, responders versus non-responders.

COMPUTE IN THE pertpy ENVIRONMENT; DRAW ANYWHERE (since 2026-09-15).

The Milo run - subsampling, PCA, neighbours, UMAP, neighbourhoods, edgeR -
is the computing half and needs pertpy, whose scvi/jax imports do not resolve
in stad_ceacam. Its result is two tables under data/ (cnsfig.cache):

    nhoods.csv        one neighbourhood per row: its graph position, size,
                      logFC, SpatialFDR and the colour value pertpy's
                      plot_nhood_graph would give it (NaN above alpha)
    ceacam_cells.csv  the UMAP position of every Epi_CEACAM5/6 cell, for the
                      dashed outline

The drawing half reads those two tables and nothing else, so it runs in the
cnsplots environment like every other panel, in seconds. To recompute:

    conda run -n pertpy_milo python create_epithelial_milo_pre_rvsnr.py --recompute

The drawing is pertpy's plot_nhood_graph re-expressed in matplotlib from the
table (RdBu_r, symmetric limits at the largest |logFC|, extreme |logFC| drawn
on top, non-significant neighbourhoods light grey, area = Nhood_size * min_size)
- proven against the pertpy drawing's harvest with compare_panel_content.py
--against-json on 2026-09-15.

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

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
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
from cnsfig import cache                                  # noqa: E402

# Milo parameters
K_NEIGHBORS = 15
NHOOD_PROP = 0.1
SOLVER = "edger"
ALPHA = 0.1
PLOT_ALPHA = 1.0                    # plot_nhood_graph(alpha=1.0): colour every nhood
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
#: The embedding fills the slot since 2026-09-15 (the author's third reading:
#: "the circles could be larger still"): the axes had left 40% of the slot's
#: width to margins and a scanpy colour bar. The colour bar is drawn here at
#: CBAR_W_MM in the right margin; the axes takes the rest.
MARGIN = dict(left=2.0, right=11.0, top=4.5, bottom=2.5)
CBAR_W_MM = 2.0
CBAR_GAP_MM = 1.0

SCALE = 4                           # the earlier canvas multiplier, for MARK
SMALL_PT = 5.0                      # the earlier smallest body type

#: The neighbourhood circle at the 90th percentile of neighbourhood size,
#: measured off the published raster at 600 dpi: the isolated components run
#: 0.5-1.6 mm across and their 90th percentile is 1.55 mm. Pinned at the
#: percentile rather than the maximum because one outsized neighbourhood
#: would otherwise shrink every other circle.
#: 2026-09-14 (evening), the author's ruling: the circles read small at the
#: page's size, so every circle is scaled up by the same factor (1.55 -> 2.3
#: mm at the 90th percentile). Proportional: the size ordering and the ratio
#: between any two neighbourhoods are unchanged. 2026-09-15, third reading:
#: "a little larger still" - 2.3 -> 3.0 mm, again proportional.
NHOOD_P90_MM = 3.0
#: The cluster the published page rings, as its obs value; 2A prints it as
#: CEACAM5/6 and the page as Epi_CEACAM5/6.
CEACAM_STATE = "C2_Epi_CEACAM6"
OUTLINE_LABEL = "Epi_CEACAM5/6"
OUTLINE_FRACTION = 0.90


_STATE = {}


def _run_milo():
    """The computing half. pertpy and scanpy are imported here and nowhere
    else, so the drawing half never needs them."""
    if _STATE:
        return _STATE["mdata"]
    try:
        import scanpy as sc
        import pertpy as pt
        import mudata as mu
    except ImportError as exc:
        raise SystemExit(
            f"data/ is missing a table and this environment cannot compute it "
            f"({exc}). Run\n  conda run -n pertpy_milo python "
            f"{Path(__file__).name} --recompute") from exc
    print("=" * 60)
    print("Milo Differential Abundance Analysis")
    print("=" * 60)

    print("\n[1/6] Loading epithelial data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)                 # noqa: F405
    print(f"  Loaded {adata.n_obs} cells")

    print("\n[2/6] Filtering for Pre-treatment cells...")
    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs[GROUP_COL].isin([GROUP_A, GROUP_B])
    adata = adata[pre_mask & valid_mask].copy()
    print(f"  Pre-treatment cells: {adata.n_obs}")

    if adata.n_obs > MAX_CELLS:
        print(f"\n[3/6] Subsampling to {MAX_CELLS} cells...")
        sc.pp.subsample(adata, n_obs=MAX_CELLS, random_state=42)
    else:
        print("\n[3/6] No subsampling needed")

    if 'X_pca' not in adata.obsm:
        print("\n[4/6] Computing PCA...")
        sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    else:
        print("\n[4/6] Using existing PCA")

    print(f"\n[5/6] Computing neighbors (k={K_NEIGHBORS})...")
    sc.pp.neighbors(adata, n_neighbors=K_NEIGHBORS, n_pcs=50, key_added='neighbors')

    if 'X_umap' not in adata.obsm:
        print("  Computing UMAP...")
        sc.tl.umap(adata)

    print(f"\n[6/6] Running Milo analysis...")
    milo = pt.tl.Milo()
    mdata = mu.MuData({"rna": adata})
    milo.make_nhoods(mdata, neighbors_key="neighbors", prop=NHOOD_PROP)
    milo.count_nhoods(mdata, sample_col=SAMPLE_COL)
    milo.da_nhoods(mdata, design=f"~ {GROUP_COL}", solver=SOLVER)
    milo.build_nhood_graph(mdata)
    _STATE["mdata"] = mdata
    return mdata


def compute_nhoods():
    """One row per neighbourhood, with the colour value plot_nhood_graph
    assigns: logFC where SpatialFDR <= PLOT_ALPHA (1.0: all), else NaN."""
    mdata = _run_milo()
    nh = mdata["milo"].T
    obs = nh.obs
    xy = nh.obsm["X_milo_graph"]
    df = pd.DataFrame({
        "nhood": np.arange(nh.n_obs),
        "x": xy[:, 0], "y": xy[:, 1],
        "Nhood_size": obs["Nhood_size"].astype(float).to_numpy(),
        "logFC": obs["logFC"].astype(float).to_numpy(),
        "SpatialFDR": obs["SpatialFDR"].astype(float).to_numpy(),
    })
    # plot_nhood_graph(alpha=1.0, min_logFC=0): every neighbourhood coloured
    # by its logFC, as the drawing has always called it; ALPHA (0.1) is the
    # test's threshold, not the colouring's.
    color = df["logFC"].copy()
    color[df["SpatialFDR"] > PLOT_ALPHA] = np.nan
    df["graph_color"] = color
    return df


def compute_ceacam_cells():
    mdata = _run_milo()
    rna = mdata["rna"]
    in_cluster = (rna.obs["minor_cell_state"].astype(str)
                  == CEACAM_STATE).to_numpy()
    umap = rna.obsm["X_umap"][in_cluster]
    return pd.DataFrame({"x": umap[:, 0], "y": umap[:, 1]})


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    nhoods = cache.table(BASE_DIR, "nhoods", compute_nhoods)
    cells = cache.table(BASE_DIR, "ceacam_cells", compute_ceacam_cells)

    print("\nCreating visualization...")
    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    style.margins_mm(fig, **MARGIN)

    # THE CIRCLES ARE THE PUBLISHED PAGE'S SIZE  (2026-09-14)
    #   pertpy sizes a neighbourhood at Nhood_size * min_size points squared
    #   and defaults min_size to 10, whatever the canvas. On this 1:1 canvas
    #   that drew the largest neighbourhoods 5 mm across and the graph fused
    #   into one opaque blob. The published panel is a raster; measured at
    #   600 dpi its isolated circles run 0.5-1.6 mm across, so the largest
    #   neighbourhood is drawn NHOOD_MAX_MM across and every other one scales
    #   with it. Nothing about which neighbourhoods exist or what colour they
    #   take moves; only the marker area.
    sizes = nhoods["Nhood_size"].astype(float)
    nhood_p90 = float(np.percentile(sizes, 90))
    min_size = (NHOOD_P90_MM * style.PT_PER_MM) ** 2 / nhood_p90
    print(f"  neighbourhoods {sizes.min():.0f}-{sizes.max():.0f} cells, p90 "
          f"{nhood_p90:.0f} -> {NHOOD_P90_MM} mm across (min_size "
          f"{min_size:.4f}); largest {np.sqrt(sizes.max() * min_size) / style.PT_PER_MM:.2f} mm")

    # plot_nhood_graph, from the table: extreme |logFC| on top (NaN first),
    # RdBu_r on symmetric limits, NaN in scanpy's light grey.
    order = nhoods.assign(abs_logFC=nhoods["graph_color"].abs()) \
                  .sort_values("abs_logFC", na_position="first").index
    nh = nhoods.loc[order]
    vmax = float(np.nanmax(np.abs(nh["graph_color"])))
    ax.scatter(nh["x"], nh["y"], c=nh["graph_color"], cmap="RdBu_r",
               vmin=-vmax, vmax=vmax, s=nh["Nhood_size"] * min_size,
               plotnonfinite=True, edgecolors="none", zorder=2)
    ax.collections[-1].cmap.set_bad("lightgray")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    # The colour bar, at millimetres in the right margin.
    box = ax.get_position()
    cax = fig.add_axes([box.x1 + CBAR_GAP_MM / PANEL_W_MM, box.y0,
                        CBAR_W_MM / PANEL_W_MM, box.height])
    fig.colorbar(ax.collections[-1], cax=cax)
    cax.tick_params(length=2.0, width=style.RULE_PT)
    cax.set_frame_on(True)
    for sp in cax.spines.values():
        sp.set_linewidth(style.RULE_PT)

    umap_pts = cells[["x", "y"]].to_numpy()
    # THE OUTLINE AND ITS LABEL ARE THE PUBLISHED PAGE'S  (2026-09-14)
    #   The submitted panel rings the Epi_CEACAM5/6 cluster with a dashed
    #   outline and names it beside the ring. Both were placed at assembly and
    #   no script drew them, which is why the redraw had neither. The ring is
    #   drawn here as the density contour enclosing OUTLINE_FRACTION of that
    #   cluster's cells in the same embedding the neighbourhoods sit on, so it
    #   follows the cluster wherever the embedding puts it; the label sits at
    #   the ring's lower-right, where the page prints it. Declared in
    #   00_Config/shared/labels.py and PROVENANCE.csv.
    if len(umap_pts) < 10:
        raise RuntimeError(f"{CEACAM_STATE!r}: {len(umap_pts)} cells; the "
                           f"outline cannot be drawn")
    pts = umap_pts
    kde = gaussian_kde(pts.T)
    dens = kde(pts.T)
    level = float(np.quantile(dens, 1.0 - OUTLINE_FRACTION))
    pad = 1.0
    gx = np.linspace(pts[:, 0].min() - pad, pts[:, 0].max() + pad, 200)
    gy = np.linspace(pts[:, 1].min() - pad, pts[:, 1].max() + pad, 200)
    GX, GY = np.meshgrid(gx, gy)
    Z = kde(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
    ax.contour(GX, GY, Z, levels=[level], colors="black", linestyles="--",
               linewidths=style.RULE_PT, zorder=4)
    inside = Z >= level
    x_right = GX[inside].max()
    y_low = GY[inside].min()
    ax.text(x_right, y_low, OUTLINE_LABEL, ha="right", va="top",
            fontsize=style.tick_pt(), zorder=5)

    ax.set_title('Pre-R vs Pre-NR')
    # The axis labels UMAP1 / UMAP2 are set as the pertpy drawing set them
    # and, as there, hidden with the frame (frameon=False): the page prints
    # neither. compare_panel_content.py sees the strings either way.
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)
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
