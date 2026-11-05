#!/usr/bin/env python3
"""
Figure 2E: Milo differential abundance on the pre-treatment epithelial
neighbourhood graph, responders versus non-responders.

Milo is run on every pre-treatment tumour sample with a response label
(stomach and liver, all cells, design ~ site + response). SpatialFDR is
computed within the neighbourhoods where more than half of the cells are
Epi_CEACAM5/6. Neighbourhoods below ALPHA are drawn in their logFC colour with
a thin outline; the rest as a pale tint of their colour. A dashed contour
marks the Epi_CEACAM5/6 cluster.

Reads data/nhoods.csv and data/ceacam_cells.csv beside this script; with
--recompute, runs Milo on the epithelial h5ad (EPITHELIAL_H5AD) and rewrites
them. Milo needs pertpy, so recompute in the pertpy environment:

    conda run -n pertpy_milo python create_epithelial_milo_pre_rvsnr.py --recompute

Drawing from the cached tables runs in the plotting environment.
Writes epithelial_milo_pre_rvsnr.[svg|pdf|png] beside this script.
"""

import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import os
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD                          # noqa: E402
from cnsfig import Panel, style  # noqa: E402
from cnsfig import cache                                  # noqa: E402
from cnsfig.layout import umap_corner_axis                 # noqa: E402

# Milo parameters
K_NEIGHBORS = 15
NHOOD_PROP = 0.1
SOLVER = "edger"
ALPHA = 0.1
PLOT_ALPHA = 1.0                    # plot_nhood_graph(alpha=1.0): colour every nhood
TINT = 0.3                          # non-significant: 30% of the logFC colour on white
# Outline of a significant neighbourhood; thinner than usual so the many
# outlines in the CEACAM cluster stay distinct.
SIG_EDGE_PT = 0.1

# Samples: every pre-treatment tumour sample with a response label, stomach
# and liver. SpatialFDR within the Epi_CEACAM5/6 neighbourhoods uses pertpy's
# weighted BH on the subset (P values outside it set to NaN); the global
# SpatialFDR is kept in the table as well.
SITES = {"Stomach": "stomach_pre_grouping", "Liver": "liver_pre_grouping"}
GROUP_COL = "response"
GROUP_A = "Responsed"
GROUP_B = "No-response"
SAMPLE_COL = "sample"
DESIGN = f"~ site + {GROUP_COL}"
FAMILY_FRACTION = 0.5

SCALE = 4                           # sets MARK relative to the tick size
SMALL_PT = 5.0

PANEL = Panel(2, "E", __file__, small_pt=SMALL_PT)
# The colour bar sits in the right margin; the graph takes the rest.
MARGIN = dict(left=2.0, right=11.0, top=4.5, bottom=2.5)
CBAR_W_MM = 2.0
CBAR_GAP_MM = 1.0

# Diameter (mm) of a neighbourhood at the 90th percentile of size; all
# circles scale by the same factor, so sizes stay proportional.
NHOOD_P90_MM = 1.0
# The cluster outlined, as its obs value and as printed.
CEACAM_STATE = "C2_Epi_CEACAM6"
OUTLINE_LABEL = "Epi_CEACAM5/6"
OUTLINE_FRACTION = 0.90


_STATE = {}


def _run_milo():
    """Run Milo. pertpy and scanpy are imported only here, so drawing from
    the cached tables does not need them."""
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

    adata = sc.read_h5ad(EPITHELIAL_H5AD)

    obs = adata.obs
    group = pd.Series("Not selected", index=obs.index, dtype=object)
    for site, col in SITES.items():
        on_site = (obs['Sample site'] == site).to_numpy()
        group[on_site] = obs.loc[on_site, col].astype(str)
    pre_mask = (obs['Treatment phase'] == 'Pre').to_numpy()
    valid_mask = group.isin([GROUP_A, GROUP_B]).to_numpy()
    adata = adata[pre_mask & valid_mask].copy()
    adata.obs[GROUP_COL] = pd.Categorical(group[pre_mask & valid_mask],
                                          categories=[GROUP_B, GROUP_A])
    adata.obs['site'] = adata.obs['Sample site'].astype(str)
    # Response is a patient property: a patient in both groups would mean
    # the labels disagree between sites.
    if (adata.obs.groupby('Patient ID', observed=True)[GROUP_COL].nunique() > 1).any():
        raise RuntimeError("a patient carries both response labels")

    if 'X_pca' not in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

    sc.pp.neighbors(adata, n_neighbors=K_NEIGHBORS, n_pcs=50, key_added='neighbors')

    if 'X_umap' not in adata.obsm:
        sc.tl.umap(adata)

    milo = pt.tl.Milo()
    mdata = mu.MuData({"rna": adata})
    milo.make_nhoods(mdata, neighbors_key="neighbors", prop=NHOOD_PROP)
    milo.count_nhoods(mdata, sample_col=SAMPLE_COL)
    milo.da_nhoods(mdata, design=DESIGN, solver=SOLVER)
    milo.build_nhood_graph(mdata)

    # SpatialFDR within the Epi_CEACAM5/6 neighbourhoods.
    rna, var = mdata["rna"], mdata["milo"].var
    nh = rna.obsm["nhoods"].tocsc()
    in_c = (rna.obs["minor_cell_state"].astype(str) == CEACAM_STATE).to_numpy()
    frac = np.array([in_c[nh[:, j].indices].mean() for j in range(nh.shape[1])])
    family = frac > FAMILY_FRACTION
    p_all, fdr_all = var["PValue"].copy(), var["SpatialFDR"].copy()
    var.loc[~family, "PValue"] = np.nan
    milo._graph_spatial_fdr(mdata["milo"])
    var["SpatialFDR_ceacam"] = var["SpatialFDR"].to_numpy()
    var["PValue"], var["SpatialFDR"] = p_all, fdr_all
    var["ceacam_fraction"] = frac
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
        "ceacam_fraction": obs["ceacam_fraction"].astype(float).to_numpy(),
        "SpatialFDR_ceacam": obs["SpatialFDR_ceacam"].astype(float).to_numpy(),
    })
    # Every neighbourhood is coloured by logFC (PLOT_ALPHA = 1.0); ALPHA is
    # the significance threshold, not the colouring's.
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
    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)

    nhoods = cache.table(PANEL.dir, "nhoods", compute_nhoods)
    cells = cache.table(PANEL.dir, "ceacam_cells", compute_ceacam_cells)

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)
    style.margins_mm(fig, **MARGIN)

    # Marker area is Nhood_size * min_size (as in pertpy); min_size is set
    # from NHOOD_P90_MM instead of pertpy's default of 10.
    sizes = nhoods["Nhood_size"].astype(float)
    nhood_p90 = float(np.percentile(sizes, 90))
    min_size = (NHOOD_P90_MM * style.PT_PER_MM) ** 2 / nhood_p90

    # RdBu_r on symmetric limits, largest |logFC| on top. Significant
    # neighbourhoods at full colour with an outline; the rest as an opaque
    # tint (not transparency, so overlaps do not darken).
    vmax = float(np.nanmax(np.abs(nhoods["graph_color"])))
    norm = mpl.colors.Normalize(-vmax, vmax)
    sig = (nhoods["SpatialFDR_ceacam"] < ALPHA).to_numpy()
    by_abs = nhoods.assign(abs_logFC=nhoods["logFC"].abs()).sort_values("abs_logFC")
    bg, fg = by_abs[~sig[by_abs.index]], by_abs[sig[by_abs.index]]
    tint = mpl.colormaps["RdBu_r"](norm(bg["logFC"].to_numpy()))[:, :3]
    ax.scatter(bg["x"], bg["y"], c=TINT * tint + (1 - TINT),
               s=bg["Nhood_size"] * min_size, edgecolors="none", zorder=2)
    colour_layer = ax.scatter(fg["x"], fg["y"], c=fg["logFC"], cmap="RdBu_r", norm=norm,
                              s=fg["Nhood_size"] * min_size, edgecolors="black",
                              linewidths=SIG_EDGE_PT, zorder=3)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    # The colour bar, at millimetres in the right margin.
    box = ax.get_position()
    cax = fig.add_axes([box.x1 + CBAR_GAP_MM / PANEL.w_mm, box.y0,
                        CBAR_W_MM / PANEL.w_mm, box.height])
    fig.colorbar(colour_layer, cax=cax)
    cax.tick_params(length=2.0, width=style.RULE_PT)
    cax.set_frame_on(True)
    for sp in cax.spines.values():
        sp.set_linewidth(style.RULE_PT)

    umap_pts = cells[["x", "y"]].to_numpy()
    # Dashed outline: the density contour enclosing OUTLINE_FRACTION of the
    # Epi_CEACAM5/6 cells, with the cluster name centred below it.
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
    x_mid = (GX[inside].min() + GX[inside].max()) / 2
    y_low = GY[inside].min()
    from matplotlib.transforms import offset_copy
    ax.text(x_mid, y_low, OUTLINE_LABEL, ha="center", va="top",
            fontsize=style.tick_pt(), zorder=5,
            transform=offset_copy(ax.transData, fig=fig, y=-0.8 / 25.4,
                                  units="inches"))

    ax.set_title('Pre-R vs Pre-NR')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    # UMAP1/UMAP2 arrows in the lower-left margin, shorter than on the other
    # maps because the graph fills the axes.
    umap_corner_axis(ax, inset_mm=(1.1, 0.6), length_mm=4.5)

    os.makedirs(PANEL.dir, exist_ok=True)
    PANEL.finish(fig, 'epithelial_milo_pre_rvsnr')


if __name__ == '__main__':
    main()
