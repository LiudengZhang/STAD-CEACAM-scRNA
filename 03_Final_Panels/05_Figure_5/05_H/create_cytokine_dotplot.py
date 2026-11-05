#!/usr/bin/env python3
"""Figure 5F: dot plot of TNF, IL6, IL1B and IL1A across thirteen major cell
types (hepatocytes excluded), cmap='Reds', standard_scale='var'.

Expression is recomputed from layers['counts'] as log1p(counts / row sum *
1e4). Per-group values are cached in data/dot_frames.csv; the h5ad is read
only when that table is missing or --recompute is passed.

Input:  FULL_DATASET_DEPOSIT_H5AD (00_Config/paths.py)
Output: cytokine_dotplot.{svg,pdf,png} in this directory
Run:    python create_cytokine_dotplot.py [--recompute]
"""
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FULL_DATASET_DEPOSIT_H5AD
from cnsfig import Panel, style  # noqa: E402
from cnsfig import legend as cnslegend
from cnsfig import cache
import anndata as ad
import numpy as np
import pandas as pd

TARGET_SUM = 1e4

PANEL = Panel(5, "F", __file__, small_pt=10.0)

# Layout in mm: row-name column, size-key width, gene-name band, and the
# band under the matrix that holds the gene names and the keys.
LEFT_MM = 25.0
SIZE_KEY_W_MM = 16.0
GENE_MM = 5.6
BOTTOM_MM = 16.5

GENES = ['TNF', 'IL6', 'IL1B', 'IL1A']

# Cell type order, as in the radar plot (panel H)
CELL_TYPE_ORDER = [
    'B cells',
    'DC cells',
    'Endothelial cells',
    'Epithelial cells',
    'Fibroblast',
    'Mast cells',
    'Monocytes/Macrophages',
    'Neutrophils',
    'NK cells',
    'Pericyte',
    'Plasma cells',
    'CD4+ T cells',
    'CD8+ T cells',
]


def compute_dot_frames():
    """Scaled mean expression and fraction expressing per cell type and gene,
    as a long table."""
    adata = sc.read_h5ad(FULL_DATASET_DEPOSIT_H5AD)

    mask = adata.obs['major_cell_type'].isin(CELL_TYPE_ORDER)
    adata = adata[mask].copy()

    # Normalise from the counts layer (per cell, so unaffected by the filter).
    if 'counts' not in adata.layers:
        raise RuntimeError(
            "layers['counts'] is absent; this panel is drawn from the counts "
            "layer and must not silently fall back to .X")
    adata.X = adata.layers['counts'].copy()
    adata.uns.pop('log1p', None)        # or sc.pp.log1p refuses to record its base
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    del adata.layers['counts']

    var_names = adata.raw.var_names if adata.raw is not None else adata.var_names
    available = [g for g in GENES if g in var_names]

    dp = sc.pl.DotPlot(adata, var_names=available, groupby='major_cell_type',
                       categories_order=CELL_TYPE_ORDER,
                       use_raw=(adata.raw is not None), standard_scale='var')
    color, size = dp.dot_color_df, dp.dot_size_df
    rows = []
    for ct in CELL_TYPE_ORDER:
        for g in available:
            rows.append({"cell_type": ct, "gene": g,
                         "mean_scaled": float(color.loc[ct, g]),
                         "fraction": float(size.loc[ct, g])})
    return pd.DataFrame(rows)


def main():
    PANEL.apply()

    frames = cache.table(PANEL.dir, "dot_frames", compute_dot_frames)
    available = [g for g in GENES if g in set(frames["gene"])]
    color_df = frames.pivot(index="cell_type", columns="gene", values="mean_scaled") \
                     .loc[CELL_TYPE_ORDER, available]
    size_df = frames.pivot(index="cell_type", columns="gene", values="fraction") \
                    .loc[CELL_TYPE_ORDER, available]
    scaffold = ad.AnnData(
        X=np.zeros((len(CELL_TYPE_ORDER), len(available)), dtype=np.float32),
        obs=pd.DataFrame({"major_cell_type": pd.Categorical(
            CELL_TYPE_ORDER, categories=CELL_TYPE_ORDER)},
            index=[f"scaffold_{i}" for i in range(len(CELL_TYPE_ORDER))]),
        var=pd.DataFrame(index=available))

    dotplot = sc.pl.dotplot(
        scaffold,
        var_names=available,
        groupby='major_cell_type',
        categories_order=CELL_TYPE_ORDER,
        use_raw=False,
        cmap='Reds',
        show=False,
        save=None,
        return_fig=True,
        figsize=style.figsize_mm(PANEL.w_mm, PANEL.h_mm),
        dot_color_df=color_df,
        dot_size_df=size_df,
    )
    # Largest dot 2.63 mm across. cmap is passed again because DotPlot.style()
    # resets every style field it takes.
    dotplot.style(largest_dot=(2.63 * style.PT_PER_MM) ** 2, cmap='Reds')
    # No scanpy legend column; the keys are drawn under the matrix.
    dotplot.legend(show=False)
    dotplot.make_figure()
    axes = dotplot.get_axes()

    main_ax = axes['mainplot_ax']
    main_ax.set_xticklabels([t.get_text() for t in main_ax.get_xticklabels()],
                            rotation=90, style='italic')
    main_ax.tick_params(labelsize=style.tick_pt())

    fig = dotplot.fig
    style.margins_mm(fig, left=LEFT_MM, right=0.6, top=4.3, bottom=BOTTOM_MM)
    # scanpy's gridspec keeps its own columns beside the matrix, so place the
    # matrix by hand from the row-name column to the right edge.
    _pos = main_ax.get_position()
    main_ax.set_position([LEFT_MM / PANEL.w_mm, _pos.y0,
                          (PANEL.w_mm - LEFT_MM - 0.6) / PANEL.w_mm, _pos.height])
    for name in ("size_legend_ax", "color_legend_ax"):
        if axes.get(name) is not None:
            axes[name].remove()
    df = dotplot.dot_color_df
    vb = dotplot.vboundnorm
    vmin = float(vb.vmin) if vb.vmin is not None else float(df.min().min())
    vmax = float(vb.vmax) if vb.vmax is not None else float(df.max().max())
    areas, labels = cnslegend.scanpy_dot_areas(dotplot, n=4)
    keep = [(a, t) for a, t in zip(areas, labels) if (a or 0) ** 0.5 >= 1.0]
    areas, labels = [a for a, _ in keep], [t for _, t in keep]
    # Size key at the left, colour bar in the rest of the width.
    x_left = 0.6
    y_top = PANEL.h_mm - BOTTOM_MM + GENE_MM
    W, H = PANEL.w_mm, PANEL.h_mm
    tp = style.tick_pt()

    def band_axes(x_mm, w_mm, y_mm, h_mm):
        return fig.add_axes([x_mm / W, 1 - (y_mm + h_mm) / H, w_mm / W, h_mm / H])

    head_h = (tp * 1.15 + 1.0) / style.PT_PER_MM
    hax = band_axes(x_left, SIZE_KEY_W_MM, y_top, head_h); hax.set_axis_off()
    hax.text(0.5, 1.0, 'Cells in group (%)', ha='center', va='top', fontsize=tp,
             transform=hax.transAxes)
    circ_h = (2 * max(cnslegend._radii_pt(areas)) + 1.6 * tp + 2.0) / style.PT_PER_MM
    sax = band_axes(x_left, SIZE_KEY_W_MM, y_top + head_h + 0.3, circ_h)
    cnslegend.dot_size_key(sax, areas, labels, label_pt=tp, colour='gray',
                           edge_colour='black', edge_lw=style.EDGE_PT,
                           orientation='h')
    cx = x_left + SIZE_KEY_W_MM + 1.5
    cw = W - 0.6 - cx
    hax2 = band_axes(cx, cw, y_top, 2 * head_h); hax2.set_axis_off()
    hax2.text(0.5, 1.0, 'Mean expression\nin group', ha='center', va='top',
              fontsize=tp, linespacing=1.15, transform=hax2.transAxes)
    cax = band_axes(cx, cw, y_top + 2 * head_h + 0.6, 1.5)
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize
    cb = ColorbarBase(cax, cmap=plt.get_cmap(dotplot.cmap), norm=Normalize(vmin, vmax),
                      orientation='horizontal', ticks=[vmin, (vmin + vmax) / 2, vmax])
    cb.set_ticklabels([f"{t:g}" for t in (vmin, (vmin + vmax) / 2, vmax)])
    cax.tick_params(labelsize=tp, width=style.RULE_PT, length=1.5, pad=1.0)
    cb.outline.set_linewidth(style.RULE_PT)
    labs = cax.get_xticklabels()
    labs[0].set_ha('left'); labs[-1].set_ha('right')
    key = {'size_legend_ax': sax, 'color_legend_ax': cax}
    cnslegend.require_size_key(key, dot_areas=areas,
                               panel='Figure 5 panel F')

    # Keep everything vector (scanpy rasterises some elements).
    for one in dotplot.fig.get_axes():
        for coll in one.collections:
            coll.set_rasterized(False)
        for img in one.images:
            img.set_rasterized(False)

    PANEL.finish(fig, 'cytokine_dotplot')


if __name__ == "__main__":
    main()
