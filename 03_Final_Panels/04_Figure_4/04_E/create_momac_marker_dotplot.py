#!/usr/bin/env python3
"""Figure 4E: canonical markers of the seven MoMac states, three genes per
state, as a scanpy dot plot (cmap='Reds', standard_scale='var').

The per-state means and fractions are cached in data/dot_frames.csv; the h5ad
is read only when that table is missing or --recompute is passed. The plot is
then built from the cached frames through DotPlot's dot_color_df/dot_size_df.

Input:  MOMAC_H5AD (00_Config/paths.py)
Output: momac_marker_dotplot.{svg,pdf,png} and momac_marker_summary.csv
        in this directory
Run:    python create_momac_marker_dotplot.py [--recompute]
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
from collections import OrderedDict
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD                # noqa: E402
from cnsfig import Panel, style  # noqa: E402
from cnsfig import legend as cnslegend      # noqa: E402
from cnsfig import cache                    # noqa: E402
import anndata as ad                        # noqa: E402

PANEL = Panel(4, "E", __file__, small_pt=10.0)

#: Width of the key column at the right; the matrix keeps about 85 mm.
KEY_COLUMN_MM = 14.0

sc.settings.verbosity = 1

SHORT_NAMES = OrderedDict([
    ('C0_Mac_Classic_TREM2',           'Mac_TREM2'),
    ('C1_Mono_Classic_CD14',           'Mono_CD14'),
    ('C2_MoMac_Intermediate_HLA-DRA',  'MoMac_Inter'),
    ('C3_Mac_Inflam_IL1B',             'MoMac_IL1B'),
    ('C4_Mono_Alternative_CD16',       'Mono_CD16'),
    ('C5_Mac_Prolif_MKI67',            'Mac_Prolif'),
    ('C6_Mac_Metallothionein_MT1G',    'Mac_MT1G'),
])

CANONICAL_MARKERS = OrderedDict([
    ('Mac_TREM2',    ['TREM2', 'C1QA', 'APOE']),
    ('Mono_CD14',    ['CD14', 'S100A8', 'VCAN']),
    ('MoMac_Inter',  ['HLA-DRA', 'CST3', 'CLEC10A']),
    ('MoMac_IL1B',   ['IL1B', 'TNF', 'CXCL8']),
    ('Mono_CD16',    ['FCGR3A', 'CDKN1C', 'LST1']),
    ('Mac_Prolif',   ['MKI67', 'TOP2A', 'STMN1']),
    ('Mac_MT1G',     ['MT1G', 'MT2A', 'MT1X']),
])


def compute_dot_frames():
    """Scaled mean expression and fraction expressing for the markers, as a
    long table (cell_state, gene, mean_scaled, fraction)."""
    adata = sc.read_h5ad(MOMAC_H5AD)

    if 'minor_cell_state' not in adata.obs.columns:
        raise ValueError("minor_cell_state column not found in adata.obs")

    adata.obs['momac_short'] = adata.obs['minor_cell_state'].replace(SHORT_NAMES)

    all_markers = []
    for cell_type, markers in CANONICAL_MARKERS.items():
        all_markers.extend(markers)

    if adata.raw is not None:
        var_names = adata.raw.var_names
    else:
        var_names = adata.var_names

    cell_type_order = list(CANONICAL_MARKERS.keys())
    plot_genes = [m for ct in cell_type_order for m in CANONICAL_MARKERS[ct]
                  if m in var_names]

    dp = sc.pl.DotPlot(adata, var_names=plot_genes, groupby='momac_short',
                       categories_order=cell_type_order,
                       use_raw=(adata.raw is not None), standard_scale='var')
    color, size = dp.dot_color_df, dp.dot_size_df
    rows = []
    for ct in cell_type_order:
        for g in plot_genes:
            rows.append({"cell_state": ct, "gene": g,
                         "mean_scaled": float(color.loc[ct, g]),
                         "fraction": float(size.loc[ct, g])})
    return pd.DataFrame(rows)


def main():
    PANEL.apply()

    frames = cache.table(PANEL.dir, "dot_frames", compute_dot_frames)
    cell_type_order = list(CANONICAL_MARKERS.keys())
    present = set(frames["gene"])
    filtered_markers = OrderedDict()
    for cell_type, markers in CANONICAL_MARKERS.items():
        avail = [m for m in markers if m in present]
        if avail:
            filtered_markers[cell_type] = avail
    plot_genes = [g for ct in cell_type_order for g in filtered_markers.get(ct, [])]
    color_df = frames.pivot(index="cell_state", columns="gene", values="mean_scaled") \
                     .loc[cell_type_order, plot_genes]
    size_df = frames.pivot(index="cell_state", columns="gene", values="fraction") \
                    .loc[cell_type_order, plot_genes]
    scaffold = ad.AnnData(
        X=np.zeros((len(cell_type_order), len(plot_genes)), dtype=np.float32),
        obs=pd.DataFrame({"momac_short": pd.Categorical(
            cell_type_order, categories=cell_type_order)},
            index=[f"scaffold_{i}" for i in range(len(cell_type_order))]),
        var=pd.DataFrame(index=plot_genes))

    dotplot = sc.pl.dotplot(
        scaffold,
        var_names=plot_genes,
        groupby='momac_short',
        categories_order=cell_type_order,
        use_raw=False,
        cmap='Reds',
        show=False,
        save=None,
        return_fig=True,
        figsize=style.figsize_mm(PANEL.w_mm, PANEL.h_mm),
        dot_color_df=color_df,
        dot_size_df=size_df,
    )
    # Largest dot 2.58 mm across. cmap is passed again because DotPlot.style()
    # resets every style field it takes.
    dotplot.style(largest_dot=(2.58 * style.PT_PER_MM) ** 2, cmap='Reds')
    # A narrow key column; cnsfig.legend replaces scanpy's key axes with two
    # stacked keys at the top of it.
    dotplot.legend(width=KEY_COLUMN_MM / 25.4)
    dotplot.make_figure()
    axes = dotplot.get_axes()
    main_ax = axes['mainplot_ax']
    # Gene symbols italic, at 55 degrees so neighbours do not overlap.
    main_ax.set_xticklabels([t.get_text() for t in main_ax.get_xticklabels()],
                            rotation=55, ha='right', style='italic',
                            fontsize=style.tick_pt())
    main_ax.tick_params(labelsize=style.tick_pt())

    # De-rasterize all axes (scanpy dotplot may rasterize internally)
    for ax in dotplot.fig.get_axes():
        for coll in ax.collections:
            coll.set_rasterized(False)
        for img in ax.images:
            img.set_rasterized(False)

    fig = dotplot.fig
    style.margins_mm(fig, left=14.5, right=KEY_COLUMN_MM + 1.5, top=4.0,
                     bottom=10.5)
    # scanpy's gridspec keeps its own legend column, so place the main axes
    # by hand, ending 2.5 mm short of the key.
    _pos = main_ax.get_position()
    _W = PANEL.w_mm
    main_ax.set_position([14.5 / _W, _pos.y0,
                          (_W - 14.5 - KEY_COLUMN_MM - 2.5) / _W, _pos.height])
    key, areas = cnslegend.scanpy_compact_key(
        dotplot, axes, KEY_COLUMN_MM,
        size_title='Cells in\ngroup (%)',
        cbar_title='Mean\nexpression\nin group',
        label_pt=style.tick_pt(), title_pt=style.tick_pt())
    # Raises if the dot area varies and nothing decodes it.
    cnslegend.require_size_key(key, dot_areas=areas, panel='Figure 4 panel E')

    PANEL.finish(fig, 'momac_marker_dotplot')

    summary_data = []
    for cell_type in cell_type_order:
        if cell_type in filtered_markers:
            markers = filtered_markers[cell_type]
            summary_data.append({
                'Cell_State': cell_type,
                'Markers': ', '.join(markers),
                'N_Markers': len(markers)
            })
    pd.DataFrame(summary_data).to_csv(
        PANEL.dir / 'momac_marker_summary.csv', index=False)


if __name__ == "__main__":
    main()
