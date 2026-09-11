#!/usr/bin/env python3
"""
Figure 4 panel E - canonical markers of the seven MoMac states, three genes per
state, twenty-one genes in the printed order.

  printed panel  Figure 4 E       (PROVENANCE.csv - the directory is "04_E",
                                   which also holds the drawing for panel D;
                                   do NOT read the directory as the letter)

The marker list, the state order, the short state names, the FCGR3A/CD16
display name, `use_raw`, `cmap='Reds'`, `standard_scale='var'` and the marker
summary table are unchanged. `sc.pl.dotplot` still computes and draws
everything. Only the canvas and the type change.

ONE STATE NAME AND ONE GENE SYMBOL ARE RESTORED
    The page names the C3 state MoMac_IL1B and names the gene by its HGNC
    symbol, FCGR3A. The script carried the earlier Mac_IL1B and displayed
    FCGR3A as CD16. The same row and the same gene are plotted, in the same
    order, from the same values; only the printed names change.

MARK
    The earlier drawing set no font size of its own, so every string it drew
    came out at matplotlib's default 10 pt on a canvas scanpy sized itself.
    That reads as SCALE = 1 and SMALL_PT = 10, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    and it is unused: the script specifies no marker area, line width or cap
    size of its own - scanpy sizes the dots from the axes it is given. The
    factor is stated anyway so the panel is auditable like the rest.

Two deviations from the usual recipe, forced by scanpy drawing the figure:

  1. `sc.pl.dotplot` builds its own figure, so the panel box is passed to it as
     `figsize=` rather than through `style.subplots_mm`. scanpy's `figsize` is
     the whole figure: asking for w x h mm gives a w x h mm canvas.
  2. scanpy lays the dotplot out flush against the left edge and lets the y
     tick labels hang outside the canvas. `style.fit_margins` moves scanpy's
     gridspec inwards until the ink is inside the canvas and the panel-letter
     corner is clear, and raises if it cannot. `tight_layout` did this job when
     the save cropped to the ink, and a 1:1 save does not crop.

The colour bar's heading is the only string on the panel set over more than one
line. It is given a line spacing of 1.15 rather than matplotlib's 1.2, which is
what the published page sets it at and what keeps the heading clear of its own
scale.

The dot-size key scanpy draws beside the colour bar is suppressed, because the
published page does not carry one.

Input : MOMAC_H5AD (00_Config/paths.py)
Output: this directory / momac_marker_dotplot.{svg,pdf,png}
        this directory / momac_marker_summary.csv
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from collections import OrderedDict
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import MOMAC_H5AD                # noqa: E402
import panel_style_cns as style             # noqa: E402
import slots                                # noqa: E402

# Output directory
OUTPUT_DIR = Path(__file__).resolve().parent

SCALE = 1                       # the earlier drawing had no SCALE: 1:1 at 10 pt
SMALL_PT = 10.0                 # matplotlib's default font.size, left unchanged

PANEL_LETTER = "E"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(4, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(4, PANEL_LETTER)

#: Matplotlib's own default is 1.2; the published headings are set tighter.
LEGEND_LINESPACING = 1.15

# Scanpy settings
sc.settings.verbosity = 1

# ============================================================================
# CANONICAL MARKERS BY MOMAC SUBSET
# ============================================================================

# Short display names for dotplot y-axis
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

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    print("=" * 80)
    print("MoMac Subset Marker Dotplot Generator")
    print("=" * 80)

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f} (unused here)")

    # 1. Load Data
    print("\n[1/4] Loading data...")
    print(f"  Path: {MOMAC_H5AD}")

    adata = sc.read_h5ad(MOMAC_H5AD)
    print(f"  Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    if 'minor_cell_state' not in adata.obs.columns:
        raise ValueError("minor_cell_state column not found in adata.obs")

    # Rename to short display names
    adata.obs['momac_short'] = adata.obs['minor_cell_state'].replace(SHORT_NAMES)

    # 2. Validate Markers
    print("\n[2/4] Validating marker genes...")
    all_markers = []
    for cell_type, markers in CANONICAL_MARKERS.items():
        all_markers.extend(markers)

    if adata.raw is not None:
        var_names = adata.raw.var_names
    else:
        var_names = adata.var_names

    available_markers = [m for m in all_markers if m in var_names]
    missing_markers = [m for m in all_markers if m not in var_names]

    print(f"  Available: {len(available_markers)}/{len(all_markers)} markers")
    if missing_markers:
        print(f"  Missing markers: {missing_markers}")

    filtered_markers = OrderedDict()
    for cell_type, markers in CANONICAL_MARKERS.items():
        avail = [m for m in markers if m in var_names]
        if avail:
            filtered_markers[cell_type] = avail

    # 3. Prepare
    print("\n[3/4] Preparing visualization...")
    cell_type_order = list(CANONICAL_MARKERS.keys())

    plot_genes = []
    for cell_type in cell_type_order:
        if cell_type in filtered_markers:
            plot_genes.extend(filtered_markers[cell_type])

    # 4. Create Dotplot
    print("\n[4/4] Creating dotplot...")
    dotplot = sc.pl.dotplot(
        adata,
        var_names=plot_genes,
        groupby='momac_short',
        categories_order=cell_type_order,
        use_raw=(adata.raw is not None),
        cmap='Reds',
        show=False,
        save=None,
        standard_scale='var',
        return_fig=True,
        figsize=style.figsize_mm(PANEL_W_MM, PANEL_H_MM),
    )
    # The published page carries the colour bar and no dot-size key, so the key
    # scanpy draws by default is suppressed: a legend the paper never printed
    # would be new content, and this is a visualisation-only change.
    dotplot.legend(show_size_legend=False)
    dotplot.make_figure()
    axes = dotplot.get_axes()

    # scanpy asks matplotlib for the relative size 'small', which is 0.833 of
    # the base size. The sizes are set explicitly instead, from the same two
    # the rest of the figure uses.
    for one in axes.values():
        one.tick_params(labelsize=style.tick_pt())
        if one.get_title():
            one.set_title(one.get_title(), fontsize=style.body_pt(),
                          linespacing=LEGEND_LINESPACING)

    # De-rasterize all axes (scanpy dotplot may rasterize internally)
    for ax in dotplot.fig.get_axes():
        for coll in ax.collections:
            coll.set_rasterized(False)
        for img in ax.images:
            img.set_rasterized(False)

    fig = dotplot.fig
    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    stem = 'momac_marker_dotplot'
    style.save_panel(fig, OUTPUT_DIR / stem)
    print(f"  Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")

    # Summary table
    summary_data = []
    for cell_type in cell_type_order:
        if cell_type in filtered_markers:
            markers = filtered_markers[cell_type]
            summary_data.append({
                'Cell_State': cell_type,
                'Markers': ', '.join(markers),
                'N_Markers': len(markers)
            })

    summary_df = pd.DataFrame(summary_data)
    summary_path = OUTPUT_DIR / 'momac_marker_summary.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f"  Saved: {summary_path}")
    print("\nComplete!")


if __name__ == "__main__":
    main()
