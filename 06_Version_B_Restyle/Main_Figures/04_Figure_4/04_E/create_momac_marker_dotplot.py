#!/usr/bin/env python3
"""
Figure 4 panel E, RESTYLED (Version B) - canonical marker dotplot for the seven
MoMac minor cell states, three markers each, 21 genes.

Note the letter. This script lives in the directory `04_E` and PROVENANCE.csv
does record it as printed panel **E** - but the same directory also holds
`create_momac_umap.py`, which is printed panel **D**, and there is no directory
`04_D`. The agreement here is a coincidence, not a rule; CLAUDE.md rule 2 still
applies.

Version A is
`03_Revised_Panels/Main_Figures/04_Figure_4/04_E/create_momac_marker_dotplot.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, margins are millimetres, and the save is
`style.save_panel`.

Every gene, every cell state, every ordering, `standard_scale='var'`, the raw
layer it reads from and every string is Version A's. The drawing code is the
same code.

  printed panel  Figure 4 E     (PROVENANCE.csv)
  printed rect   114.9 x 39.2 mm   (panel_rects.csv)
  Version B box  160.0 x 50.0 mm

MARK
    Version A sets no font size of its own. Its type came from
    `sc.settings.set_figure_params()`, which puts `font.size = 14` into
    rcParams, and scanpy's DotPlot then labels its rows and columns at
    matplotlib's relative `'small'`, i.e. 0.833 x 14 = 11.66 pt. The assembler
    fitted the saved 618.2 x 262.0 pt SVG into 114.9 x 39.2 mm, a fit of
    0.4240, so that type printed at 4.94 pt. There is no SCALE constant, so
    SCALE = 1, SMALL_PT = 11.66 and

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 11.66 = 0.600

    MARK is recorded for the audit but is not applied anywhere: this panel sets
    no marker size, no line width and no pad in points. Every length in it is
    scanpy's, and scanpy scales them with the canvas.

THE CANVAS
    `sc.pl.dotplot` builds its own figure; Version A's
    `plt.figure(figsize=(fig_width, fig_height))` was discarded by scanpy and
    never reached the page, which is why those two lines are gone rather than
    converted. The 1:1 box is passed through scanpy's own `figsize` argument -
    the supported way to size a DotPlot - and the axes region inside it is then
    placed with `style.margins_mm`. scanpy builds the layout with a plain
    `GridSpec`, so it honours the figure's subplot parameters and the
    millimetre margins take effect.

THE 7 pt FLOOR
    scanpy labels the dotplot with the relative size `'small'`. Against
    cnsplots' `font.size = 8` that resolves to 6.66 pt, not 7 - the one place in
    this figure where the type system does not reach the floor the restyle
    exists to establish. Every text artist below `style.tick_pt()` is raised to
    it after the plot is built. Only sizes move; not one string, position or
    value does.

Input : MOMAC_H5AD (00_Config/paths.py)
Output: this directory / momac_marker_dotplot.{svg,pdf,png}
        this directory / momac_marker_summary.csv
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import OrderedDict
import sys

# Central config
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import MOMAC_H5AD                # noqa: E402
import panel_style_cns as style             # noqa: E402

# Output directory
OUTPUT_DIR = Path(__file__).resolve().parent

# Version A's canvas convention, kept only so MARK can be derived from it.
SCALE = 1
SMALL_PT = 11.66            # scanpy's 'small' against its own font.size = 14

PRINTED_MM = (114.9, 39.2)  # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 160.0, 50.0
MARGIN = dict(left=16.0, right=2.5, top=2.0, bottom=12.2)

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=300, facecolor='white')

# ============================================================================
# CANONICAL MARKERS BY MOMAC SUBSET
# ============================================================================

# Short display names for dotplot y-axis
SHORT_NAMES = OrderedDict([
    ('C0_Mac_Classic_TREM2',           'Mac_TREM2'),
    ('C1_Mono_Classic_CD14',           'Mono_CD14'),
    ('C2_MoMac_Intermediate_HLA-DRA',  'MoMac_Inter'),
    ('C3_Mac_Inflam_IL1B',             'Mac_IL1B'),
    ('C4_Mono_Alternative_CD16',       'Mono_CD16'),
    ('C5_Mac_Prolif_MKI67',            'Mac_Prolif'),
    ('C6_Mac_Metallothionein_MT1G',    'Mac_MT1G'),
])

CANONICAL_MARKERS = OrderedDict([
    ('Mac_TREM2',    ['TREM2', 'C1QA', 'APOE']),
    ('Mono_CD14',    ['CD14', 'S100A8', 'VCAN']),
    ('MoMac_Inter',  ['HLA-DRA', 'CST3', 'CLEC10A']),
    ('Mac_IL1B',     ['IL1B', 'TNF', 'CXCL8']),
    ('Mono_CD16',    ['CD16', 'CDKN1C', 'LST1']),
    ('Mac_Prolif',   ['MKI67', 'TOP2A', 'STMN1']),
    ('Mac_MT1G',     ['MT1G', 'MT2A', 'MT1X']),
])


def _text_artists(fig):
    """Every text artist scanpy put on the figure, wherever it lives."""
    out = list(fig.texts)
    for ax in fig.axes:
        out += [ax.title, ax.xaxis.label, ax.yaxis.label] + list(ax.texts)
        out += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
    return [t for t in out if t.get_text()]


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    print("=" * 80)
    print("MoMac Subset Marker Dotplot Generator")
    print("=" * 80)

    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f} (unused - this "
          f"panel sets no lengths in points)")

    # 1. Load Data
    print("\n[1/4] Loading data...")
    print(f"  Path: {MOMAC_H5AD}")

    adata = sc.read_h5ad(MOMAC_H5AD)
    print(f"  Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    if 'minor_cell_state' not in adata.obs.columns:
        raise ValueError("minor_cell_state column not found in adata.obs")

    # Rename to short display names
    adata.obs['momac_short'] = adata.obs['minor_cell_state'].replace(SHORT_NAMES)

    # Rename genes for display (FCGR3A → CD16)
    GENE_DISPLAY = {'FCGR3A': 'CD16'}
    adata.var_names = pd.Index([GENE_DISPLAY.get(g, g) for g in adata.var_names])
    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        raw_adata.var_names = pd.Index([GENE_DISPLAY.get(g, g) for g in raw_adata.var_names])
        adata.raw = raw_adata

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

    sc.pl.dotplot(
        adata,
        var_names=plot_genes,
        groupby='momac_short',
        categories_order=cell_type_order,
        use_raw=(adata.raw is not None),
        cmap='Reds',
        show=False,
        save=None,
        standard_scale='var',
        figsize=style.figsize_mm(PANEL_W_MM, PANEL_H_MM),
    )

    fig = plt.gcf()

    # De-rasterize all axes (scanpy dotplot may rasterize internally)
    for ax in fig.get_axes():
        for coll in ax.collections:
            coll.set_rasterized(False)
        for img in ax.images:
            img.set_rasterized(False)

    # The 7 pt floor - see the module docstring. scanpy labels the plot at the
    # relative size 'small', which against cnsplots' 8 pt base is 6.66 pt.
    for t in _text_artists(fig):
        if t.get_fontsize() < style.tick_pt():
            t.set_fontsize(style.tick_pt())

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUTPUT_DIR / 'momac_marker_dotplot')
    print(f"  Saved: {OUTPUT_DIR / 'momac_marker_dotplot'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")

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
