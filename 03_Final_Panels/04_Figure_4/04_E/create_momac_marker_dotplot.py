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
from cnsfig import legend as cnslegend      # noqa: E402
from cnsfig import cache                    # noqa: E402
import anndata as ad                        # noqa: E402

# Output directory
OUTPUT_DIR = Path(__file__).resolve().parent

SCALE = 1                       # the earlier drawing had no SCALE: 1:1 at 10 pt
SMALL_PT = 10.0                 # matplotlib's default font.size, left unchanged

PANEL_LETTER = "E"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(4, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(4, PANEL_LETTER)

#: Matplotlib's own default is 1.2; the published headings are set tighter.
LEGEND_LINESPACING = 1.15
#: The size-and-colour key column at the right of the matrix (2026-09-14).
KEY_COLUMN_MM = 14.0    # the matrix keeps 85 mm, a 4.0 mm gene pitch, so 45 degrees clears

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

def compute_dot_frames():
    """The computing half: scanpy's own standard_scale='var' means and
    fraction-expressing matrices for the canonical markers, as a long table
    (cell_state, gene, mean_scaled, fraction)."""
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

    cell_type_order = list(CANONICAL_MARKERS.keys())
    plot_genes = [m for ct in cell_type_order for m in CANONICAL_MARKERS[ct]
                  if m in var_names]

    print("\n[3/4] Summarising...")
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
    print("=" * 80)
    print("MoMac Subset Marker Dotplot Generator")
    print("=" * 80)

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f} (unused here)")

    # THE MATRICES ARE READ FROM data/, since 2026-09-15 (cnsfig.cache).
    #   Opening the MoMac h5ad and summarising means and fractions is the
    #   computing half (compute_dot_frames); it runs only when
    #   data/dot_frames.csv is absent or --recompute is passed. The drawing
    #   builds scanpy's DotPlot from the two 7 x 26 frames through its
    #   documented dot_color_df / dot_size_df inputs; the AnnData it is given
    #   is a scaffold carrying the category order and no expression.
    frames = cache.table(OUTPUT_DIR, "dot_frames", compute_dot_frames)
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
    print(f"  dot frames: {color_df.shape[0]} states x {color_df.shape[1]} genes "
          f"(from data/dot_frames.csv)")
    scaffold = ad.AnnData(
        X=np.zeros((len(cell_type_order), len(plot_genes)), dtype=np.float32),
        obs=pd.DataFrame({"momac_short": pd.Categorical(
            cell_type_order, categories=cell_type_order)},
            index=[f"scaffold_{i}" for i in range(len(cell_type_order))]),
        var=pd.DataFrame(index=plot_genes))

    # 4. Create Dotplot
    print("\n[4/4] Creating dotplot...")
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
        figsize=style.figsize_mm(PANEL_W_MM, PANEL_H_MM),
        dot_color_df=color_df,
        dot_size_df=size_df,
    )
    # The size key is DRAWN, since the author's ruling of 2026-09-11.
    # SUPERSEDED.md beside this script records the divergence from the printed
    # panel. That marker is scoped `Applies-to: E`, because this directory also
    # draws printed panel D, which is untouched.
    #
    # It used to be suppressed with legend(show_size_legend=False), on the
    # reasoning that the published page carried no dot-size key and adding one
    # would be new content. PROVENANCE.csv recorded what that cost and left it
    # standing: "the dot area encodes the fraction of expressing cells in each
    # group, and with no size key nothing on the page decodes it." The ruling
    # is that a key decoding dots already drawn adds no data and is therefore
    # not new content; the dots themselves do not move.
    # THE BIGGEST DOT IS THE PUBLISHED PAGE'S  (2026-09-11)
    #   scanpy sizes its largest dot at 200 pt2 - 14.1 pt, 4.99 mm across.
    #   Counted off 00_GROUND_TRUTH/figures/Figure 4.pdf, the published panel's
    #   largest dot is 2.58 mm, which is 7.31 pt and so 42.1 pt2. At 200 the
    #   dots in the dense columns run into one another and the panel reads as a
    #   block of red rather than as a graded matrix.
    #
    #   cmap is passed again because DotPlot.style() rewrites every style field
    #   it takes, and calling it for largest_dot alone puts the colour map back
    #   to scanpy's default. Figure 5 panel F came out blue-green that way.
    #
    #   Every dot scales by the same factor and the size key is built from
    #   these same areas, so what a dot means does not change.
    # scatter's s is the diameter squared, not pi r^2 - the 42.1 set on
    # 2026-09-11 drew the largest dot 2.29 mm, 11% under the page (2026-09-14).
    dotplot.style(largest_dot=(2.58 * style.PT_PER_MM) ** 2, cmap='Reds')
    # THE KEY COLUMN IS KEY_COLUMN_MM, NOT SCANPY'S 38 MM  (2026-09-14)
    #   scanpy reserves 1.5 in at the right and sets the size key and the
    #   colour bar at opposite ends of it. The author read the paper between
    #   them as waste. The column is reserved at KEY_COLUMN_MM, scanpy's own
    #   key axes are removed, and cnsfig.legend stacks the two keys at the top
    #   of the column, the same way Figures 2G and 5F decode their dots.
    dotplot.legend(width=KEY_COLUMN_MM / 25.4)
    dotplot.make_figure()
    axes = dotplot.get_axes()
    main_ax = axes['mainplot_ax']
    # Gene symbols italic and angled, as the published page sets them. At 45
    # degrees on this 4.0 mm pitch neighbours still overlap by about 1 mm2
    # at 6 pt and at 50 clear by 0.05 pt; 55 is the nearest angle that
    # clears by the 0.5 pt the gate asks for.
    main_ax.set_xticklabels([t.get_text() for t in main_ax.get_xticklabels()],
                            rotation=55, ha='right', style='italic',
                            fontsize=style.tick_pt())

    # scanpy asks matplotlib for the relative size 'small', which is 0.833 of
    # the base size. The sizes are set explicitly instead, from the same two
    # the rest of the figure uses.
    main_ax.tick_params(labelsize=style.tick_pt())

    # De-rasterize all axes (scanpy dotplot may rasterize internally)
    for ax in dotplot.fig.get_axes():
        for coll in ax.collections:
            coll.set_rasterized(False)
        for img in ax.images:
            img.set_rasterized(False)

    fig = dotplot.fig
    # MARGINS ARE SET, NOT FITTED  (2026-09-14)
    #   fit_margins shrank scanpy's grid from the right until the matrix was
    #   53 mm wide and 30 mm of paper stood between it and the key; at that
    #   width 21 gene names at 45 degrees overlap. The margins are the
    #   millimetres the row names, the angled gene names and the key column
    #   need, and the matrix takes the rest: about 85 mm, a 4 mm gene pitch.
    #   The top margin clears the letter cell for the first row name.
    style.margins_mm(fig, left=14.5, right=KEY_COLUMN_MM + 1.5, top=4.0,
                     bottom=10.5)
    # THE MATRIX FILLS THE ROOM  (2026-09-14, evening)
    #   scanpy's gridspec keeps a legend column of its own, so after the
    #   margins were set the matrix stopped 20 mm short of the key ("the
    #   labels and the plot are far apart" - the author). The main axes is
    #   placed by hand: from the left margin to 2.5 mm short of the key.
    _pos = main_ax.get_position()
    _W = PANEL_W_MM
    main_ax.set_position([14.5 / _W, _pos.y0,
                          (_W - 14.5 - KEY_COLUMN_MM - 2.5) / _W, _pos.height])
    key, areas = cnslegend.scanpy_compact_key(
        dotplot, axes, KEY_COLUMN_MM,
        size_title='Cells in\ngroup (%)',
        cbar_title='Mean\nexpression\nin group',   # 'expression' is 9.9 mm; the column 13
        label_pt=style.tick_pt(), title_pt=style.tick_pt())
    print(f"  size key: {len(areas)} steps")
    # Raises if the dot area varies and nothing decodes it. This panel is the
    # reason the check exists.
    cnslegend.require_size_key(key, dot_areas=areas, panel='Figure 4 panel E')
    axes = {'mainplot_ax': main_ax}
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
