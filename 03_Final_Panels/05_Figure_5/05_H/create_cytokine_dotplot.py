#!/usr/bin/env python3
"""
Figure 5 panel F - dotplot of TNF, IL6, IL1B and IL1A across the thirteen major
cell types.

The gene list, the cell-type order, the `major_cell_type` filter (Hepatocyte
excluded), `use_raw`, `cmap='Reds'` and `standard_scale='var'` are unchanged.
`sc.pl.dotplot` still computes and draws everything.

  printed panel  Figure 5 F       (PROVENANCE.csv - the directory is "05_H";
                                   do NOT read the directory as the letter)

The thirteen cell-type names are printed in the short forms the radar panel of
this same figure already uses, so that one figure names a cell type one way.
The rows and their order do not move; the mapping is declared in
00_Config/shared/labels.py.

VALUES: recomputed from `layers['counts']`, author's ruling 2026-09-09
    `full_dataset.h5ad` has no `.raw`, because build_clean_h5ad.py promoted
    `.raw.X` to `.X`; and for this file that `.raw` was the doubly normalised
    matrix of 00_Data_Audit/FINDINGS.md 12.11-12.12. So the panel as shipped
    drew its colours from a damaged matrix. That was measured, not argued: the
    52 dot fills of the previous SVG were read back and the `Reds` ramp
    inverted, and they reproduce the group means this file's `.X` gives.

    The integer counts are intact, so the sound values are recovered without
    rebuilding anything - but NOT from the file the old `FULL_DATASET_H5AD`
    resolved to, which carries no counts layer. See the INPUT note below the
    imports: the object read now is that same matrix, proved identical cell by
    cell, in the deposited container that carries the counts.
    `target_sum=1e4` with `log1p` is the project's own transform,
    read out of FINDINGS 12.12 rather than assumed from scanpy's default, and
    confirmed exactly: redoing it from the counts layer of the deposited
    `Neutrophils_sound.h5ad` reproduces that file's `.X` to max |diff| 0.000000.

    What this moves, and only this:
      dot size   percentage of cells expressing, which depends only on the zero
                 pattern. A monotone transform cannot move it, and measured it
                 does not: max |delta pct| = 0.0000000000 over 13 groups x 4 genes.
      colour     each group's mean, scaled per gene by `standard_scale='var'`.
                 Moves by at most 0.069 (TNF), 0.121 (IL6), 0.090 (IL1B),
                 0.094 (IL1A). IL6's highest group becomes Fibroblast where the
                 shipped panel printed B cells; the highest group for TNF (DC),
                 IL1B and IL1A (both MoMac) is unchanged.

    Not one string, gene, cell type, row or filter changes. The prior version
    is archived under 07_Archive/2026-09-09_fig5F_sound_redraw/.

    THIS PANEL DOES NOT REPRODUCE THE PRINTED FIGURE, AND IS NOT MEANT TO.
    PROVENANCE.csv records it as reproduces_published = no. That is the
    *superseded* sense of no, not the broken one, so the marker beside this
    script is SUPERSEDED.md and there is no KNOWN_BROKEN.md and no Retest-by
    date - there is nothing to retest. SUPERSEDED.md records what was
    superseded, why, on whose ruling and on what date; RULES.md rule 1 sets
    out the two senses and verify_panel_provenance.py checks 3 and 10 enforce
    them.

    Since 2026-09-10 the deposited object's own `.X` is the same matrix this
    recomputes: it was rebuilt from the same counts layer by the same transform
    (07_Archive/2026-09-10_before_X_recompute_from_counts/). The counts layer is
    read here anyway, and deliberately. It is the ground truth `.X` is derived
    from, so reading it is one step of provenance rather than two; it fails
    loudly if the counts are absent instead of silently drawing whatever `.X`
    happens to hold; and it is the route that survives a future re-damage of
    `.X`, which is the failure that reached this panel once already. Measured,
    the choice moves nothing: dot sizes and scaled colours agree to
    0.0000000000 either way.

MARK
    This is the one Figure 5 panel with no `SCALE` and no `N * SCALE` anywhere:
    the earlier drawing set only the font family and the font types, so every
    string it drew came out at matplotlib's default 10 pt on a canvas scanpy
    sized itself. That reads as SCALE = 1 and SMALL_PT = 10, so

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
"""
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FULL_DATASET_DEPOSIT_H5AD
import panel_style_cns as style
import slots
from cnsfig import legend as cnslegend
from cnsfig import cache
import anndata as ad
import numpy as np
import pandas as pd

# INPUT. FULL_DATASET_DEPOSIT_H5AD, not FULL_DATASET_H5AD, and it is declared
# in 00_Config/paths.py beside EPITHELIAL_DEPOSIT_H5AD and
# NEUTROPHILS_SOUND_H5AD rather than as a literal here - a Path literal into
# 06_Clean_Data/ is a working-tree address that is in no deposited record, and
# the release rebuilds the paths module instead of the panels.
# The reasoning for the repoint is in paths.py at that constant.

OUTPUT_DIR = Path(__file__).resolve().parent

SCALE = 1                       # the earlier drawing had no SCALE: 1:1 at 10 pt
SMALL_PT = 10.0                 # matplotlib's default font.size, left unchanged

# The project's own normalisation, read out of 00_Data_Audit/FINDINGS.md 12.12,
# not taken from scanpy's default. See the VALUES note in the docstring.
TARGET_SUM = 1e4

PANEL_LETTER = "F"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(5, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(5, PANEL_LETTER)

# scanpy reserves a fixed-width column at the right of the figure for the size
# key and the colour bar, 1.5 inches by default - 38.1 mm of this 38.5 mm
# panel, which leaves the dot matrix a single column of pixels. Measured on
# this panel, the matrix and the key are both clean between 12 and 16 mm and
# neither is outside that band: below 12 mm the key's thresholds print on top
# of one another, and above 16 mm IL1B and IL1A collide as the matrix narrows.
#
# The NARROWEST clean setting is taken, since 2026-09-11, not the widest. The
# four millimetres it gives back are what let the rows carry the cell-type
# names the published page prints instead of five-letter abbreviations. The
# key still lays its circles out on a pitch taken from their own radii, so
# nothing in it merges; it simply has less room to spread into.
#: The row names' column, the key column under the matrix, and the band
#: below the matrix (gene names GENE_MM, then the key).
LEFT_MM = 25.0          # Monocytes/Macrophages on one line: 22 mm at 6 pt, ticks and pad
SIZE_KEY_W_MM = 16.0
GENE_MM = 5.6
BOTTOM_MM = 16.5

# THE ROW NAMES ARE THE PUBLISHED ONES AGAIN  (2026-09-11)
#   Twelve of the thirteen were abbreviated on 2026-09-10 - Endothelial cells
#   to Endo, Monocytes/Macrophages to MoMac - because the row labels, the dot
#   matrix and a 16 mm legend column would not all fit across 38.5 mm at 6 pt.
#   The panel is the same width today; what changed is that the legend column
#   is 12 mm, which is its measured lower bound (see LEGEND_W_MM), and that
#   frees the 4 mm the long names needed.
#
#   Exactly one name still does not fit on a line: Monocytes/Macrophages sets
#   23.7 mm at 6 pt. It is wrapped at the solidus rather than cut, so every
#   word the page prints is still printed. Row order and row contents are
#   untouched, and no name is replaced by a different word.
#   The three names wider than the 12.5 mm the row column can hold are
#   wrapped, at their own space or solidus. Everything else prints on one line.
#   Measured at 6 pt: Monocytes/Macrophages 23.7 mm, Endothelial cells 15.4,
#   Epithelial cells 13.5; the next widest, CD4+ T cells, is 12.2 and stays.
SHORT_NAME = {
    'Monocytes/Macrophages': 'Monocytes/\nMacrophages',
    'Endothelial cells':     'Endothelial\ncells',
    'Epithelial cells':      'Epithelial\ncells',
}

# Genes of interest
GENES = ['TNF', 'IL6', 'IL1B', 'IL1A']

# Cell type order (matching 05_F radar, using original adata labels)
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
    """The computing half: scanpy's own mean (standard_scale='var') and
    fraction-expressing matrices, as one long table."""
    print("\n[1/3] Loading data...")
    adata = sc.read_h5ad(FULL_DATASET_DEPOSIT_H5AD)
    print(f"  {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # Filter to 13 cell types (exclude Hepatocyte)
    mask = adata.obs['major_cell_type'].isin(CELL_TYPE_ORDER)
    adata = adata[mask].copy()
    print(f"  After filtering: {adata.shape[0]:,} cells, {adata.obs['major_cell_type'].nunique()} types")

    # The values, from the counts layer. `.X` as loaded is the doubly
    # normalised matrix; the counts beside it are intact. normalize_total
    # divides each cell by its own row sum, so doing this after the cell
    # filter gives the same per-cell values as doing it before.
    if 'counts' not in adata.layers:
        raise RuntimeError(
            "layers['counts'] is absent; this panel is drawn from the counts "
            "layer by the author's ruling of 2026-09-09 and must not silently "
            "fall back to .X")
    adata.X = adata.layers['counts'].copy()
    adata.uns.pop('log1p', None)        # or sc.pp.log1p refuses to record its base
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    del adata.layers['counts']
    print(f"  values recomputed from layers['counts'] "
          f"(normalize_total target_sum={TARGET_SUM:g}, log1p)")

    # Validate genes
    var_names = adata.raw.var_names if adata.raw is not None else adata.var_names
    available = [g for g in GENES if g in var_names]
    missing = [g for g in GENES if g not in var_names]
    print(f"  Genes available: {available}")
    if missing:
        print(f"  WARNING - missing: {missing}")

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
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f} (unused here)")

    print("=" * 60)
    print("Panel F: Cytokine Dotplot (TNF, IL6, IL1B, IL1A)")
    print("=" * 60)

    # THE MATRICES ARE READ FROM data/, since 2026-09-15 (cnsfig.cache).
    #   The 75-second half of this script - opening the deposit h5ad,
    #   re-normalising from the counts layer, summarising 13 x 4 means and
    #   fractions - runs only when data/dot_frames.csv is absent or
    #   --recompute is passed. The drawing half builds scanpy's DotPlot from
    #   those two 13 x 4 frames through its documented dot_color_df /
    #   dot_size_df inputs (the standard way to hand DotPlot precomputed
    #   values); the AnnData it is given is a 13-row scaffold that carries
    #   the category order and nothing else - no expression is read from it.
    frames = cache.table(OUTPUT_DIR, "dot_frames", compute_dot_frames)
    available = [g for g in GENES if g in set(frames["gene"])]
    color_df = frames.pivot(index="cell_type", columns="gene", values="mean_scaled") \
                     .loc[CELL_TYPE_ORDER, available]
    size_df = frames.pivot(index="cell_type", columns="gene", values="fraction") \
                    .loc[CELL_TYPE_ORDER, available]
    print(f"  dot frames: {color_df.shape[0]} cell types x {color_df.shape[1]} genes "
          f"(from data/dot_frames.csv)")
    scaffold = ad.AnnData(
        X=np.zeros((len(CELL_TYPE_ORDER), len(available)), dtype=np.float32),
        obs=pd.DataFrame({"major_cell_type": pd.Categorical(
            CELL_TYPE_ORDER, categories=CELL_TYPE_ORDER)},
            index=[f"scaffold_{i}" for i in range(len(CELL_TYPE_ORDER))]),
        var=pd.DataFrame(index=available))

    print("\n[2/3] Creating dotplot...")
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
        figsize=style.figsize_mm(PANEL_W_MM, PANEL_H_MM),
        dot_color_df=color_df,
        dot_size_df=size_df,
    )
    # The two legend headings are centred over that column, so at full
    # length they reach back across the dot matrix and are painted over
    # by it: "Fraction of cells" sets 17.2 mm at 7 pt against a 14 mm
    # column. Shortened, they sit inside it.
    # THE BIGGEST DOT IS THE PUBLISHED PAGE'S, NOT SCANPY'S DEFAULT
    #   scanpy sizes its largest dot at 200 pt2 - 14.1 pt, 4.99 mm across -
    #   whatever width the matrix ends up with. Counted off
    #   00_GROUND_TRUTH/figures/Figure 5.pdf, the published panel's largest dot
    #   is 2.63 mm, and 200 pt2 in a four-column matrix this narrow overflows
    #   its own column and touches the dot beside it. 2.63 mm is 7.45 pt, so
    #   the area is 43.6 pt2.
    #
    #   Every dot is scaled by the same factor and the size key is built from
    #   these same areas by fit_size_key, so what a dot means is unchanged.
    #
    #   cmap is passed again here on purpose. scanpy's DotPlot.style() rewrites
    #   every style field it takes, so calling it for largest_dot alone put the
    #   colour map back to scanpy's default and the panel came out blue-green
    #   instead of Reds. Caught by looking at the rebuilt panel; nothing in the
    #   pipeline checks a colour map.
    # scatter's s is the diameter squared, not pi r^2 - the 43.6 set on
    # 2026-09-11 drew the largest dot 2.33 mm, 11% under the page (2026-09-14).
    dotplot.style(largest_dot=(2.63 * style.PT_PER_MM) ** 2, cmap='Reds')
    # The column is reserved at LEGEND_W_MM; scanpy's own key axes are then
    # removed and cnsfig.legend stacks the two keys at its top, the same way
    # Figures 2G and 4E decode their dots (2026-09-14).
    # THE KEY GOES UNDER THE MATRIX  (2026-09-14, evening)
    #   Beside it, the key column left a four-column matrix 10 mm wide and
    #   its own headings painted across the last column ("F's problem is
    #   obvious" - the author). scanpy reserves no legend column now, the
    #   matrix takes the width between the row names and the right edge,
    #   and cnsfig.legend.compact_key stacks the size key and the colour bar
    #   in a KEY_W_MM column under the gene names, flush with the matrix.
    dotplot.legend(show=False)
    dotplot.make_figure()
    axes = dotplot.get_axes()

    main_ax = axes['mainplot_ax']
    # All thirteen names on one line (2026-09-14, evening): the column is
    # 23 mm, 'Monocytes/Macrophages' sets 22, and one-line rows are what let
    # the matrix keep a 3.8 mm pitch with the key underneath.
    # Gene symbols italic, as the published page sets them.
    main_ax.set_xticklabels([t.get_text() for t in main_ax.get_xticklabels()],
                            rotation=90, style='italic')

    # scanpy asks matplotlib for the relative size 'small', which is 0.833 of
    # the base size and put every string on this panel at 5.83 pt. The sizes
    # are set explicitly instead, from the same two the rest of the figure uses.
    main_ax.tick_params(labelsize=style.tick_pt())

    # The size key is redrawn by 00_Config/cnsfig/legend.py, which lays the
    # circles out on a pitch taken from their own radii.
    #
    # What shipped before: scanpy places the key's circles at
    # `np.arange(n) + 0.5` - one data unit apart - and sizes them from
    # largest_dot = 200 pt^2, a circle 15.96 pt across. In this 16 mm column
    # one data unit is about 4.1 pt, so a 16 pt circle was drawn every 4.1 pt
    # and the six of them fused into a solid grey wedge. Measured on the
    # shipped panel, the last two circles overlapped by 3.41 pt. Thinning the
    # tick labels - what this code used to do - fixed the labels and left the
    # wedge, because the labels were never the problem.
    #
    # The column cannot hold the key across, and the circles may not be shrunk
    # to make it fit: a key drawn at a different scale from the dots it
    # decodes is worse than none. So it goes down the column, which has room.
    fig = dotplot.fig
    style.margins_mm(fig, left=LEFT_MM, right=0.6, top=4.3, bottom=BOTTOM_MM)
    # scanpy's gridspec keeps columns of its own beside the matrix; the
    # matrix is placed by hand from the row-name column to the right edge
    # (the same fix as Figure 4E), which is what gives the four gene names
    # a 3.2 mm pitch instead of 2.3.
    _pos = main_ax.get_position()
    main_ax.set_position([LEFT_MM / PANEL_W_MM, _pos.y0,
                          (PANEL_W_MM - LEFT_MM - 0.6) / PANEL_W_MM, _pos.height])
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
    # Two keys side by side under the matrix: the size key in a
    # SIZE_KEY_W_MM column at the left edge, the colour bar in the rest of
    # the width, each with its heading.
    # The band runs the canvas width: the row-name column is empty below
    # the last row, and the matrix alone is 15 mm wide.
    x_left = 0.6
    y_top = PANEL_H_MM - BOTTOM_MM + GENE_MM
    W, H = PANEL_W_MM, PANEL_H_MM
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
    print(f"  size key: {len(areas)} steps")
    cnslegend.require_size_key(key, dot_areas=areas,
                               panel='Figure 5 panel F')

    # scanpy rasterises the colour bar's solids, and a raster in a vector panel
    # is a raster on the page: the gate counts <image> elements and expects
    # none here.
    for one in dotplot.fig.get_axes():
        for coll in one.collections:
            coll.set_rasterized(False)
        for img in one.images:
            img.set_rasterized(False)

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    print("\n[3/3] Saving...")
    style.save_panel(fig, OUTPUT_DIR / 'cytokine_dotplot')
    print(f"  Saved: cytokine_dotplot.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")
    print("Done!")


if __name__ == "__main__":
    main()
