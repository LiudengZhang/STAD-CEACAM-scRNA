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
# The widest clean setting is taken, so the size key's dots separate as far as
# the panel allows.
LEGEND_W_MM = 16.0

# The short form of each cell-type name, as the radar panel of this figure
# already prints it. Row order and row contents are untouched.
SHORT_NAME = {
    'DC cells': 'DC',
    'Endothelial cells': 'Endo',
    'Epithelial cells': 'Epi',
    'Fibroblast': 'Fibro',
    'Mast cells': 'Mast',
    'Monocytes/Macrophages': 'MoMac',
    'Neutrophils': 'Neut',
    'NK cells': 'NK',
    'Pericyte': 'Peri',
    'Plasma cells': 'Plasma',
    'CD4+ T cells': 'CD4+ T',
    'CD8+ T cells': 'CD8+ T',
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


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f} (unused here)")

    print("=" * 60)
    print("Panel F: Cytokine Dotplot (TNF, IL6, IL1B, IL1A)")
    print("=" * 60)

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
        print(f"  WARNING — missing: {missing}")

    print("\n[2/3] Creating dotplot...")
    dotplot = sc.pl.dotplot(
        adata,
        var_names=available,
        groupby='major_cell_type',
        categories_order=CELL_TYPE_ORDER,
        use_raw=(adata.raw is not None),
        cmap='Reds',
        show=False,
        save=None,
        standard_scale='var',
        return_fig=True,
        figsize=style.figsize_mm(PANEL_W_MM, PANEL_H_MM),
    )
    # The two legend headings are centred over that column, so at full
    # length they reach back across the dot matrix and are painted over
    # by it: "Fraction of cells" sets 17.2 mm at 7 pt against a 14 mm
    # column. Shortened, they sit inside it.
    dotplot.legend(width=LEGEND_W_MM / 25.4,
                   size_title='Cells in\ngroup (%)',
                   colorbar_title='Mean expr.\nin group')
    dotplot.make_figure()
    axes = dotplot.get_axes()

    main_ax = axes['mainplot_ax']
    main_ax.set_yticklabels(
        [SHORT_NAME.get(t.get_text(), t.get_text())
         for t in main_ax.get_yticklabels()])

    # scanpy asks matplotlib for the relative size 'small', which is 0.833 of
    # the base size and put every string on this panel at 5.83 pt. The sizes
    # are set explicitly instead, from the same two the rest of the figure uses.
    for one in axes.values():
        one.tick_params(labelsize=style.tick_pt())
        if one.get_title():
            one.set_title(one.get_title(), fontsize=style.body_pt())

    # The size key names every other dot. Six thresholds need 14.1 mm of label
    # track and the key is 8.9 mm wide at its widest legible setting, so all
    # six print on top of one another; three name the same scale and read.
    size_ax = axes.get('size_legend_ax')
    if size_ax is not None:
        size_ax.set_xticks(list(size_ax.get_xticks())[1::2])

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

    print("\n[3/3] Saving...")
    style.save_panel(fig, OUTPUT_DIR / 'cytokine_dotplot')
    print(f"  Saved: cytokine_dotplot.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")
    print("Done!")


if __name__ == "__main__":
    main()
