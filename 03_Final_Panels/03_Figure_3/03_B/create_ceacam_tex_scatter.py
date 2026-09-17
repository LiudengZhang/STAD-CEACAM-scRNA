#!/usr/bin/env python3
"""
Figure 3 panel D - CEACAM5 and CEACAM6 against the CD8+ Tex fraction, one point
per stomach sample, coloured by treatment phase and response.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. Margins are measured from the rendered ink
rather than typed, and the panel-letter corner is left clear for the assembler.

  printed panel  Figure 3 D       (PROVENANCE.csv; the directory is 03_B, and
                                   the letter was looked up, not inferred)

Every value read, every filter, every statistic and every string is the earlier
drawing's. The Tex state selection, the per-sample percentage, the merge, the
four-group assignment, the Spearman test and the P-value formatting are
untouched, and expression is read exactly as it was read before.

THE COUNTS COME FROM THE EPITHELIAL OBJECT'S OWN LAYER
    The counts were read from a second file holding the same 149,373 cells and
    the same 56,034 genes with the genes in alphabetical order. They are
    layers['counts'] of the epithelial object here, in the gene order .X uses.
    The two agree over every one of their 203,799,816 non-zeros, so the library
    size each cell normalises by is a sum over the same values, and every gene
    below is selected by name rather than by position. The order the counts are
    stored in therefore reaches no number this panel prints.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the shared legend - at 4.5 * SCALE. MARK carries the
    non-type point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The marker area, the marker edge and the dashed regression line are scaled
    by it. Tick widths and lengths and spine widths are not: those are style,
    and cnsplots sets them. The legend's own `markerscale` is the earlier
    drawing's, so the key keeps its size relative to the point it stands for.

THE LEGEND IS INSIDE THE CANVAS, AND THE MARGINS ARE MILLIMETRES
    The earlier drawing anchored the shared legend entirely outside the canvas
    and relied on a tight bounding box at save time to enlarge the saved image
    until the legend fitted. A panel drawn at print size cannot crop to its own
    ink, so a strip of the panel is reserved for the legend and it is anchored
    at the inside of the right edge. Same handles, same labels, same order.

    That strip is why the margins here are millimetres of paper rather than a
    fit to the ink. Fitting the ink moves the axes until the drawing touches
    all four edges, and a legend that is anchored to the figure rather than to
    an axes does not move with it: the axes is then fitted straight through the
    strip. The margins are measured instead, and the two checks that matter -
    that no ink leaves the canvas and that the panel-letter corner is clear -
    are made afterwards and are fatal.

THE COHORT LEAVES THE TITLE, AND THE AXIS LABELS ARE RE-WRAPPED
    The cohort is named once in the caption and each axes keeps its own
    correlation and its own P value; on one line the pair sets 23.6 mm against
    a plotting box of about 10 mm, so they take a line each.

    The x label carries the gene alone; at 7 pt the unit sets 18.0 mm and the
    two plotting boxes are about 7 mm wide and 14 mm apart, so the two axes'
    units print over each other whatever the margins are. The y label keeps the
    page's word order and one word fewer: the label is rotated, so its length
    is vertical, and the page's full form sets 27.1 mm against a 22.9 mm panel.
    Two lines would fit that height but cost 5 mm of width, and at every gutter
    wide enough to take them the right axes' title prints into the shared
    legend instead. Both are declared in RENAMES_FIGURE_3, and the x unit is
    stated in the caption.

MARKER AND RULE WIDTHS ARE MEASURED OFF THE PUBLISHED PAGE  (2026-09-11)
    Every one of these scatter panels carried a stray factor of SCALE on its
    non-type sizes - `s=30*SCALE*AREA`, `linewidth=0.8*SCALE*MARK` - on top of
    AREA and MARK, which already carry the 4x canvas across. The markers came
    out about twice as wide as the page prints them and ran together.

    Counted out of `00_GROUND_TRUTH/figures/Figure 3.pdf` geometry:

        panel A   0.520 mm across, 67 marks      panel B   0.421 mm, 445
        panel D   0.518 mm, 69                   panel E   0.424 mm, 773
        the dashed regression rule               0.595 pt

    The constants below are those numbers converted through this panel's own
    MARK. Not one coordinate, colour or statistic moves.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import TCD8_H5AD, EPITHELIAL_DEPOSIT_H5AD
import panel_style_cns as style  # noqa: E402
import slots  # noqa: E402
from cnsfig import layout as cnslayout, corr_stats, rich_xlabel, rich_ylabel  # noqa: E402
from cnsfig import cache, group_key  # noqa: E402

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 4.5                  # the earlier smallest body type, before * SCALE

MARK = style.tick_pt() / (SMALL_PT * SCALE)   # length multiplier
AREA = MARK ** 2                              # area multiplier

PANEL_LETTER = "D"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(3, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(3, PANEL_LETTER)

PAD_MM = 0.6                    # paper left between the ink and every edge
# Millimetres of paper. left holds the rotated y label and its tick labels,
# bottom the two-line x label and its tick labels, top the two-line title, and
# right the strip the shared legend is drawn in. wspace leaves the right axes
# its own y label and tick column between the two plotting boxes.
#
# wspace was 0.8 until 2026-09-11 and that was 0.13 pt too little. The right
# axes' rotated y label sits in this gutter and the LEFT axes' x label
# ("CEACAM5", wider than its own plotting box) overhangs into it from the
# other side; measured on the shipped panel the two came within 0.13 pt of
# each other and printed as one string. Nothing caught it, because two strings
# that abut share no area and the collision test measured shared area - see
# 10_Reproduction/check_restyled_panel.py, which now measures clearance and
# convicts this panel at its old setting.
# right 14.0 -> 11.5 and wspace 0.8 -> 1.0 on 2026-09-14: the legend is set
# tighter (below) and the two titles' third lines, 17 mm each, need the
# plots 19 mm apart centre to centre.
#: The right column's margin (D over E): E's two-line y label sets it.
LEFT_MM = 12.5   # E uses the same; frames align (2026-09-15)

EPI = EPITHELIAL_DEPOSIT_H5AD
OUTPUT_DIR = Path(__file__).parent

# Standard 4-group palette (blue=R, warm=NR; lighter=Pre, darker=Post)
COLOR_MAP = {
    'Pre-R':  '#74b9ff',
    'Pre-NR': '#e17055',
    'Post-R': '#0984e3',
    'Post-NR':'#d63031',
    'Other':  '#999999',
}

DRAW_ORDER = ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']


def load_primary():
    """Load ALL 32 stomach samples: epithelial CEACAM + CD8 Tex fraction."""
    # --- Epithelial: CEACAM expression ---
    print("[Primary] Loading epithelial counts...")
    epi = sc.read_h5ad(EPI)
    epi.X = epi.layers['counts']
    del epi.layers['counts']
    # uns['log1p'] describes the matrix just replaced. Left in place
    # it makes the normalisation below report the counts as already
    # log-transformed.
    epi.uns.pop('log1p', None)

    # Filter: stomach only (all treatment phases)
    stomach_mask = epi.obs['Sample site'].astype(str).str.lower().str.contains('stomach')
    epi = epi[stomach_mask].copy()

    # Normalize raw counts
    sc.pp.normalize_total(epi, target_sum=1e4)
    sc.pp.log1p(epi)

    epi_df = pd.DataFrame({
        'sample': epi.obs['sample'].values,
        'CEACAM5': epi[:, 'CEACAM5'].X.toarray().flatten() if hasattr(epi[:, 'CEACAM5'].X, 'toarray') else epi[:, 'CEACAM5'].X.flatten(),
        'CEACAM6': epi[:, 'CEACAM6'].X.toarray().flatten() if hasattr(epi[:, 'CEACAM6'].X, 'toarray') else epi[:, 'CEACAM6'].X.flatten(),
        'treatment_phase': epi.obs['Treatment phase'].values,
        'pre_group': epi.obs['stomach_pre_grouping'].values,
        'post_group': epi.obs['stomach_post_grouping'].values,
    })
    epi_sample = epi_df.groupby('sample').agg({
        'CEACAM5': 'mean', 'CEACAM6': 'mean',
        'treatment_phase': 'first',
        'pre_group': 'first',
        'post_group': 'first',
    }).reset_index()
    stomach_samples = set(epi_sample['sample'].values)
    print(f"  {len(epi_sample)} stomach samples")
    del epi

    # --- CD8: Tex fraction ---
    print("[Primary] Loading CD8+ T cells...")
    cd8 = sc.read_h5ad(TCD8_H5AD)
    cd8_df = pd.DataFrame({
        'sample': cd8.obs['sample'].values,
        'minor_cell_state': cd8.obs['minor_cell_state'].values,
    })
    del cd8

    # Filter to stomach samples
    cd8_df = cd8_df[cd8_df['sample'].isin(stomach_samples)]

    cd8_list = []
    for sample, grp in cd8_df.groupby('sample'):
        n_total = len(grp)
        n_tex = (grp['minor_cell_state'] == 'C6_CD8_Tex_PDCD1').sum()
        cd8_list.append({
            'sample': sample,
            'tex_fraction': (n_tex / n_total * 100) if n_total > 0 else 0,
        })
    cd8_sample = pd.DataFrame(cd8_list)
    print(f"  CD8 data for {len(cd8_sample)} samples")

    # Merge
    merged = epi_sample.merge(cd8_sample, on='sample', how='inner')

    # 4-group assignment
    def assign_4group(row):
        phase = str(row['treatment_phase'])
        if phase == 'Pre':
            grp = str(row['pre_group'])
            if grp == 'Responsed': return 'Pre-R'
            if grp == 'No-response': return 'Pre-NR'
        elif phase == 'Post':
            grp = str(row['post_group'])
            if grp == 'Responsed': return 'Post-R'
            if grp == 'No-response': return 'Post-NR'
        return 'Other'

    merged['group'] = merged.apply(assign_4group, axis=1)

    counts = merged['group'].value_counts()
    print(f"  Merged: {len(merged)} samples: {counts.to_dict()}")
    # The cached table carries what the drawing reads and no specimen number.
    return merged[['group', 'CEACAM5', 'CEACAM6', 'tex_fraction']].reset_index(drop=True)


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    # From data/samples.csv (cnsfig.cache, 2026-09-15): the two h5ads are
    # read only when the table is absent or --recompute is passed.
    data = cache.table(OUTPUT_DIR, 'samples', load_primary)

    # ONE GEOMETRY FOR THE FOUR SCATTER PAIRS  (2026-09-14, evening)
    #   The author's ruling: A, B, D and E are the same size, square, the
    #   dataset name alone in the title, rho inside the box, P in the legend
    #   (cnsfig.corr_stats writes it; edits.py reads it), and the two rows
    #   2 mm apart. cnsfig.layout.scatter_pair_mm places the boxes at
    #   millimetres, so the four panels print one geometry by construction.
    fig = style.figure_mm(PANEL_W_MM, PANEL_H_MM)
    axes = cnslayout.scatter_pair_mm(fig, left_mm=LEFT_MM)

    rows = []
    for col_idx, ceacam in enumerate(['CEACAM5', 'CEACAM6']):
        ax = axes[col_idx]
        x = data[ceacam].values.astype(float)
        y = data['tex_fraction'].values.astype(float)
        groups = data['group'].values

        r_val, p_val = stats.spearmanr(x, y)

        for grp in DRAW_ORDER:
            mask = groups == grp
            if mask.sum() > 0:
                ax.scatter(x[mask], y[mask], c=COLOR_MAP[grp], s=19.4*AREA,
                           alpha=0.85, edgecolors='white', linewidths=style.EDGE_PT,
                           label=grp, zorder=3)

        # Regression line
        valid = np.isfinite(x) & np.isfinite(y)
        xv, yv = x[valid], y[valid]
        if len(xv) >= 3 and xv.std() > 0:
            slope, intercept = np.polyfit(xv, yv, 1)
            x_line = np.linspace(xv.min(), xv.max(), 100)
            ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=style.RULE_PT, alpha=0.6, zorder=2)

        p_str = corr_stats.p_string(p_val)
        print(f"    In house: rho = {r_val:.2f}, {p_str}")
        rows.append((ceacam, r_val, p_val, int(valid.sum())))
        ax.set_title('In house', fontsize=style.tick_pt())
        cnslayout.corr_annotate(ax, r_val)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for spine in ['bottom', 'left']:
            ax.spines[spine].set_linewidth(style.RULE_PT)
        ax.tick_params(axis='both', width=style.RULE_PT, length=3*SCALE*MARK)
        rich_xlabel(ax, f"*{ceacam}*")
        # One y label, on the left axes: SUPERSEDED.md beside this script
        # records that as a deliberate divergence from the printed panel,
        # under the author's ruling of 2026-09-11.
        if col_idx == 0:
            # The plus is a plain character, not $^+$: mathtext draws a
            # superscript at 70% of its base, below the 6 pt floor.
            # Two lines (2026-09-14, evening): on one it sets 21.6 mm against
            # a 13.5 mm box and leaves the canvas. Declared in labels.py.
            rich_ylabel(ax, 'CD8+ Tex\nfraction (%)')
    corr_stats.write(OUTPUT_DIR, rows)

    # The shared key, right of the second box. Until 2026-09-15 it was built
    # from the scatter handles, so each key circle printed at the data
    # marker's 0.5 mm and could not be seen; cnsfig.legend.group_key draws
    # them at a fixed 1.3 mm.
    group_key(fig, [(g, COLOR_MAP[g]) for g in DRAW_ORDER],
              x_mm=LEFT_MM + 2 * cnslayout.SCATTER_BOX_MM + cnslayout.SCATTER_GAP_MM + 0.8,
              y_mm=cnslayout.SCATTER_TOP_MM)

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    style.save_panel(fig, OUTPUT_DIR / 'ceacam_tex_scatter')
    print(f"\nSaved: {OUTPUT_DIR / 'ceacam_tex_scatter'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
