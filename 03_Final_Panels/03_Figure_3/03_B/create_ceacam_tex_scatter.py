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
MARGIN = dict(left=10.0, right=14.0, top=6.5, bottom=7.4, wspace=0.8)

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
    return merged


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    data = load_primary()

    # 1x2 layout
    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, 1, 2)

    for col_idx, ceacam in enumerate(['CEACAM5', 'CEACAM6']):
        ax = axes[col_idx]
        x = data[ceacam].values.astype(float)
        y = data['tex_fraction'].values.astype(float)
        groups = data['group'].values

        r_val, p_val = stats.spearmanr(x, y)

        for grp in DRAW_ORDER:
            mask = groups == grp
            if mask.sum() > 0:
                ax.scatter(x[mask], y[mask], c=COLOR_MAP[grp], s=30*SCALE*AREA,
                           alpha=0.85, edgecolors='white', linewidths=0.3*SCALE*MARK,
                           label=grp, zorder=3)

        # Regression line
        valid = np.isfinite(x) & np.isfinite(y)
        xv, yv = x[valid], y[valid]
        if len(xv) >= 3 and xv.std() > 0:
            slope, intercept = np.polyfit(xv, yv, 1)
            x_line = np.linspace(xv.min(), xv.max(), 100)
            ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=0.8*SCALE*MARK, alpha=0.6, zorder=2)

        # Stats text — 1 sig digit (floor). A very small P is written out rather
        # than set as a mathtext power of ten: mathtext draws a superscript at
        # 70% of its base, so an exponent on a 7 pt title prints at 4.9 pt,
        # below the floor this figure set is set to.
        import math
        _e = math.floor(math.log10(p_val)); _c = int(p_val / 10**_e)
        if _e >= -3:
            p_str = f'P = {_c * 10**_e:.{-_e}f}'
        else:
            p_str = f'P = {_c}e{_e}'

        # The cohort is named in the caption; see THE COHORT LEAVES THE TITLE.
        print(f"    Primary Cohort (scRNA-seq): rho = {r_val:.2f}, {p_str}")
        ax.set_title(f'ρ = {r_val:.2f}\n{p_str}', linespacing=1.4)

        ax.set_xlabel(ceacam)
        ax.set_ylabel('Tex in CD8+ (%)')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for spine in ['bottom', 'left']:
            ax.spines[spine].set_linewidth(0.5*SCALE*MARK)
        ax.tick_params(axis='both', width=0.5*SCALE*MARK, length=3*SCALE*MARK)

        ax.set_box_aspect(1)

    # The shared legend, in the reserved strip at the right edge.
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right',
               bbox_to_anchor=(1.0 - PAD_MM / PANEL_W_MM, 0.5),
               framealpha=0, edgecolor='none', markerscale=0.8)

    style.margins_mm(fig, **MARGIN)
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
