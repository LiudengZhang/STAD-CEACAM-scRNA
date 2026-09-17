#!/usr/bin/env python3
"""
Figure 2 panel J - immune-checkpoint gene expression in pre-treatment stomach
epithelium, non-responders against responders.

Each dot is one checkpoint gene: its height is the log2 fold change of the
sample-mean expression, its area the significance band of a two-sided
Mann-Whitney U test over the samples, and its colour the same fold change.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed.

  printed panel  Figure 2 J       (PROVENANCE.csv; NOT inferred from "02_I")

LAYOUT ORDER
    The margins are set before the colour bar is made, because matplotlib sizes
    the bar from the axes box it finds and takes its own share of it. They are
    millimetres of paper, not fractions of the canvas, and they are not fitted
    to the ink afterwards: moving the subplot parameters once the bar exists
    would leave the bar behind.

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the legend - at 5 * SCALE. MARK carries the non-type
    point sizes across to the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The three significance dot areas take AREA; the dot edge width, the zero
    rule, the grid rule and the legend frame take MARK. Tick widths and lengths
    and spine widths do not: those are style, and cnsplots sets them.

LEGEND KEYS
    cnsplots sets `legend.markerscale = 0.5`, which a scatter handler applies
    as an *area* factor of 0.25. Each key is therefore divided by
    `markerscale ** 2`, so that after the legend applies its factor the key
    prints at exactly the area of the dot it labels.

THE GENE LABELS ARE SET AT 80 DEGREES  (2026-09-14)
    The published page sets them at 45. Twenty-nine ticks across this plotting
    box stand 3.3 mm apart; a 6 pt line box is 2.5 mm deep, and at 45 degrees
    neighbours are 2.3 mm apart and overlap whatever their length. At 60
    degrees they are 2.9 mm apart and still graze by 0.9 mm2, at 70 they
    clear by 0.29 pt, at 75 by 0.49, and at 80 by the 0.5 pt the gate asks for. Upright (the 2026-09-11 form) was
    further from the page than the angle the page's own geometry allows.

Every gene, every filter, every fold change, every P value and every string is
the earlier drawing's. The drawing code is the same code.

DRAWING READS A TABLE, since 2026-09-15 (cnsfig.cache): the per-gene fold
change and P value live in data/checkpoint_stats.csv; the h5ad is read only
when the table is absent or --recompute is passed.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

# Central config
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402
import slots                                             # noqa: E402
from cnsfig import cache                                  # noqa: E402

SCALE = 4                           # the earlier canvas multiplier
SMALL_PT = 5.0                      # the earlier smallest body type

PANEL_LETTER = "J"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(2, PANEL_LETTER)
MARGIN = dict(left=6.5, right=0.8, top=2.5, bottom=12.0)

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = BASE_DIR

# Checkpoint genes
INHIBITORY_GENES = [
    'PDCD1', 'CD274', 'CTLA4', 'LAG3',
    'HAVCR2', 'TIGIT', 'CD276', 'CD47', 'IDO1',
    'VSIR', 'BTLA', 'SIGLEC15', 'KLRC1', 'CD96', 'PVRIG', 'VTCN1',
    'NT5E', 'ENTPD1', 'ADORA2A',
    'CEACAM1', 'CEACAM5', 'CEACAM6',
]

COSTIMULATORY_GENES = [
    'ICOS', 'CD27', 'CD70', 'CD40',
    'TNFRSF9', 'TNFRSF4', 'TNFRSF18',
]


#: THE THREE DOT SIZES  (2026-09-14, then 2026-09-15)
#:   Counted off 'Figure 2.pdf' geometry: 42 circles at 1.40 mm (28 of
#:   them), 1.81 mm (10) and 2.21 mm (2) across - the three significance
#:   classes, at the earlier drawing's area ratio 150 : 100 : 60. On the
#:   author's third reading (2026-09-15) the two significant classes did not
#:   stand out enough from the third, so the areas are 1 : 0.55 : 0.30 - the
#:   largest 2.5 mm across, then 1.85 and 1.37 mm. The class boundaries
#:   (0.10, 0.15) and which gene falls in which class are unchanged; only
#:   the three areas, which compare_panel_content.py reports as size_ratios.
#:   On the fifth reading (2026-09-16: "the P > 0.15 dots could be smaller")
#:   the third class halved again, 0.30 -> 0.15 of the largest area (0.97 mm
#:   across); the other two and the class boundaries are unchanged.
#:   scatter's s is the diameter squared in points.
DOT_D_MAX_MM = 2.5
DOT_S_MAX = (DOT_D_MAX_MM * style.PT_PER_MM) ** 2
DOT_AREA_RATIO = {"<=0.10": 1.0, "<=0.15": 0.55, ">0.15": 0.15}


def assign_dot_size(p_value):
    if p_value <= 0.10:
        return DOT_S_MAX * DOT_AREA_RATIO["<=0.10"]
    elif p_value <= 0.15:
        return DOT_S_MAX * DOT_AREA_RATIO["<=0.15"]
    else:
        return DOT_S_MAX * DOT_AREA_RATIO[">0.15"]


def compute_checkpoint_statistics(adata):
    print("Computing checkpoint gene statistics...")
    pre_mask = adata.obs['Treatment phase'] == 'Pre'
    valid_mask = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
    adata_filtered = adata[pre_mask & valid_mask].copy()
    print(f"  Pre-treatment cells: {adata_filtered.n_obs}")

    all_checkpoint_genes = INHIBITORY_GENES + COSTIMULATORY_GENES
    gene_names = adata_filtered.raw.var_names if adata_filtered.raw is not None else adata_filtered.var_names
    available_genes = [g for g in all_checkpoint_genes if g in gene_names]

    results = []
    for gene in available_genes:
        if adata_filtered.raw is not None and gene in adata_filtered.raw.var_names:
            gene_idx = adata_filtered.raw.var_names.get_loc(gene)
            expr = adata_filtered.raw.X[:, gene_idx].toarray().flatten() if hasattr(adata_filtered.raw.X, 'toarray') else adata_filtered.raw.X[:, gene_idx].flatten()
        else:
            gene_idx = adata_filtered.var_names.get_loc(gene)
            expr = adata_filtered.X[:, gene_idx].toarray().flatten() if hasattr(adata_filtered.X, 'toarray') else adata_filtered.X[:, gene_idx].flatten()
            expr = np.nan_to_num(expr, nan=0.0)

        adata_filtered.obs['_gene_expr'] = expr
        sample_means = adata_filtered.obs.groupby('sample', observed=True).agg({
            '_gene_expr': 'mean',
            'stomach_pre_grouping': 'first'
        }).reset_index()

        responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'Responsed']['_gene_expr'].values
        non_responder_vals = sample_means[sample_means['stomach_pre_grouping'] == 'No-response']['_gene_expr'].values

        mean_r = np.mean(responder_vals)
        mean_nr = np.mean(non_responder_vals)
        epsilon = 1e-10
        log2fc = np.log2((mean_nr + epsilon) / (mean_r + epsilon))

        if len(responder_vals) > 1 and len(non_responder_vals) > 1:
            stat, pval = stats.mannwhitneyu(non_responder_vals, responder_vals, alternative='two-sided')
        else:
            pval = 1.0

        category = 'Inhibitory' if gene in INHIBITORY_GENES else 'Costimulatory'
        results.append({'Gene': gene, 'Category': category, 'Merged_Log2FC': log2fc, 'Merged_P_Value': pval})

    return pd.DataFrame(results)


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    def compute():
        print("Loading epithelial data...")
        adata = sc.read_h5ad(EPITHELIAL_H5AD)             # noqa: F405
        print(f"Loaded {adata.n_obs} cells")
        return compute_checkpoint_statistics(adata)

    data = cache.table(BASE_DIR, "checkpoint_stats", compute)
    # Sort by Log2FC for better visualization
    data = data.sort_values('Merged_Log2FC', ascending=False)

    genes = data['Gene'].values
    fold_changes = data['Merged_Log2FC'].values
    p_values = data['Merged_P_Value'].values
    sizes = [assign_dot_size(p) for p in p_values]

    # FLIPPED: genes on X-axis, Log2FC on Y-axis
    x_positions = np.arange(len(genes))

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    # Before the colorbar - see LAYOUT ORDER in the header.
    style.margins_mm(fig, **MARGIN)

    cmap = plt.cm.RdBu_r

    # FLIPPED scatter: x=gene positions, y=fold_changes
    scatter = ax.scatter(
        x_positions, fold_changes, s=sizes, c=fold_changes, cmap=cmap, alpha=0.8,
        edgecolors='black', linewidths=style.EDGE_PT, vmin=-2, vmax=2
    )

    # X-axis: genes
    ax.set_xticks(x_positions)
    # Angled, as the published page sets them (2026-09-14). At 6 pt on this
    # 3.3 mm pitch 45 degrees cannot clear - a 6 pt line box is 2.5 mm deep
    # and 45 degrees leaves 2.3 mm between neighbours - so the angle is 60,
    # 80, the nearest that clears by the 0.5 pt the gate asks for (75: 0.49).
    # Same genes, same order, same strings.
    ax.set_xticklabels(genes, rotation=80, ha='right', style='italic',
                       rotation_mode='anchor')
    ax.set_xlim(-1, len(genes))

    # Y-axis: Log2FC
    y_min, y_max = fold_changes.min(), fold_changes.max()
    y_range = y_max - y_min
    ax.set_ylim(y_min - 0.15 * y_range, y_max + 0.15 * y_range)
    ax.set_ylabel('Log2 FC (NR/R)')

    # Horizontal reference line at 0
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=style.RULE_PT, alpha=0.7)
    # No y grid: the published page draws the zero rule alone (2026-09-14).

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend for significance. See LEGEND KEYS in the header for the
    # markerscale correction.
    key = 1.0 / plt.rcParams["legend.markerscale"] ** 2
    legend_elements = [
        plt.scatter([], [], s=assign_dot_size(0.10)*key, c='gray', alpha=0.6, edgecolors='black', linewidths=style.EDGE_PT, label='p ≤ 0.10'),
        plt.scatter([], [], s=assign_dot_size(0.15)*key, c='gray', alpha=0.6, edgecolors='black', linewidths=style.EDGE_PT, label='p ≤ 0.15'),
        plt.scatter([], [], s=assign_dot_size(0.50)*key, c='gray', alpha=0.6, edgecolors='black', linewidths=style.EDGE_PT, label='p > 0.15')
    ]

    # handletextpad 0.2 -> 0.8 and labelspacing 0.6 (2026-09-15): the keys
    # touched their labels.
    legend = ax.legend(handles=legend_elements, title='Significance', loc='upper right',
                       frameon=True, handletextpad=0.8, borderpad=0.5,
                       labelspacing=0.6, edgecolor='black', framealpha=0.9)
    legend.get_frame().set_linewidth(style.RULE_PT)

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('Log2 FC')

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    style.save_panel(fig, Path(OUTPUT_DIR) / "checkpoint_dotplot_flipped")
    print(f"\nSaved: {Path(OUTPUT_DIR) / 'checkpoint_dotplot_flipped'}"
          f".[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
