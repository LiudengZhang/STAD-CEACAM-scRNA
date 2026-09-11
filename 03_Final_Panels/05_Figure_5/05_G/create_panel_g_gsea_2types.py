#!/usr/bin/env python3
"""
Figure 5 panel I - GSEA enrichment plots for the TNF-alpha/NF-kB Hallmark set
in fibroblasts and epithelial cells.

The two h5ads, the `stomach_post_grouping` R-vs-NR contrast, the 5000-cell
subsample at `random_state=42`, `t-test_overestim_var`, the
`logfoldchanges * -log10(p)` ranking metric, `min_size=5`, `max_size=500`,
`permutation_num=1000`, `seed=42` and the "p < 0.001" bound formatting are
unchanged. Not one of them is touched by a restyle.

  printed panel  Figure 5 I       (PROVENANCE.csv - the directory is "05_G";
                                   do NOT read the directory as the letter)

WHAT THIS PANEL DRAWS, AND WHY IT IS TWO AXES AND NOT FOUR
----------------------------------------------------------
The published panel prints two enrichment curves and nothing else. Measured off
00_GROUND_TRUTH/figures/Figure 5.pdf at 600 dpi, the ink inside panel I's
rectangle falls in four bands - 81.5-92.2 mm (the fibroblast curve), 92.3-94.1
(its gene-rank tick labels), 96.9-107.9 (the epithelial curve), 108.1-109.9
(its tick labels) - with 94.1-96.9 mm empty, and the extracted text holds no
"Hits" and no "Gene Rank" anywhere in the rectangle. There is no gene-hit
strip on the page.

The panel script drew four axes: a curve and a hit strip for each cell type.
The published figure is the ground truth and the code is the suspect, so the
two strips are not drawn here. They are declared in
00_Config/shared/labels.py as strings the published panel does not print.

The two curves are drawn over their own gene rank, not a shared one: the two
ranked lists are 57058 and 56034 genes long, so one shared x limit would
redraw the fibroblast curve against the epithelial rank. The score is named
once, beside both curves, and that label is set in full - the published panel
clips it to "Enrichment Scor", a glyph short, and the string the script has
always set is "Enrichment Score".

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the tick labels and the boxed NES/p/FDR annotation -
    at 5 * SCALE, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The ES curve, the zero rule, the peak rule, the annotation box edge and the
    hit rules are scaled by it. Tick widths and lengths and spine widths are
    not: those are style, and cnsplots sets them.

NOTE FOR THE GATE
    `gp.prerank` is given an `outdir`, so gseapy writes its own tables and
    report plots there. `compare_panel_content.py` blocks disk writes, so this
    panel must be gated with `--divert-writes`, which swallows and lists them.
"""
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gseapy as gp
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *
import panel_style_cns as style
import slots

BASE_DIR = Path(__file__).parent

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 5.0                  # the earlier smallest body type, before * SCALE
MARK = style.tick_pt() / (SMALL_PT * SCALE)

PANEL_LETTER = "I"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(5, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(5, PANEL_LETTER)

PATHWAY_NAME = "TNF-alpha Signaling via NF-kB"

CELL_TYPE_CONFIG = {
    "Fibroblast": {"path": FIBROBLAST_H5AD, "label": "Fibroblast"},
    "Epithelial": {"path": EPITHELIAL_H5AD, "label": "Epithelial"},
}

COMPARISON = {"column": "stomach_post_grouping", "group1": "Responsed", "group2": "No-response"}


def run_deg_ttest(adata, comparison_info, max_cells=5000):
    col = comparison_info["column"]
    g1, g2 = comparison_info["group1"], comparison_info["group2"]
    if col not in adata.obs.columns:
        return None
    valid_mask = adata.obs[col].isin([g1, g2])
    adata_sub = adata[valid_mask].copy()
    n_g1 = (adata_sub.obs[col] == g1).sum()
    n_g2 = (adata_sub.obs[col] == g2).sum()
    print(f"    R: {n_g1}, NR: {n_g2}")
    if n_g1 < 10 or n_g2 < 10:
        return None
    if adata_sub.n_obs > max_cells:
        sc.pp.subsample(adata_sub, n_obs=max_cells, random_state=42)
    try:
        sc.tl.rank_genes_groups(adata_sub, groupby=col, groups=[g2], reference=g1,
                                method='t-test_overestim_var', pts=True)
        return sc.get.rank_genes_groups_df(adata_sub, group=g2)
    except Exception as e:
        print(f"    DEG error: {e}")
        return None


def run_gsea(deg_df, cell_type):
    if deg_df is None or len(deg_df) == 0:
        return None
    deg_df = deg_df.copy()
    deg_df['rank_metric'] = deg_df['logfoldchanges'] * (-np.log10(deg_df['pvals'].clip(lower=1e-300)))
    rnk = deg_df.set_index('names')['rank_metric'].dropna()
    rnk = rnk[~rnk.index.duplicated(keep='first')].sort_values(ascending=False)
    if len(rnk) < 15:
        return None
    gsea_outdir = BASE_DIR / f"gsea_{cell_type}"
    gsea_outdir.mkdir(parents=True, exist_ok=True)
    try:
        pre_res = gp.prerank(rnk=rnk, gene_sets=str(HALLMARK_GMT),
                             outdir=str(gsea_outdir), min_size=5, max_size=500,
                             permutation_num=1000, seed=42, verbose=False)
        return pre_res
    except Exception as e:
        print(f"    GSEA error: {e}")
        return None


def draw_gsea_subplot(ax_top, pre_res, label):
    """Draw GSEA enrichment on two axes (ES curve + gene hits)."""
    term_data = pre_res.results[PATHWAY_NAME]
    RES = term_data['RES']
    hits = term_data['hits']
    nes = term_data['nes']
    pval = term_data['pval']
    fdr = term_data['fdr']

    x = np.arange(len(RES))
    color = '#C62828'

    # Top: running ES
    ax_top.plot(x, RES, color=color, linewidth=1.5 * MARK)
    ax_top.axhline(y=0, color='gray', linestyle='--', linewidth=0.5 * MARK)
    ax_top.fill_between(x, 0, RES, where=(np.array(RES) >= 0), color='#EF9A9A', alpha=0.3)
    ax_top.fill_between(x, 0, RES, where=(np.array(RES) < 0), color='#90CAF9', alpha=0.3)

    peak_idx = np.argmax(np.abs(RES))
    ax_top.axvline(x=peak_idx, color='red', linestyle=':', linewidth=0.8 * MARK, alpha=0.7)

    ax_top.set_xlim(0, len(RES))
    ax_top.set_title(label)

    # The printed panel reads "p < 0.001", not "p = 0.000". A permutation P
    # value below the resolution of the test is bounded, never zero, and three
    # decimal places on their own turn the bound into a false exact value.
    def fmt(name, v):
        return f"{name} < 0.001" if v < 0.001 else f"{name} = {v:.3f}"

    stats_text = f"NES = {nes:.2f}\n{fmt('p', pval)}\n{fmt('FDR', fdr)}"
    ax_top.text(0.98, 0.95, stats_text, transform=ax_top.transAxes,
                fontsize=style.tick_pt(), va='top', ha='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray',
                          linewidth=0.5 * MARK))

    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)

    print(f"    {label}: NES={nes:.2f}, p={pval:.3f}, FDR={fdr:.3f}, "
          f"{len(hits)} leading-edge hits")


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("=" * 60)
    print("Panel I: GSEA Enrichment — Fibroblast + Epithelial")
    print("  (t-test DEGs, min_size=5)")
    print("=" * 60)

    # Run GSEA for each cell type
    gsea_results = {}
    for cell_type, config in CELL_TYPE_CONFIG.items():
        print(f"\n=== {cell_type} ===")
        if not config["path"].exists():
            print(f"  Not found: {config['path']}")
            continue
        adata = sc.read_h5ad(config["path"])
        adata.obs_names_make_unique()
        print(f"  {adata.n_obs} cells")

        deg_df = run_deg_ttest(adata, COMPARISON)
        if deg_df is None:
            continue
        pre_res = run_gsea(deg_df, cell_type)
        if pre_res is not None and PATHWAY_NAME in pre_res.results:
            gsea_results[cell_type] = pre_res

    if not gsea_results:
        print("No GSEA results! Exiting.")
        return

    # Draw directly in subplots — 2 cell types × (ES + hits) = 4 rows
    n_types = len(gsea_results)
    # The axes are NOT shared. The two ranked gene lists are different
    # lengths, so one x limit for both would redraw the fibroblast curve
    # against the epithelial rank; the label is stated once, the scales are not.
    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, n_types, 1)
    axes = np.atleast_1d(axes)
    # The gap between the two curves has to hold the upper one's rank labels
    # and the lower one's title, 2.11 and 2.47 mm at these sizes. Measured, the
    # two clear each other from 0.85 of an axes height and touch below it.
    # fit_margins moves the four outer edges and keeps this spacing.
    fig.subplots_adjust(hspace=0.85)

    for idx, (cell_type, pre_res) in enumerate(gsea_results.items()):
        label = CELL_TYPE_CONFIG[cell_type]["label"]
        draw_gsea_subplot(axes[idx], pre_res, label)

    # One label for two axes. At 7 pt it stands 20.5 mm tall against an 11 mm
    # axes, which is why the published panel clips it to "Enrichment Scor"; it
    # is therefore anchored at the boundary between the two curves and left to
    # span both. Its `y` is the position along the lower axes, so matplotlib
    # still places it to the left of the tick labels itself and the fit still
    # reserves the column it needs. `y` is walked down until the label clears
    # the corner the panel letter is drawn into.
    shared = axes[-1]
    for step in range(10):
        # `ha` is what centres a rotated label along the axis; `va`
        # would move it sideways, over the tick labels the axis
        # places it clear of.
        shared.set_ylabel('Enrichment Score', y=1.0 - 0.05 * step)
        try:
            style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
        except RuntimeError as exc:
            if "letter" not in str(exc) or step == 9:
                raise
            continue
        break

    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    style.save_panel(fig, BASE_DIR / 'panel_g_gsea_4types')
    print(f"\nSaved: panel_g_gsea_4types.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
