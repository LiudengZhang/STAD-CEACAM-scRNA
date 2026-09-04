#!/usr/bin/env python3
"""
Figure 5 panel I, RESTYLED (Version B) - GSEA enrichment plots for the
TNF-alpha/NF-kB Hallmark set in fibroblasts and epithelial cells.

Version A is
`03_Revised_Panels/Main_Figures/05_Figure_5/05_G/create_panel_g_gsea_2types.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

The two h5ads, the `stomach_post_grouping` R-vs-NR contrast, the 5000-cell
subsample at `random_state=42`, `t-test_overestim_var`, the
`logfoldchanges * -log10(p)` ranking metric, `min_size=5`, `max_size=500`,
`permutation_num=1000`, `seed=42` and the "p < 0.001" bound formatting are
Version A's, unchanged. Not one of them is touched by a restyle.

  printed panel  Figure 5 I        (PROVENANCE.csv - the directory is "05_G";
                 do NOT read the directory as the letter)
  printed rect   44.6 x 29.6 mm    (panel_rects.csv)
  Version B box  62.0 x 104.0 mm

    Four stacked axes - a running-ES curve and a hit strip for each of two cell
    types - each need their own y axis label at 8 pt, and the ES axes carry a
    three-line boxed statistics annotation. 29.6 mm of height cannot hold four
    axes at 7 pt. The inter-axes gap (`hspace`) also has to grow: at 7/8 pt
    the first hit strip's "Gene Rank" label and the second curve's "Epithelial"
    title collided at Version A's spacing.

MARK
    Version A drew 24 x 28 cm at SCALE = 4 and set its smallest body type - the
    tick labels and the boxed NES/p/FDR annotation - at `5 * SCALE`, so
    SMALL_PT = 5 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 20 = 0.35

    The ES curve, the zero rule, the peak rule, the annotation box edge and the
    hit rules are scaled by it. Tick widths/lengths and spine widths are not:
    those are style, and cnsplots sets them.

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
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *
import panel_style_cns as style

BASE_DIR = Path(__file__).parent

SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 5.0                  # Version A's smallest body type, before * SCALE
MARK = style.tick_pt() / (SMALL_PT * SCALE)

PRINTED_MM = (44.6, 29.6)       # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 62.0, 104.0
MARGIN = dict(left=14.0, right=2.0, top=6.0, bottom=10.0, hspace=1.25)

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


def draw_gsea_subplot(ax_top, ax_bot, pre_res, label):
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

    ax_top.set_ylabel('Enrichment Score')
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

    # Bottom: gene hits
    for hit in hits:
        ax_bot.axvline(x=hit, color='black', linewidth=0.5 * MARK, alpha=0.5)
    ax_bot.set_xlim(0, len(RES))
    ax_bot.set_xlabel('Gene Rank')
    ax_bot.set_ylabel('Hits')
    ax_bot.set_yticks([])
    ax_bot.spines['top'].set_visible(False)
    ax_bot.spines['right'].set_visible(False)

    print(f"    {label}: NES={nes:.2f}, p={pval:.3f}, FDR={fdr:.3f}")


def main():
    family = style.apply()
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("=" * 60)
    print("Panel G: GSEA Enrichment — Fibroblast + Epithelial")
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
    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, n_types * 2, 1,
                                  gridspec_kw={'height_ratios': [3, 1] * n_types})

    for idx, (cell_type, pre_res) in enumerate(gsea_results.items()):
        ax_top = axes[idx * 2]
        ax_bot = axes[idx * 2 + 1]
        label = CELL_TYPE_CONFIG[cell_type]["label"]
        draw_gsea_subplot(ax_top, ax_bot, pre_res, label)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, BASE_DIR / 'panel_g_gsea_4types')
    print(f"\nSaved: panel_g_gsea_4types.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
