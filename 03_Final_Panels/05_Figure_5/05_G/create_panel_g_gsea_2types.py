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

DRAWING READS TABLES, since 2026-09-15 (cnsfig.cache): data/es_curves.csv
(cell_type, rank, RES) and data/gsea_stats.csv (cell_type, nes, pval, fdr,
n_hits, peak_idx). The h5ads, the t-test and gseapy run only when a table is
absent or --recompute is passed.

THE STATS BLOCK STANDS BESIDE THE CURVE, since 2026-09-15: the framed box of
the earlier drawing covered the tail of both curves (the author's third
reading), and on a 9 mm axes a three-line block at 6 pt is two thirds of the
height - it cannot stand inside the axes without flattening the curve to a
third (measured: the y limit would have to reach 3.5 for a peak of 0.9). So
the plotting box gives up STATS_W_MM at its right and the block stands there,
unframed, top-aligned with the box. Same three strings.
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
from cnsfig import cache
import pandas as pd

BASE_DIR = Path(__file__).parent
#: The column at the right of each plotting box that holds the three-line
#: stats block ("FDR < 0.001" at 6 pt is 10.6 mm), and the paper between.
STATS_W_MM = 11.6
STATS_GAP_MM = 1.0

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


_STATE = {}


def _run_all():
    """The computing half: t-test DEGs and gseapy prerank per cell type."""
    if _STATE:
        return _STATE["results"]
    gsea_results = {}
    for cell_type, config in CELL_TYPE_CONFIG.items():
        print(f"\n=== {cell_type} ===")
        if not config["path"].exists():
            raise SystemExit(f"  Not found: {config['path']}")
        adata = sc.read_h5ad(config["path"])
        adata.obs_names_make_unique()
        print(f"  {adata.n_obs} cells")

        deg_df = run_deg_ttest(adata, COMPARISON)
        if deg_df is None:
            raise SystemExit(f"{cell_type}: no DEGs")
        pre_res = run_gsea(deg_df, cell_type)
        if pre_res is None or PATHWAY_NAME not in pre_res.results:
            raise SystemExit(f"{cell_type}: {PATHWAY_NAME!r} not in the GSEA results")
        gsea_results[cell_type] = pre_res
    _STATE["results"] = gsea_results
    return gsea_results


def compute_stats():
    rows = []
    for cell_type, pre_res in _run_all().items():
        t = pre_res.results[PATHWAY_NAME]
        RES = np.asarray(t['RES'], dtype=float)
        rows.append({"cell_type": cell_type, "nes": float(t['nes']),
                     "pval": float(t['pval']), "fdr": float(t['fdr']),
                     "n_hits": len(t['hits']), "peak_idx": int(np.argmax(np.abs(RES)))})
    return pd.DataFrame(rows)


def compute_curves():
    out = []
    for cell_type, pre_res in _run_all().items():
        RES = np.asarray(pre_res.results[PATHWAY_NAME]['RES'], dtype=float)
        out.append(pd.DataFrame({"cell_type": cell_type,
                                 "rank": np.arange(len(RES)), "RES": RES}))
    return pd.concat(out, ignore_index=True)


def draw_gsea_subplot(ax_top, RES, stats, label):
    """Draw the running ES and its stats block on one axes."""
    nes, pval, fdr = stats['nes'], stats['pval'], stats['fdr']
    x = np.arange(len(RES))
    color = '#C62828'

    # Top: running ES
    ax_top.plot(x, RES, color=color, linewidth=style.RULE_PT)
    ax_top.axhline(y=0, color='gray', linestyle='--', linewidth=style.RULE_PT)
    ax_top.fill_between(x, 0, RES, where=(np.array(RES) >= 0), color='#EF9A9A', alpha=0.3)
    ax_top.fill_between(x, 0, RES, where=(np.array(RES) < 0), color='#90CAF9', alpha=0.3)

    peak_idx = int(stats['peak_idx'])
    ax_top.axvline(x=peak_idx, color='red', linestyle=':', linewidth=style.RULE_PT, alpha=0.7)

    ax_top.set_xlim(0, len(RES))
    ax_top.set_title(label)

    # The printed panel reads "p < 0.001", not "p = 0.000". A permutation P
    # value below the resolution of the test is bounded, never zero, and three
    # decimal places on their own turn the bound into a false exact value.
    def fmt(name, v):
        return f"{name} < 0.001" if v < 0.001 else f"{name} = {v:.3f}"

    stats_text = f"NES = {nes:.2f}\n{fmt('p', pval)}\n{fmt('FDR', fdr)}"
    # Beside the box, not over the curve: drawn by place_stats_block once
    # fit_margins has set the box, so the margins are fitted to the curve's
    # own furniture and the block takes a column of its own.
    txt = stats_text

    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)

    print(f"    {label}: NES={nes:.2f}, p={pval:.3f}, FDR={fdr:.3f}, "
          f"{int(stats['n_hits'])} leading-edge hits")
    return txt


def place_stats_block(fig, ax, stats_text):
    """Narrow the plotting box by the stats column and stand the block in
    it, STATS_GAP_MM right of the box's right spine, top-aligned."""
    box = ax.get_position()
    w_mm = fig.get_figwidth() * 25.4
    new_w = box.width - (STATS_W_MM + STATS_GAP_MM) / w_mm
    if new_w * w_mm < 20.0:
        raise RuntimeError(f"the plotting box would be {new_w * w_mm:.1f} mm "
                           f"wide beside the stats column")
    ax.set_position([box.x0, box.y0, new_w, box.height])
    ax_w_mm = new_w * w_mm
    return ax.text(1.0 + STATS_GAP_MM / ax_w_mm, 1.0, stats_text,
                   transform=ax.transAxes, fontsize=style.tick_pt(),
                   va='top', ha='left', linespacing=1.2)


def main():
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    print("=" * 60)
    print("Panel I: GSEA Enrichment — Fibroblast + Epithelial")
    print("  (t-test DEGs, min_size=5)")
    print("=" * 60)

    stats = cache.table(BASE_DIR, "gsea_stats", compute_stats).set_index("cell_type")
    curves = cache.table(BASE_DIR, "es_curves", compute_curves)
    gsea_results = {ct: curves.loc[curves["cell_type"] == ct, "RES"].to_numpy()
                    for ct in CELL_TYPE_CONFIG if ct in set(curves["cell_type"])}
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

    blocks = []
    for idx, (cell_type, RES) in enumerate(gsea_results.items()):
        label = CELL_TYPE_CONFIG[cell_type]["label"]
        blocks.append((axes[idx], draw_gsea_subplot(axes[idx], RES, stats.loc[cell_type], label), RES))

    # Each curve names its own axis (2026-09-14, evening): a label shared
    # across the two was clipped to "Enrichment Scor" on the page and hung
    # in the gap here. 'ES' is the standard short form of the enrichment
    # score, declared in labels.py; the legend spells it out.
    for one in axes:
        one.set_ylabel('ES')
    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    for ax, txt, RES in blocks:
        place_stats_block(fig, ax, txt)
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
