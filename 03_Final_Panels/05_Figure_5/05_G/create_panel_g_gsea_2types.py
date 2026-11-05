#!/usr/bin/env python3
"""Figure 5I: GSEA running enrichment score for the Hallmark TNF-alpha/NF-kB
set in fibroblasts and epithelial cells, post-treatment non-responders vs
responders.

Genes are ranked by logFC * -log10(p) from a t-test (t-test_overestim_var)
on up to 5000 subsampled cells, then run through gseapy prerank (1000
permutations, seed 42). NES, P and FDR are printed beside each curve.

Curves and statistics are cached in data/es_curves.csv and data/gsea_stats.csv;
the h5ads are read only when a table is missing or --recompute is passed.

Input:  FIBROBLAST_H5AD, EPITHELIAL_H5AD, HALLMARK_GMT (00_Config/paths.py)
Output: panel_g_gsea_4types.{svg,pdf,png} in this directory
Run:    python create_panel_g_gsea_2types.py [--recompute]
"""
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import gseapy as gp
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *
from cnsfig import Panel, style  # noqa: E402
from cnsfig import cache
import pandas as pd

#: Width of the stats column right of each plot, and the gap before it (mm).
STATS_W_MM = 11.6
STATS_GAP_MM = 1.0

PANEL = Panel(5, "I", __file__, small_pt=5.0)

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
    gsea_outdir = PANEL.dir / f"gsea_{cell_type}"
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
    """t-test DEGs and gseapy prerank for each cell type (computed once)."""
    if _STATE:
        return _STATE["results"]
    gsea_results = {}
    for cell_type, config in CELL_TYPE_CONFIG.items():
        if not config["path"].exists():
            raise SystemExit(f"  Not found: {config['path']}")
        adata = sc.read_h5ad(config["path"])
        adata.obs_names_make_unique()

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
    """Draw the running ES on one axes and return its stats block's text."""
    nes, pval, fdr = stats['nes'], stats['pval'], stats['fdr']
    x = np.arange(len(RES))
    color = '#C62828'

    ax_top.plot(x, RES, color=color, linewidth=style.RULE_PT)
    ax_top.axhline(y=0, color='gray', linestyle='--', linewidth=style.RULE_PT)
    ax_top.fill_between(x, 0, RES, where=(np.array(RES) >= 0), color='#EF9A9A', alpha=0.3)
    ax_top.fill_between(x, 0, RES, where=(np.array(RES) < 0), color='#90CAF9', alpha=0.3)

    peak_idx = int(stats['peak_idx'])
    ax_top.axvline(x=peak_idx, color='red', linestyle=':', linewidth=style.RULE_PT, alpha=0.7)

    ax_top.set_xlim(0, len(RES))
    ax_top.set_title(label)

    # P below the permutation resolution prints as "< 0.001".
    def fmt(name, v):
        return f"{name} < 0.001" if v < 0.001 else f"{name} = {v:.3f}"

    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    return f"NES = {nes:.2f}\n{fmt('P', pval)}\n{fmt('FDR', fdr)}"


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
    PANEL.apply()

    stats = cache.table(PANEL.dir, "gsea_stats", compute_stats).set_index("cell_type")
    curves = cache.table(PANEL.dir, "es_curves", compute_curves)
    gsea_results = {ct: curves.loc[curves["cell_type"] == ct, "RES"].to_numpy()
                    for ct in CELL_TYPE_CONFIG if ct in set(curves["cell_type"])}
    if not gsea_results:
        raise SystemExit("no GSEA results")

    # x axes are not shared: the two ranked lists differ in length.
    n_types = len(gsea_results)
    fig, axes = style.subplots_mm(PANEL.w_mm, PANEL.h_mm, n_types, 1)
    axes = np.atleast_1d(axes)
    fig.subplots_adjust(hspace=0.85)

    blocks = []
    for idx, (cell_type, RES) in enumerate(gsea_results.items()):
        label = CELL_TYPE_CONFIG[cell_type]["label"]
        blocks.append((axes[idx], draw_gsea_subplot(axes[idx], RES, stats.loc[cell_type], label), RES))

    for one in axes:
        one.set_ylabel('ES')
    PANEL.fit(fig)
    for ax, txt, RES in blocks:
        place_stats_block(fig, ax, txt)
    PANEL.finish(fig, 'panel_g_gsea_4types')


if __name__ == "__main__":
    main()
