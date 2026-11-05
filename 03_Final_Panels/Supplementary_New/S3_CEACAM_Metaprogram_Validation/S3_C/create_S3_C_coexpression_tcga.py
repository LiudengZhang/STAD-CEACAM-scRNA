#!/usr/bin/env python3
"""Supplementary Figure S3 C: CEACAM5 against CEACAM6 in the BayesPrism-
deconvolved epithelial expression of TCGA-STAD (407 tumour-only samples), with
Spearman rho and a least-squares line.

Reads TCGA_BAYESPRISM_EPI and TCGA_BULK_TUMOR_ONLY; writes the panel through
_driver_base.save. Run: python create_S3_C_coexpression_tcga.py"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import TCGA_BAYESPRISM_EPI, TCGA_BULK_TUMOR_ONLY  # noqa: E402
import panel_style_cns as style                           # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_C"

# Panel size in mm.
W, H = 36.0, 44.0
MARGIN = dict(left=9.0, right=4.0, top=8.0, bottom=8.0)
DOT_MM = 1.0

def load():
    """Deconvolved epithelial CEACAM5/6 (log2+1) of the 407 tumour-only samples.

    """
    epi = pd.read_csv(TCGA_BAYESPRISM_EPI, sep='\t', index_col=0)
    tumor_only = pd.read_csv(TCGA_BULK_TUMOR_ONLY, sep='\t', index_col=0,
                             nrows=1).columns
    keep = sorted(set(epi.index) & set(tumor_only))
    if len(keep) != 407:
        raise SystemExit(f"expected the 407 tumour-only TCGA-STAD samples the "
                         f"Methods describe, got {len(keep)}")
    epi = epi.loc[keep]
    return pd.DataFrame({'CEACAM5': np.log2(epi['CEACAM5'].values + 1),
                         'CEACAM6': np.log2(epi['CEACAM6'].values + 1)})

def draw(df):
    base.apply_style()
    fig, ax = style.subplots_mm(W, H)
    style.margins_mm(fig, **MARGIN)
    ax.scatter(df['CEACAM5'], df['CEACAM6'], c='#666666',
               s=(DOT_MM * style.PT_PER_MM) ** 2, alpha=0.5, edgecolors='white',
               linewidths=style.EDGE_PT, zorder=3)
    r, p = stats.spearmanr(df['CEACAM5'], df['CEACAM6'])
    slope, intercept, _, _, _ = stats.linregress(df['CEACAM5'], df['CEACAM6'])
    xs = np.linspace(df['CEACAM5'].min(), df['CEACAM5'].max(), 100)
    ax.plot(xs, slope * xs + intercept, 'k--', linewidth=style.RULE_PT, alpha=0.7)
    p_str = 'P < 0.001' if p < 0.001 else style.p_label(p)
    ax.text(0.05, 0.95, f'ρ = {r:.2f}, {p_str}\nn = {len(df)}',
            transform=ax.transAxes, fontsize=style.tick_pt(), va='top')
    ax.set_xlabel('CEACAM5 (log2+1)')
    ax.set_ylabel('CEACAM6 (log2+1)')
    ax.set_title('CEACAM5/6 Co-expression\n(TCGA-STAD)')
    base.despine(ax)
    return fig

if __name__ == "__main__":
    base.save(draw(load()), FIG, PANEL, "panel_S3_C")

