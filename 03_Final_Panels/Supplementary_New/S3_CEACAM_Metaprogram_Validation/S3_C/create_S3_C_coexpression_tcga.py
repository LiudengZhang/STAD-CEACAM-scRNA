#!/usr/bin/env python3
"""
S2 panel C - CEACAM5 against CEACAM6 in the BayesPrism-deconvolved epithelial
expression of TCGA-STAD, drawn at the size it prints at.

  printed panel  Supplementary Figure S2 C   (PROVENANCE.csv; the submission-tree
                 directory was S3_F - the letter is looked up, not inferred)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/
                 S2_CEACAM_Metaprogram_Validation/S2_F/create_S2_F_coexpression_tcga.py

The drawing reads the two small deconvolution tables directly, as the
predecessor did; nothing is recomputed beyond the Spearman rho and the
regression line the predecessor also computed from them. The white plate
behind the statistics is gone (the author's ruling on plates, 2026-09-15).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import TCGA_BAYESPRISM_EPI                    # noqa: E402
import panel_style_cns as style                           # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_C"

# The printed box, millimetres: one of four in row 2.
W, H = 36.0, 44.0
MARGIN = dict(left=9.0, right=4.0, top=8.0, bottom=8.0)
DOT_MM = 1.0


def load():
    epi = pd.read_csv(TCGA_BAYESPRISM_EPI, sep='\t', index_col=0)
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
    # No mathtext: a superscript exponent or a subscript 2 renders at 0.7 of
    # the size and falls under the 6 pt floor (sweep_pages reads 4.2 pt), so
    # P is printed as the paper prints a P below 0.001 and log2 is set plain.
    # Both declared in labels.py RENAMES_S1_S6.
    p_str = 'P < 0.001' if p < 0.001 else style.p_label(p)
    ax.text(0.05, 0.95, f'ρ = {r:.2f}, {p_str}\nn = {len(df)}',
            transform=ax.transAxes, fontsize=style.tick_pt(), va='top')
    ax.set_xlabel('CEACAM5 (log2+1)')
    ax.set_ylabel('CEACAM6 (log2+1)')
    ax.set_title('CEACAM5/6 Co-expression\n(TCGA-STAD)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"ink outside the {W} x {H} mm canvas: {over}")
    return fig


def main():
    base.save(draw(load()), FIG, PANEL, "panel_S3_C")
    return 0


if __name__ == "__main__":
    sys.exit(main())
