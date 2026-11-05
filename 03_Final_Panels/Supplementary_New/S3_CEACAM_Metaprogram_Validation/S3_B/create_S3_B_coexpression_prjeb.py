#!/usr/bin/env python3
"""Supplementary Figure S3 B: CEACAM5 against CEACAM6 in the BayesPrism-
deconvolved epithelial expression of PRJEB25780, coloured by response, with
Spearman rho and a least-squares line.

Reads TIGER_BAYESPRISM_EPI and TIGER_META; writes the panel through
_driver_base.save. Run: python create_S3_B_coexpression_prjeb.py"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "_drivers"),
                str(Path(__file__).resolve().parents[4] / "00_Config")]
from paths import TIGER_BAYESPRISM_EPI, TIGER_META        # noqa: E402
import panel_style_cns as style                           # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_B"
RESPONSE_COLORS = {'R': '#B2182B', 'NR': '#2166AC'}

# Panel size in mm.
W, H = 36.0, 44.0
MARGIN = dict(left=9.0, right=4.0, top=8.0, bottom=8.0)
DOT_MM = 1.4

def load():
    epi = pd.read_csv(TIGER_BAYESPRISM_EPI, sep='\t', index_col=0)
    if 'CEACAM5' in epi.index:
        epi = epi.T
    meta = pd.read_csv(TIGER_META, sep='\t')
    resp = meta.set_index('sample_id')['response_NR'].to_dict()
    resp = {k: ('R' if v == 'R' else 'NR') for k, v in resp.items() if pd.notna(v)}
    common = [s for s in epi.index if s in resp]
    df = pd.DataFrame({
        'CEACAM5': np.log2(epi.loc[common, 'CEACAM5'].values + 1),
        'CEACAM6': np.log2(epi.loc[common, 'CEACAM6'].values + 1),
        'Response': [resp[s] for s in common]}, index=common)
    return df[df['Response'].isin(['R', 'NR'])]

def draw(df):
    base.apply_style()
    fig, ax = style.subplots_mm(W, H)
    style.margins_mm(fig, **MARGIN)
    for response in ['NR', 'R']:
        sub = df[df['Response'] == response]
        ax.scatter(sub['CEACAM5'], sub['CEACAM6'], c=RESPONSE_COLORS[response],
                   s=(DOT_MM * style.PT_PER_MM) ** 2, alpha=0.7, edgecolors='white',
                   linewidths=style.EDGE_PT, label=f'{response} (n={len(sub)})', zorder=3)
    r, p = stats.spearmanr(df['CEACAM5'], df['CEACAM6'])
    slope, intercept, _, _, _ = stats.linregress(df['CEACAM5'], df['CEACAM6'])
    xs = np.linspace(df['CEACAM5'].min(), df['CEACAM5'].max(), 100)
    ax.plot(xs, slope * xs + intercept, 'k--', linewidth=style.RULE_PT, alpha=0.7)
    p_str = 'P < 0.001' if p < 0.001 else style.p_label(p)
    ax.text(0.05, 0.95, f'ρ = {r:.2f}, {p_str}\nn = {len(df)}',
            transform=ax.transAxes, fontsize=style.tick_pt(), va='top')
    ax.set_xlabel('Epi. CEACAM5 (log2+1)')
    ax.set_ylabel('Epi. CEACAM6 (log2+1)')
    ax.set_title('CEACAM5/6 Co-expression\n(PRJEB25780)')
    ax.legend(loc='lower right', frameon=False, handletextpad=0.3)
    base.despine(ax)
    return fig

if __name__ == "__main__":
    base.save(draw(load()), FIG, PANEL, "panel_S3_B")

