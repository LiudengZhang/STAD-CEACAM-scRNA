#!/usr/bin/env python3
"""
S2 panel H - TCGA-STAD overall survival by BayesPrism-deconvolved epithelial
CEACAM5 and CEACAM6 expression (median split), drawn at the size it prints at.

  printed panel  Supplementary Figure S2 H   (PROVENANCE.csv; the submission-tree
                 directory was S3_D - the letter is looked up, not inferred)
  predecessor    submission-tree/03_Final_Panels/10_Supplementaries/
                 S2_CEACAM_Metaprogram_Validation/S2_D/create_S2_D_survival.py

Each box is titled "TCGA-STAD / <gene> (Epi. deconvolved)" as the printed
page sets it; the predecessor's figure-level suptitle is declared as a
removal in labels.py.

DRAWING READS A TABLE (cnsfig.cache): data/km_curves.csv holds, per gene and
group, lifelines' Kaplan-Meier survival function (timeline in months,
estimate), the group size and the log-rank P; the TCGA tables are read and
the estimator fitted only when the table is absent or with --recompute.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parents[1] / "_drivers"))
from paths import TCGA_BAYESPRISM_EPI, TCGA_CLINICAL          # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
import _driver_base as base                               # noqa: E402

FIG, PANEL = "S3_CEACAM_Metaprogram_Validation", "S3_H"
COLORS = {'High': '#B2182B', 'Low': '#2166AC'}
GENES = ['CEACAM5', 'CEACAM6']

# The printed box, millimetres: two boxes in the last row, beside G since
# 2026-09-16 (the author's fifth reading: G moves down into letter order,
# H gives it the room; 110 -> 95 mm).
W, H = 95.0, 48.0
MARGIN = dict(left=10.0, right=1.5, top=8.5, bottom=8.0)


def compute_curves():
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    epi = pd.read_csv(TCGA_BAYESPRISM_EPI, sep='\t', index_col=0)
    if 'CEACAM5' in epi.index:
        epi = epi.T
    clin = pd.read_csv(TCGA_CLINICAL, sep='\t')
    clin['os_time'] = clin['days_to_death'].fillna(clin['days_to_last_follow_up'])
    clin['os_event'] = (clin['vital_status'] == 'Dead').astype(int)
    clin = clin.dropna(subset=['os_time'])
    clin = clin[clin['os_time'] > 0].set_index('submitter_id')
    common = sorted(set(epi.index) & set(clin.index))
    rows = []
    for gene in GENES:
        expr = np.log2(epi.loc[common, gene].values + 1)
        high = expr >= np.median(expr)
        ids = {'High': [c for c, h in zip(common, high) if h],
               'Low': [c for c, h in zip(common, high) if not h]}
        t = {g: clin.loc[ids[g], 'os_time'].values / 30.44 for g in ids}
        e = {g: clin.loc[ids[g], 'os_event'].values for g in ids}
        p = logrank_test(t['High'], t['Low'], event_observed_A=e['High'],
                         event_observed_B=e['Low']).p_value
        for g in ('High', 'Low'):
            kmf = KaplanMeierFitter().fit(t[g], event_observed=e[g])
            sf = kmf.survival_function_
            for tl, val in zip(sf.index.to_numpy(dtype=float), sf.iloc[:, 0].to_numpy(dtype=float)):
                rows.append({'gene': gene, 'group': g, 'n': len(ids[g]),
                             'timeline': tl, 'survival': val, 'p_logrank': float(p)})
    return pd.DataFrame(rows)


def draw(tab):
    base.apply_style()
    fig, axes = style.subplots_mm(W, H, 1, 2)
    style.margins_mm(fig, **MARGIN)
    fig.subplots_adjust(wspace=0.35)
    for ax, gene in zip(axes, GENES):
        sub = tab[tab['gene'] == gene]
        for g in ('High', 'Low'):
            s = sub[sub['group'] == g]
            ax.plot(s['timeline'].to_numpy(), s['survival'].to_numpy(),
                    drawstyle='steps-post', color=COLORS[g], linewidth=style.RULE_PT,
                    label=f"{g} (n={int(s['n'].iloc[0])})")
        p = float(sub['p_logrank'].iloc[0])
        p_str = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        ax.text(0.95, 0.95, f'Log-rank\n{p_str}', transform=ax.transAxes, ha='right',
                va='top', fontsize=style.tick_pt())
        ax.set_title(f'TCGA-STAD\n{gene} (Epi. deconvolved)', pad=2)
        ax.set_xlabel('Time (months)')
        ax.set_ylabel('Overall Survival' if gene == GENES[0] else '')
        ax.set_ylim(0, 1.05)
        ax.legend(loc='lower left', frameon=False, handlelength=1.2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"ink outside the {W} x {H} mm canvas: {over}")
    return fig


def main():
    tab = cache.table(HERE, "km_curves", compute_curves)
    base.save(draw(tab), FIG, PANEL, "panel_S3_H")
    return 0


if __name__ == "__main__":
    sys.exit(main())
