#!/usr/bin/env python3
"""Supplementary Figure S8 A-E: proportions of the five immune modules in stomach
samples, R against NR before and after treatment (two-sided Mann-Whitney U,
exact P on each bracket). One script draws all five panels into S8_A ... S8_E.

Reads data/module_samples.csv and data/tests.csv; the module-proportion and
clinical tables are joined only if those are missing or with --recompute.
Run: python create_S8_module_boxplots.py [--recompute]

"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / "_drivers"),
                str(Path(__file__).resolve().parents[3] / "00_Config")]
from paths import FIG4_MODULE_PROPORTIONS, FIG4_CLINICAL_METADATA   # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
from cnsfig import boxes                                  # noqa: E402
from cnsfig.boxes import box_xlim                         # noqa: E402
import _driver_base as base                               # noqa: E402

FIG = "S8_Immune_Modules"
from shared.figure_config import GROUP_COLORS  # noqa: E402
COLORS = dict(GROUP_COLORS)
MODULE_NAMES = {1: 'IM-T/NK/DC', 2: 'IM-MoMac', 3: 'IM-Mixed', 4: 'IM-Neutrophil',
                5: 'IM-B/Plasma'}
MODULE_LETTERS = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E'}
GROUPS = ['Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']
POSITIONS = [0, 1, 2.5, 3.5]

# Panel size in mm.
W, H = 53.0, 44.5
MARGIN = dict(left=9.5, right=1.5, top=4.5, bottom=9.0)

def compute_samples():
    props = pd.read_csv(FIG4_MODULE_PROPORTIONS, index_col=0)
    clin = pd.read_csv(FIG4_CLINICAL_METADATA)
    if 'Sample ID' in clin.columns:
        clin = clin.set_index('Sample ID')
    merged = props.join(clin, how='inner')
    stomach = merged['Sample site'] == 'Stomach'
    pre = merged[stomach & (merged['Treatment phase'] == 'Pre')
                 & merged['stomach_pre_grouping'].isin(['No-response', 'Responsed'])].copy()
    pre['group'] = pre['stomach_pre_grouping'].map({'Responsed': 'Pre-R', 'No-response': 'Pre-NR'})
    post = merged[stomach & (merged['Treatment phase'] == 'Post')
                  & merged['stomach_post_grouping'].isin(['No-response', 'Responsed'])].copy()
    post['group'] = post['stomach_post_grouping'].map({'Responsed': 'Post-R', 'No-response': 'Post-NR'})
    both = pd.concat([pre, post])
    rows = []
    for m in range(1, 6):
        for g in GROUPS:
            for v in both[both['group'] == g][f'Module_{m}'].to_numpy(dtype=float):
                rows.append({'module': m, 'group': g, 'value': v})
    return pd.DataFrame(rows)

def compute_tests():
    from scipy.stats import mannwhitneyu
    tab = cache.table(HERE, "module_samples", compute_samples)
    rows = []
    for m in range(1, 6):
        sub = tab[tab['module'] == m]
        for phase, a, b in (('pre', 'Pre-R', 'Pre-NR'), ('post', 'Post-R', 'Post-NR')):
            x = sub[sub['group'] == a]['value'].to_numpy()
            y = sub[sub['group'] == b]['value'].to_numpy()
            p = mannwhitneyu(x, y, alternative='two-sided')[1] if len(x) >= 2 and len(y) >= 2 else np.nan
            rows.append({'module': m, 'phase': phase, 'p_two_sided': float(p)})
    return pd.DataFrame(rows)

def draw(m, samples, tests):
    fig, ax = style.subplots_mm(W, H)
    style.margins_mm(fig, **MARGIN)
    sub = samples[samples['module'] == m]
    data = [sub[sub['group'] == g]['value'].to_numpy() for g in GROUPS]
    p_pre = float(tests[(tests['module'] == m) & (tests['phase'] == 'pre')]['p_two_sided'].iloc[0])
    p_post = float(tests[(tests['module'] == m) & (tests['phase'] == 'post')]['p_two_sided'].iloc[0])
    bxp = boxes.draw_boxes(ax, data, POSITIONS, [COLORS[g] for g in GROUPS], width=0.6)
    allv = np.concatenate([d for d in data if len(d)])
    y_max, y_rng = allv.max(), allv.max() - allv.min()
    y_rng = y_rng or 0.1
    labels = []
    for (x0, x1), p in (((0, 1), p_pre), ((2.5, 3.5), p_post)):
        _, txt, _ = boxes.bracket(fig, ax, x0, x1, y_max, y_rng, p, kind="pair",
                                  lift=0.08, arm=0.02, exact=True)
        labels.append(txt)
    ax.set_xticks(POSITIONS)
    ax.set_xticklabels(list(GROUPS))
    ax.set_xlim(*box_xlim(POSITIONS, 0.6))
    ax.set_ylabel('Module Proportion')
    ax.set_title(MODULE_NAMES[m], pad=3)
    base.despine(ax)
    for txt in labels:
        boxes.ylim_above(ax, txt)
    boxes.assert_no_points(ax, bxp)
    return fig

if __name__ == "__main__":
    base.apply_style()
    samples = cache.table(HERE, "module_samples", compute_samples)
    tests = cache.table(HERE, "tests", compute_tests)
    for m in range(1, 6):
        letter = MODULE_LETTERS[m]
        base.save(draw(m, samples, tests), FIG, f"S8_{letter}", f"panel_S8_{letter}")
