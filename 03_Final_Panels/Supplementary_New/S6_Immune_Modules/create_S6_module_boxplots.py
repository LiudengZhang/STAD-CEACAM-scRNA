#!/usr/bin/env python3
"""
S5 panels A-E - the five immune-module proportions, pre-treatment R against
NR and post-treatment R against NR, drawn at the size they print at.

  printed panels  Supplementary Figure S6 A-E   (PROVENANCE.csv; one row per
                  panel, each pointing at this script)
  predecessor     submission-tree/03_Final_Panels/10_Supplementaries/S5_Immune_Modules/
                  create_S5_module_boxplots.py  (five 5.5 x 5 cm canvases at 4x)

One script draws the five panels, as before, each into its own directory
S6_<letter>/ so the assembler finds one SVG per panel. Since 2026-09-16 the
boxes and brackets are cnsfig.boxes' (one box style through the paper, no
individual points - the author's ruling of that day) and the P labels go
through panel_style_cns.p_text_kw like every other pairwise bracket: a star
below 0.05, the two-decimal value at or above it. The predecessor and the
printed page set "ns" there; the legend says which.

DRAWING READS A TABLE (cnsfig.cache): data/module_samples.csv holds each
sample's proportion in each module with its group, and data/tests.csv the two
Mann-Whitney P values per module; the module-proportion and clinical tables
are joined only when absent or with --recompute.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(HERE.parent / "_drivers"))
from paths import FIG4_MODULE_PROPORTIONS, FIG4_CLINICAL_METADATA   # noqa: E402
import panel_style_cns as style                           # noqa: E402
from cnsfig import cache                                  # noqa: E402
from cnsfig import boxes                                  # noqa: E402
from cnsfig.boxes import box_xlim                         # noqa: E402
import _driver_base as base                               # noqa: E402

FIG = "S6_Immune_Modules"
COLORS = {'Pre-R': '#bde0fe', 'Pre-NR': '#a2d2ff', 'Post-R': '#ffcfd2', 'Post-NR': '#f1c0e8'}
MODULE_NAMES = {1: 'IM-T/NK/DC', 2: 'IM-MoMac', 3: 'IM-Mixed', 4: 'IM-Neutrophil',
                5: 'IM-B/Plasma'}
MODULE_LETTERS = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E'}
GROUPS = ['Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']
POSITIONS = [0, 1, 2.5, 3.5]

# The printed box, millimetres: three to a row.
W, H = 53.0, 48.0
MARGIN = dict(left=9.5, right=1.5, top=4.5, bottom=12.5)


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
    # ONE BOX, ONE BRACKET, NO POINTS (cnsfig.boxes, 2026-09-16, the author's
    # fifth reading: "box plots should look alike throughout; don't show every
    # point"). Until then: 0.7-alpha faces, the samples jittered on top (rng
    # 42) and a local star ladder that printed "ns" where every other pairwise
    # bracket in the paper prints the two-decimal value (p_label).
    bxp = boxes.draw_boxes(ax, data, POSITIONS, [COLORS[g] for g in GROUPS], width=0.6)
    allv = np.concatenate([d for d in data if len(d)])
    y_max, y_rng = allv.max(), allv.max() - allv.min()
    y_rng = y_rng or 0.1
    labels = []
    for (x0, x1), p in (((0, 1), p_pre), ((2.5, 3.5), p_post)):
        _, txt, _ = boxes.bracket(fig, ax, x0, x1, y_max, y_rng, p, kind="pair",
                                  lift=0.08, arm=0.02)
        labels.append(txt)
    ax.set_xticks(POSITIONS)
    ax.set_xticklabels([f"{g}\n(n={len(d)})" for g, d in zip(GROUPS, data)])
    for x, t in ((0.5, 'Pre-treatment'), (3.0, 'Post-treatment')):
        # Regular weight since 2026-09-16: only panel letters and UMAP cluster
        # labels are bold (RULES.md, line weight and type weight).
        ax.text(x, -0.30, t, ha='center', va='top', fontsize=style.tick_pt(),
                transform=ax.get_xaxis_transform())
    ax.set_xlim(*box_xlim(POSITIONS, 0.6))
    ax.set_ylabel('Module Proportion')
    ax.set_title(MODULE_NAMES[m], pad=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for txt in labels:
        boxes.ylim_above(ax, txt)
    boxes.assert_no_points(ax, bxp)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(f"S5 {MODULE_LETTERS[m]}: ink outside the {W} x {H} mm canvas: {over}")
    return fig


def main():
    base.apply_style()
    samples = cache.table(HERE, "module_samples", compute_samples)
    tests = cache.table(HERE, "tests", compute_tests)
    for m in range(1, 6):
        letter = MODULE_LETTERS[m]
        base.save(draw(m, samples, tests), FIG, f"S6_{letter}", f"panel_S6_{letter}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
