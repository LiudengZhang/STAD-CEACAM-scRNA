#!/usr/bin/env python3
"""
S5: Immune Module Proportion Boxplots (M1-M5)
Each panel shows TWO side-by-side pairs:
  LEFT: Pre-treatment (Pre-R vs Pre-NR)
  RIGHT: Post-treatment (Post-R vs Post-NR)
4-group colors: Pre-R #bde0fe, Pre-NR #a2d2ff, Post-R #ffcfd2, Post-NR #f1c0e8
Mann-Whitney U test for each pair, sample-level.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FIG4_MODULE_PROPORTIONS, FIG4_CLINICAL_METADATA

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

SCALE = 4
DPI = 300
CM = 1 / 2.54

# 4-group colors (poster scheme)
COLORS = {
    'Pre-R': '#bde0fe',
    'Pre-NR': '#a2d2ff',
    'Post-R': '#ffcfd2',
    'Post-NR': '#f1c0e8',
}

MODULE_NAMES = {
    1: 'IM-T/NK/DC',
    2: 'IM-MoMac',
    3: 'IM-Mixed',
    4: 'IM-Neutrophil',
    5: 'IM-B/Plasma',
}

MODULE_LETTERS = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E'}

OUT_DIR = Path(__file__).parent


def create_module_panel(module_num, module_props, clinical_data):
    """Create a single module boxplot panel."""
    letter = MODULE_LETTERS[module_num]
    name = MODULE_NAMES[module_num]
    col = f'Module_{module_num}'
    panel_dir = OUT_DIR / f'S5_{letter}'
    panel_dir.mkdir(exist_ok=True)

    print(f"\n  Module {module_num}: {name}")

    merged = module_props.join(clinical_data, how='inner')

    # Pre-treatment groups
    pre = merged[
        (merged['Sample site'] == 'Stomach') &
        (merged['Treatment phase'] == 'Pre') &
        (merged['stomach_pre_grouping'].isin(['No-response', 'Responsed']))
    ].copy()
    pre['group'] = pre['stomach_pre_grouping'].map({'Responsed': 'Pre-R', 'No-response': 'Pre-NR'})

    # Post-treatment groups
    post = merged[
        (merged['Sample site'] == 'Stomach') &
        (merged['Treatment phase'] == 'Post') &
        (merged['stomach_post_grouping'].isin(['No-response', 'Responsed']))
    ].copy()
    post['group'] = post['stomach_post_grouping'].map({'Responsed': 'Post-R', 'No-response': 'Post-NR'})

    pre_r = pre[pre['group'] == 'Pre-R'][col].values
    pre_nr = pre[pre['group'] == 'Pre-NR'][col].values
    post_r = post[post['group'] == 'Post-R'][col].values
    post_nr = post[post['group'] == 'Post-NR'][col].values

    print(f"    Pre-R: n={len(pre_r)}, Pre-NR: n={len(pre_nr)}")
    print(f"    Post-R: n={len(post_r)}, Post-NR: n={len(post_nr)}")

    # Mann-Whitney tests
    if len(pre_r) >= 2 and len(pre_nr) >= 2:
        _, p_pre = mannwhitneyu(pre_r, pre_nr, alternative='two-sided')
    else:
        p_pre = np.nan

    if len(post_r) >= 2 and len(post_nr) >= 2:
        _, p_post = mannwhitneyu(post_r, post_nr, alternative='two-sided')
    else:
        p_post = np.nan

    print(f"    Pre P={p_pre:.4f}" if not np.isnan(p_pre) else "    Pre P=N/A")
    print(f"    Post P={p_post:.4f}" if not np.isnan(p_post) else "    Post P=N/A")

    # Create figure
    fig, ax = plt.subplots(figsize=(5.5 * SCALE * CM, 5.0 * SCALE * CM))

    groups_data = [pre_r, pre_nr, post_r, post_nr]
    group_names = ['Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']
    positions = [0, 1, 2.5, 3.5]

    bp = ax.boxplot(groups_data, positions=positions, widths=0.6,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.5 * SCALE),
                    whiskerprops=dict(linewidth=1.0 * SCALE),
                    capprops=dict(linewidth=1.0 * SCALE),
                    boxprops=dict(linewidth=1.0 * SCALE))

    for i, (box, gname) in enumerate(zip(bp['boxes'], group_names)):
        box.set_facecolor(COLORS[gname])
        box.set_alpha(0.7)

    # Jitter points
    rng = np.random.default_rng(42)
    for i, (data, pos, gname) in enumerate(zip(groups_data, positions, group_names)):
        if len(data) > 0:
            jitter = rng.uniform(-0.1, 0.1, len(data))
            ax.scatter([pos] * len(data) + jitter, data, c=COLORS[gname],
                       s=20 * SCALE, edgecolors='black', linewidths=0.3 * SCALE,
                       alpha=0.85, zorder=3)

    # Significance brackets
    all_vals = np.concatenate([d for d in groups_data if len(d) > 0])
    if len(all_vals) > 0:
        y_max = np.max(all_vals)
        y_range = np.max(all_vals) - np.min(all_vals)
        if y_range == 0:
            y_range = 0.1

        # Pre bracket
        bh_pre = y_max + 0.08 * y_range
        ax.plot([0, 0, 1, 1], [bh_pre - 0.02 * y_range, bh_pre, bh_pre, bh_pre - 0.02 * y_range],
                color='black', linewidth=0.8 * SCALE)
        p_str_pre = _format_pval(p_pre)
        ax.text(0.5, bh_pre + 0.02 * y_range, p_str_pre,
                ha='center', fontsize=5 * SCALE)

        # Post bracket
        bh_post = y_max + 0.08 * y_range
        ax.plot([2.5, 2.5, 3.5, 3.5],
                [bh_post - 0.02 * y_range, bh_post, bh_post, bh_post - 0.02 * y_range],
                color='black', linewidth=0.8 * SCALE)
        p_str_post = _format_pval(p_post)
        ax.text(3.0, bh_post + 0.02 * y_range, p_str_post,
                ha='center', fontsize=5 * SCALE)

        ax.set_ylim(ax.get_ylim()[0], bh_pre + 0.20 * y_range)

    # Labels
    ax.set_xticks(positions)
    labels = [f"{gname}\n(n={len(d)})" for gname, d in zip(group_names, groups_data)]
    ax.set_xticklabels(labels, fontsize=4.5 * SCALE)

    # Add Pre/Post group labels
    ax.text(0.5, -0.18, 'Pre-treatment', ha='center', fontsize=5 * SCALE,
            fontweight='bold', transform=ax.get_xaxis_transform())
    ax.text(3.0, -0.18, 'Post-treatment', ha='center', fontsize=5 * SCALE,
            fontweight='bold', transform=ax.get_xaxis_transform())

    ax.set_ylabel('Module Proportion', fontsize=5 * SCALE)
    ax.set_title(name, fontsize=6 * SCALE, fontweight='normal')

    ax.tick_params(labelsize=4.5 * SCALE, width=1.0 * SCALE, length=4 * SCALE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(1.0 * SCALE)

    plt.tight_layout()

    stem = f'panel_S5_{letter}'
    for ext in ['svg', 'png']:
        kw = {'bbox_inches': 'tight', 'facecolor': 'white'}
        if ext == 'png':
            kw['dpi'] = DPI
        fig.savefig(panel_dir / f'{stem}.{ext}', **kw)
    print(f"    Saved: {stem}")
    plt.close()


def _format_pval(p):
    if np.isnan(p):
        return 'N/A'
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def main():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

    print("=" * 60)
    print("S5: Immune Module Proportion Boxplots")
    print("=" * 60)

    # Load data
    module_props = pd.read_csv(FIG4_MODULE_PROPORTIONS, index_col=0)
    clinical_data = pd.read_csv(FIG4_CLINICAL_METADATA)
    if 'Sample ID' in clinical_data.columns:
        clinical_data = clinical_data.set_index('Sample ID')

    print(f"Module proportions: {module_props.shape}")
    print(f"Clinical metadata: {clinical_data.shape}")

    for mod_num in range(1, 6):
        create_module_panel(mod_num, module_props, clinical_data)

    print("\nAll 5 module panels complete!")


if __name__ == '__main__':
    main()
