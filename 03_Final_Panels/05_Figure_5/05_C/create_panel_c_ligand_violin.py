#!/usr/bin/env python3
"""
Panel C: Violin plots of top ligand activities from IL1B+ Macrophages.
4× scaling method for Nature Cancer.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import gseapy as gp
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *
from shared.figure_config import use_panel_style

BASE_DIR = Path(__file__).parent

# 4× scaling
DPI = 300
SCALE = 4
CM_TO_INCH = 1 / 2.54
PANEL_WIDTH_CM = 4.0 * SCALE
PANEL_HEIGHT_CM = 8.3 * SCALE

TOP_N = 15

COLOR_NFKB = '#d62728'
COLOR_OTHER = '#7f7f7f'


def _read_gmt(path):
    """The pinned Hallmark sets, in the shape gseapy.get_library returns."""
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def get_nfkb_genes():
    try:
        hallmark = _read_gmt(HALLMARK_GMT)
        return set(hallmark['HALLMARK_TNFA_SIGNALING_VIA_NFKB'])
    except Exception:
        return {'IL1B', 'TNF', 'IL6', 'CSF3', 'IL1A', 'LTA', 'LTB'}


def load_all_ligand_activities():
    print("Loading ligand activities...")
    # Use pre-aggregated raw CSV
    raw_csv = LIGAND_ACTIVITIES_RAW
    if raw_csv.exists():
        combined = pd.read_csv(raw_csv)
        print(f"  Loaded from raw CSV: {len(combined)} records, "
              f"{combined['test_ligand'].nunique()} ligands, {combined['receiver'].nunique()} receivers")
        return combined

    # Fallback: load from individual receiver directories
    all_activities = []
    pre_dir = NICHENET_DIR / "pre_responder"
    for receiver_dir in pre_dir.iterdir():
        if not receiver_dir.is_dir():
            continue
        csv_file = receiver_dir / "ligand_activities.csv"
        if csv_file.exists():
            try:
                df = pd.read_csv(csv_file)
                df['receiver'] = receiver_dir.name
                all_activities.append(df)
            except Exception:
                continue
    if not all_activities:
        raise ValueError("No ligand activity files found!")
    combined = pd.concat(all_activities, ignore_index=True)
    print(f"  {len(combined)} records, {combined['test_ligand'].nunique()} ligands, {combined['receiver'].nunique()} receivers")
    return combined


def main():
    # The overlaid points are placed with random jitter. Left unseeded it made
    # this panel the only kind in the figure that could not reproduce itself:
    # two consecutive runs of the unchanged script gave three different SVGs
    # (baseline, run 1 and run 2 all differed). The jitter is decoration - no
    # statistic depends on it - but a panel that redraws differently every time
    # cannot be checked, so it is pinned here.
    np.random.seed(0)

    use_panel_style(font_pt=7)

    nfkb_genes = get_nfkb_genes()
    df = load_all_ligand_activities()

    mean_auroc = df.groupby('test_ligand')['auroc'].mean().sort_values(ascending=False)
    top_ligands = mean_auroc.head(TOP_N).index.tolist()
    plot_df = df[df['test_ligand'].isin(top_ligands)].copy()
    plot_df['test_ligand'] = pd.Categorical(plot_df['test_ligand'], categories=top_ligands, ordered=True)

    is_nfkb = {lig: lig in nfkb_genes for lig in top_ligands}
    palette = [COLOR_NFKB if is_nfkb.get(lig, False) else COLOR_OTHER for lig in top_ligands]

    fig, ax = plt.subplots(figsize=(PANEL_WIDTH_CM * CM_TO_INCH, PANEL_HEIGHT_CM * CM_TO_INCH))

    sns.violinplot(data=plot_df, y='test_ligand', x='auroc', ax=ax,
                   palette=palette, orient='h', cut=0, inner='box',
                   linewidth=1)
    sns.stripplot(data=plot_df, y='test_ligand', x='auroc', ax=ax,
                  color='black', alpha=0.4, size=4, jitter=True)

    ax.set_xlabel('Ligand Activity\n(AUC)', fontsize=6 * SCALE)
    ax.set_ylabel('', fontsize=6 * SCALE)
    ax.tick_params(axis='both', labelsize=6 * SCALE, width=1.0, length=4)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)

    plt.tight_layout()

    output = BASE_DIR / 'macrophage_top_ligands_violin.png'
    plt.savefig(output, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.svg'), bbox_inches='tight', facecolor='white')
    plt.savefig(output.with_suffix('.pdf'), format='pdf', dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output}")


if __name__ == "__main__":
    main()
