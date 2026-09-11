#!/usr/bin/env python3
"""
Figure 5 panel G - violin plots of the top ligand activities from IL1B+
macrophages.

TOP_N, the mean-AUROC ranking, the NF-kB Hallmark membership test, the ligand
ordering and the seeded jitter are unchanged, and `np.random.seed(0)` stays
exactly where it was.

  printed panel  Figure 5 G       (PROVENANCE.csv - the directory is "05_C";
                                   do NOT read the directory as the letter)

MARK
    The earlier drawing used a canvas four times the printed size and set its
    smallest body type - the x axis label and both sets of tick labels - at
    6 * SCALE, so

        MARK = style.tick_pt() / (SMALL_PT * SCALE)

    The violin outline width and the strip-plot marker size are scaled by it.
    Tick widths and lengths and spine widths are not: those are style, and
    cnsplots sets them.
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
import panel_style_cns as style
import slots

BASE_DIR = Path(__file__).parent

SCALE = 4                       # the earlier canvas multiplier, for MARK only
SMALL_PT = 6.0                  # the earlier smallest body type, before * SCALE

PANEL_LETTER = "G"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(5, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(5, PANEL_LETTER)

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

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

    nfkb_genes = get_nfkb_genes()
    df = load_all_ligand_activities()

    mean_auroc = df.groupby('test_ligand')['auroc'].mean().sort_values(ascending=False)
    top_ligands = mean_auroc.head(TOP_N).index.tolist()
    plot_df = df[df['test_ligand'].isin(top_ligands)].copy()
    plot_df['test_ligand'] = pd.Categorical(plot_df['test_ligand'], categories=top_ligands, ordered=True)

    is_nfkb = {lig: lig in nfkb_genes for lig in top_ligands}
    palette = [COLOR_NFKB if is_nfkb.get(lig, False) else COLOR_OTHER for lig in top_ligands]

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)

    sns.violinplot(data=plot_df, y='test_ligand', x='auroc', ax=ax,
                   palette=palette, orient='h', cut=0, inner='box',
                   linewidth=1 * MARK)
    sns.stripplot(data=plot_df, y='test_ligand', x='auroc', ax=ax,
                  color='black', alpha=0.4, size=4 * MARK, jitter=True)

    ax.set_xlabel('Ligand Activity\n(AUC)')
    ax.set_ylabel('')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")
    style.save_panel(fig, BASE_DIR / 'macrophage_top_ligands_violin')
    print(f"  Saved: macrophage_top_ligands_violin.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
