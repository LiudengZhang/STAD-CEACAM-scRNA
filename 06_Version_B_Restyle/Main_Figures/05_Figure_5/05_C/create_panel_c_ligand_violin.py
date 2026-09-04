#!/usr/bin/env python3
"""
Figure 5 panel G, RESTYLED (Version B) - violin plots of the top ligand
activities from IL1B+ macrophages.

Version A is
`03_Final_Panels/05_Figure_5/05_C/create_panel_c_ligand_violin.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

TOP_N, the mean-AUROC ranking, the NF-kB Hallmark membership test, the ligand
ordering and the seeded jitter are Version A's, unchanged. `np.random.seed(0)`
stays exactly where Version A put it.

  printed panel  Figure 5 G        (PROVENANCE.csv - the directory is "05_C";
                 do NOT read the directory as the letter)
  printed rect   30.8 x 64.5 mm    (panel_rects.csv)
  Version B box  50.0 x 88.0 mm

    Fifteen ligand rows need about 5 mm each to be read at 7 pt, and the
    two-line x axis label needs about 8 mm below the axis. The published
    30.8 x 64.5 mm cannot carry either.

MARK
    Version A drew 16 x 33.2 cm at SCALE = 4 and set its smallest body type -
    the x axis label and both sets of tick labels - at `6 * SCALE`, so
    SMALL_PT = 6 and

        MARK = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 24 = 0.2917

    The violin outline width and the strip-plot marker size are scaled by it.
    Tick widths/lengths and spine widths are not: those are style, and cnsplots
    sets them.
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
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *
import panel_style_cns as style

BASE_DIR = Path(__file__).parent

SCALE = 4                       # Version A's canvas multiplier, for MARK only
SMALL_PT = 6.0                  # Version A's smallest body type, before * SCALE

PRINTED_MM = (30.8, 64.5)       # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 50.0, 88.0
MARGIN = dict(left=12.0, right=2.0, top=1.5, bottom=11.0)

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

    family = style.apply()
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

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, BASE_DIR / 'macrophage_top_ligands_violin')
    print(f"  Saved: macrophage_top_ligands_violin.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == "__main__":
    main()
