#!/usr/bin/env python3
"""Figure 5G: NicheNet ligand activity (AUROC) of the 15 top-ranked ligands
(by mean AUROC across receivers) from IL1B+ macrophages, as violins with the
individual values overlaid. Ligands in the Hallmark TNF-alpha/NF-kB set are
red.

Input:  LIGAND_ACTIVITIES_RAW, HALLMARK_GMT (00_Config/paths.py)
Output: macrophage_top_ligands_violin.{svg,pdf,png} in this directory
Run:    python create_panel_c_ligand_violin.py
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *
from cnsfig import Panel, style  # noqa: E402
from cnsfig.layout import pin_frame_mm


SCALE = 4                       # marker-size scaling constants
SMALL_PT = 6.0

PANEL = Panel(5, "G", __file__, small_pt=SMALL_PT)

TOP_N = 15

COLOR_NFKB = '#d62728'
COLOR_OTHER = '#7f7f7f'


def _read_gmt(path):
    """Read a GMT file into {set name: [genes]}."""
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
    """Read the aggregated ligand-activity table; stop if it is missing.
    """
    if not LIGAND_ACTIVITIES_RAW.exists():
        raise SystemExit(f"missing {LIGAND_ACTIVITIES_RAW}")
    return pd.read_csv(LIGAND_ACTIVITIES_RAW)


def main():
    # Seed the strip-plot jitter so the drawing is reproducible.
    np.random.seed(0)

    PANEL.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)

    nfkb_genes = get_nfkb_genes()
    df = load_all_ligand_activities()

    mean_auroc = df.groupby('test_ligand')['auroc'].mean().sort_values(ascending=False)
    top_ligands = mean_auroc.head(TOP_N).index.tolist()
    plot_df = df[df['test_ligand'].isin(top_ligands)].copy()
    plot_df['test_ligand'] = pd.Categorical(plot_df['test_ligand'], categories=top_ligands, ordered=True)

    is_nfkb = {lig: lig in nfkb_genes for lig in top_ligands}
    palette = [COLOR_NFKB if is_nfkb.get(lig, False) else COLOR_OTHER for lig in top_ligands]

    fig, ax = style.subplots_mm(PANEL.w_mm, PANEL.h_mm)

    sns.violinplot(data=plot_df, y='test_ligand', x='auroc', ax=ax,
                   palette=palette, orient='h', cut=0, inner='box',
                   linewidth=style.EDGE_PT)
    # seaborn draws the inner box at 3x `linewidth`; reset it.
    for ln in ax.lines:
        ln.set_linewidth(style.RULE_PT)
    sns.stripplot(data=plot_df, y='test_ligand', x='auroc', ax=ax,
                  color='black', alpha=0.4, size=4 * MARK, jitter=True)

    ax.set_xlabel('Ligand Activity\n(AUC)')
    ax.set_ylabel('')
    ax.set_yticklabels([t.get_text() for t in ax.get_yticklabels()],
                       style='italic')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    PANEL.fit(fig)
    # Align the frame top with the row; pad x so the violins clear the y axis.
    pin_frame_mm(fig, ax, top_mm=4.3, bottom_mm=ax.get_position().y0 * PANEL.h_mm)
    _lo, _hi = ax.get_xlim()
    ax.set_xlim(_lo - 0.08 * (_hi - _lo), _hi)
    PANEL.finish(fig, 'macrophage_top_ligands_violin')


if __name__ == "__main__":
    main()
