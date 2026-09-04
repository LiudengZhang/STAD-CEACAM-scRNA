#!/usr/bin/env python3
"""
Figure 3 panel D, RESTYLED (Version B) - CEACAM5/6 vs CD8+ Tex fraction (1x2),
ALL 32 stomach samples, 4-group colouring, no external dataset.

NOTE THE PANEL LETTER. This directory is `03_B` but PROVENANCE.csv says it
holds printed panel **D**. CLAUDE.md rule 2: never infer a panel letter from a
directory name.

Version A is
`03_Revised_Panels/Main_Figures/03_Figure_3/03_B/create_ceacam_tex_scatter.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every value read, every filter, every statistic and every string is Version A's.
The drawing code is the same code - the Tex state selection
(`C6_CD8_Tex_PDCD1`), the per-sample percentage, the merge, the 4-group
assignment, the Spearman test and the P-value formatting are untouched, and
expression is read exactly as Version A reads it.

  printed panel  Figure 3 D     (PROVENANCE.csv; NOT inferred from "03_B")
  printed rect   50.0 x 22.9 mm   (panel_rects.csv)
  Version B box  123.0 x 57.0 mm

MARK
    Version A drew a 32.0 x 16.0 cm canvas (SCALE = 4). Its smallest body type
    is the shared figure legend at `fontsize=4.5 * SCALE` - smaller than the
    5 * SCALE tick labels - so SMALL_PT = 4.5 and, by PANEL_SPEC.md,

        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 18 = 0.389
        AREA = MARK ** 2                    = 0.151

    Marker area 30*SCALE, marker edge 0.3*SCALE, the dashed regression line
    0.8*SCALE, the spines 0.5*SCALE and the tick width/length 0.5*SCALE /
    3*SCALE are each multiplied by it, so all keep their Version A size
    *relative to the type*.

THE LEGEND MOVED, AND WHY
    Version A anchored the shared legend at `bbox_to_anchor=(1.0, 0.5)` - i.e.
    entirely *outside* the canvas - and relied on `bbox_inches='tight'` at save
    time to enlarge the saved image until the legend fitted. Version B must not
    pass `bbox_inches` (it is what breaks 1:1; see panel_style_cns.save_panel),
    so a legend left outside the canvas would simply be cut off. 18 mm of the
    box is therefore reserved on the right and the legend is anchored at the
    inside edge of that reserve. Same handles, same labels, same order, same
    `markerscale=0.8`; only where it sits on the paper changed.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import TCD8_H5AD, EPITHELIAL_RAW_COUNTS_H5AD
import panel_style_cns as style  # noqa: E402

# Version A's canvas multiplier. Version B draws 1:1, so SCALE survives only to
# reproduce the exact numbers Version A set for its non-type lengths.
SCALE = 4
SMALL_PT = 4.5                       # Version A's smallest body type (legend)

# See MARK above. `style.tick_pt()` reads cnsplots' own setting, so the factor
# is derived, never a literal.
MARK = style.tick_pt() / (SMALL_PT * SCALE)   # 0.389, length multiplier
AREA = MARK ** 2                              # 0.151, area multiplier

PRINTED_MM = (50.0, 22.9)            # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 123.0, 57.0
# Millimetres of paper. right = 18 is the strip the shared legend sits in.
MARGIN = dict(left=14.0, right=18.0, top=9.0, bottom=9.0, wspace=1.0 / 3.0)
LEGEND_X = (PANEL_W_MM - MARGIN['right'] + 2.0) / PANEL_W_MM

EPI_RAW = EPITHELIAL_RAW_COUNTS_H5AD
OUTPUT_DIR = Path(__file__).parent

# Standard 4-group palette (blue=R, warm=NR; lighter=Pre, darker=Post)
COLOR_MAP = {
    'Pre-R':  '#74b9ff',
    'Pre-NR': '#e17055',
    'Post-R': '#0984e3',
    'Post-NR':'#d63031',
    'Other':  '#999999',
}

DRAW_ORDER = ['Other', 'Pre-R', 'Pre-NR', 'Post-R', 'Post-NR']


def load_primary():
    """Load ALL 32 stomach samples: epithelial CEACAM + CD8 Tex fraction."""
    # --- Epithelial: CEACAM expression ---
    print("[Primary] Loading epithelial raw counts...")
    epi = sc.read_h5ad(EPI_RAW)

    # Filter: stomach only (all treatment phases)
    stomach_mask = epi.obs['Sample site'].astype(str).str.lower().str.contains('stomach')
    epi = epi[stomach_mask].copy()

    # Normalize raw counts
    sc.pp.normalize_total(epi, target_sum=1e4)
    sc.pp.log1p(epi)

    epi_df = pd.DataFrame({
        'sample': epi.obs['sample'].values,
        'CEACAM5': epi[:, 'CEACAM5'].X.toarray().flatten() if hasattr(epi[:, 'CEACAM5'].X, 'toarray') else epi[:, 'CEACAM5'].X.flatten(),
        'CEACAM6': epi[:, 'CEACAM6'].X.toarray().flatten() if hasattr(epi[:, 'CEACAM6'].X, 'toarray') else epi[:, 'CEACAM6'].X.flatten(),
        'treatment_phase': epi.obs['Treatment phase'].values,
        'pre_group': epi.obs['stomach_pre_grouping'].values,
        'post_group': epi.obs['stomach_post_grouping'].values,
    })
    epi_sample = epi_df.groupby('sample').agg({
        'CEACAM5': 'mean', 'CEACAM6': 'mean',
        'treatment_phase': 'first',
        'pre_group': 'first',
        'post_group': 'first',
    }).reset_index()
    stomach_samples = set(epi_sample['sample'].values)
    print(f"  {len(epi_sample)} stomach samples")
    del epi

    # --- CD8: Tex fraction ---
    print("[Primary] Loading CD8+ T cells...")
    cd8 = sc.read_h5ad(TCD8_H5AD)
    cd8_df = pd.DataFrame({
        'sample': cd8.obs['sample'].values,
        'minor_cell_state': cd8.obs['minor_cell_state'].values,
    })
    del cd8

    # Filter to stomach samples
    cd8_df = cd8_df[cd8_df['sample'].isin(stomach_samples)]

    cd8_list = []
    for sample, grp in cd8_df.groupby('sample'):
        n_total = len(grp)
        n_tex = (grp['minor_cell_state'] == 'C6_CD8_Tex_PDCD1').sum()
        cd8_list.append({
            'sample': sample,
            'tex_fraction': (n_tex / n_total * 100) if n_total > 0 else 0,
        })
    cd8_sample = pd.DataFrame(cd8_list)
    print(f"  CD8 data for {len(cd8_sample)} samples")

    # Merge
    merged = epi_sample.merge(cd8_sample, on='sample', how='inner')

    # 4-group assignment
    def assign_4group(row):
        phase = str(row['treatment_phase'])
        if phase == 'Pre':
            grp = str(row['pre_group'])
            if grp == 'Responsed': return 'Pre-R'
            if grp == 'No-response': return 'Pre-NR'
        elif phase == 'Post':
            grp = str(row['post_group'])
            if grp == 'Responsed': return 'Post-R'
            if grp == 'No-response': return 'Post-NR'
        return 'Other'

    merged['group'] = merged.apply(assign_4group, axis=1)

    counts = merged['group'].value_counts()
    print(f"  Merged: {len(merged)} samples: {counts.to_dict()}")
    return merged


def main():
    family = style.apply()
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.3f}")

    data = load_primary()

    # 1x2 layout
    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, 1, 2)

    for col_idx, ceacam in enumerate(['CEACAM5', 'CEACAM6']):
        ax = axes[col_idx]
        x = data[ceacam].values.astype(float)
        y = data['tex_fraction'].values.astype(float)
        groups = data['group'].values

        r_val, p_val = stats.spearmanr(x, y)

        for grp in DRAW_ORDER:
            mask = groups == grp
            if mask.sum() > 0:
                ax.scatter(x[mask], y[mask], c=COLOR_MAP[grp], s=30*SCALE*AREA,
                           alpha=0.85, edgecolors='white', linewidths=0.3*SCALE*MARK,
                           label=grp, zorder=3)

        # Regression line
        valid = np.isfinite(x) & np.isfinite(y)
        xv, yv = x[valid], y[valid]
        if len(xv) >= 3 and xv.std() > 0:
            slope, intercept = np.polyfit(xv, yv, 1)
            x_line = np.linspace(xv.min(), xv.max(), 100)
            ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=0.8*SCALE*MARK, alpha=0.6, zorder=2)

        # Stats text — 1 sig digit (floor), scientific for very small P
        import math
        _e = math.floor(math.log10(p_val)); _c = int(p_val / 10**_e)
        if _e >= -3:
            p_str = f'P = {_c * 10**_e:.{-_e}f}'
        else:
            p_str = f'P = {_c}' + r'$\times 10^{' + str(_e) + r'}$'

        ax.set_title(f'Primary Cohort (scRNA-seq)\nρ = {r_val:.2f}, {p_str}',
                     linespacing=1.4)

        ax.set_xlabel(f'{ceacam} (mean log expr.)')
        ax.set_ylabel('CD8$^+$ Tex fraction (%)')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for spine in ['bottom', 'left']:
            ax.spines[spine].set_linewidth(0.5*SCALE*MARK)
        ax.tick_params(axis='both', width=0.5*SCALE*MARK, length=3*SCALE*MARK)

        ax.set_box_aspect(1)

    # Single shared legend on the far right — now inside the canvas, see
    # THE LEGEND MOVED above.
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center left',
               bbox_to_anchor=(LEGEND_X, 0.5),
               framealpha=0, edgecolor='none', markerscale=0.8)

    style.margins_mm(fig, **MARGIN)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, OUTPUT_DIR / 'ceacam_tex_scatter')
    print(f"\nSaved: {OUTPUT_DIR / 'ceacam_tex_scatter'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
