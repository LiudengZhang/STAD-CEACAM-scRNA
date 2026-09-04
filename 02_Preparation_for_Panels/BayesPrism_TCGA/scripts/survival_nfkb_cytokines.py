#!/usr/bin/env python3
# Paths below refer to the upstream Round_4 processing pipeline, which is
# not part of this release. This script is included as a record of how the
# input was produced; it is not called by _run_all_panels.sh.
"""Survival analysis: NF-κB-activating cytokine scores in TCGA-STAD.

Score 1: IL1A + IL1B + TNF + IL6 (4-gene, core NF-κB activators)
Score 2: IL1A + IL1B + TNF + IL6 + OSM + LTA (6-gene, extended)

Each scored as mean z-score of log2(FPKM+1), median-split KM."""
import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import zscore
import warnings
warnings.filterwarnings('ignore')

SCALE = 4
# TCGA_BASE reached the clinical table through
# Round_4/.../98_External/Bulk/02_TCGA_STAD, a symlink into
# Round_5/01_Raw_Inputs/04_External - and 04_External was renamed 02_External,
# so that route has been broken since. This is fault 3 of ../VERDICT.md, and it
# is repointed here exactly as step1_prepare_tcga_bulk_pinned.py already
# repoints it: the same file, at the name it now has.
TCGA_BASE = '/path/to/Project_4_05232025/Round_5/01_Raw_Inputs/02_External/Bulk/TCGA_STAD'

# Migrated 2026-09-03 on the author's ruling. Input and output used to be one
# directory, Round_5/98_Temp_workspace/02162026_BayesPrism, a temp workspace
# that can be cleared at any time - and the whole TCGA BayesPrism chain went
# with it when it was. The table now comes from the pipeline's declared home,
# paths.py TCGA_BAYESPRISM_DIR, which already held it. Not a byte moved: the
# md5 was checked on both sides of the copy, and the temp workspace was copied
# from, never emptied. Outputs are written beside this script rather than back
# into the deposited directory, so a re-run cannot overwrite what it is being
# compared against.
INPUT_DIR = '/path/to/Project_4_05232025/Round_5/02_Preparation_for_Panels/BayesPrism_TCGA'
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

GENES = {
    'IL1A': 'ENSG00000115008.6',
    'IL1B': 'ENSG00000125538.12',
    'TNF':  'ENSG00000232810.4',
    'IL6':  'ENSG00000136244.12',
    'OSM':  'ENSG00000099985.4',
    'LTA':  'ENSG00000226979.9',
}

COLOR_HIGH = '#B2182B'
COLOR_LOW = '#2166AC'


def run_km(merged, score_col, label_high, label_low, title, filename):
    """KM plot with median split."""
    median_val = merged[score_col].median()
    df = merged.copy()
    df['group'] = df[score_col].apply(lambda x: label_high if x >= median_val else label_low)

    fig, ax = plt.subplots(figsize=(7 * SCALE / 2.54, 6 * SCALE / 2.54))
    kmf = KaplanMeierFitter()
    colors = {label_high: COLOR_HIGH, label_low: COLOR_LOW}

    for grp in [label_high, label_low]:
        mask = df['group'] == grp
        n = mask.sum()
        kmf.fit(df.loc[mask, 'time_months'], df.loc[mask, 'event'],
                label=f'{grp} (n={n})')
        kmf.plot_survival_function(ax=ax, color=colors[grp], linewidth=0.8,
                                    ci_show=True, ci_alpha=0.15)

    high = df[df['group'] == label_high]
    low = df[df['group'] == label_low]
    lr = logrank_test(high['time_months'], low['time_months'], high['event'], low['event'])
    p_text = f'P = {lr.p_value:.1e}' if lr.p_value < 0.001 else f'P = {lr.p_value:.3f}'

    ax.text(0.95, 0.95, f'Log-rank {p_text}',
            transform=ax.transAxes, ha='right', va='top', fontsize=4 * SCALE,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, linewidth=0.3))

    ax.set_xlabel('Time (months)', fontsize=5 * SCALE)
    ax.set_ylabel('Overall survival', fontsize=5 * SCALE)
    ax.set_title(title, fontsize=5.5 * SCALE, pad=5)
    ax.legend(fontsize=4 * SCALE, loc='lower left', framealpha=0.8)
    ax.tick_params(width=1.0, length=4, labelsize=4 * SCALE)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)
    ax.set_ylim(0, 1.05)

    plt.tight_layout()
    out_path = f'{OUTPUT_DIR}/{filename}'
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  {p_text} | Saved: {out_path}")
    return lr.p_value


def main():
    print("=" * 60)
    print("NF-κB-activating cytokine survival analysis")
    print("=" * 60)

    # Load clean tumor-only FPKM
    print("\nLoading TCGA-STAD FPKM (tumor-only, n=407)...")
    fpkm = pd.read_csv(f'{INPUT_DIR}/tcga_bulk_counts_tumor_only.tsv',
                        sep='\t', index_col=0)

    # Extract genes
    ens_ids = list(GENES.values())
    names = list(GENES.keys())
    expr = fpkm.loc[ens_ids].T
    expr.columns = names
    print(f"Samples: {len(expr)}")
    for g in names:
        print(f"  {g}: mean FPKM = {expr[g].mean():.2f}")

    # Log2 transform and z-score
    expr_log2 = np.log2(expr + 1)
    expr_z = expr_log2.apply(zscore, axis=0)

    # Compute scores
    core4 = ['IL1A', 'IL1B', 'TNF', 'IL6']
    all6 = ['IL1A', 'IL1B', 'TNF', 'IL6', 'OSM', 'LTA']
    expr_z['score_4gene'] = expr_z[core4].mean(axis=1)
    expr_z['score_6gene'] = expr_z[all6].mean(axis=1)

    print(f"\n4-gene score (IL1A/IL1B/TNF/IL6): mean={expr_z['score_4gene'].mean():.3f}, std={expr_z['score_4gene'].std():.3f}")
    print(f"6-gene score (+OSM/LTA):           mean={expr_z['score_6gene'].mean():.3f}, std={expr_z['score_6gene'].std():.3f}")

    # Load clinical
    clinical = pd.read_csv(f'{TCGA_BASE}/02_Raw_Data/Clinical/TCGA_STAD_clinical_data.tsv', sep='\t')

    # Merge
    score_df = expr_z[['score_4gene', 'score_6gene']].copy()
    score_df['submitter_id'] = score_df.index
    merged = score_df.merge(clinical[['submitter_id', 'vital_status', 'days_to_death', 'days_to_last_follow_up']],
                             on='submitter_id', how='inner')

    merged['time'] = merged.apply(
        lambda r: r['days_to_death'] if pd.notna(r['days_to_death']) and r['vital_status'] == 'Dead'
        else r['days_to_last_follow_up'], axis=1)
    merged['event'] = (merged['vital_status'] == 'Dead').astype(int)
    merged = merged.dropna(subset=['time'])
    merged = merged[merged['time'] > 0]
    merged['time_months'] = merged['time'] / 30.44
    print(f"\nValid survival data: {len(merged)} samples")

    # KM plots
    print("\n--- 4-gene score (IL1A, IL1B, TNF, IL6) ---")
    p4 = run_km(merged, 'score_4gene',
                'Cytokine-high', 'Cytokine-low',
                'TCGA-STAD: survival by NF-\u03baB-activating\ncytokine score ($\\it{IL1A}$/$\\it{IL1B}$/$\\it{TNF}$/$\\it{IL6}$)',
                'survival_nfkb_cytokine_4gene.png')

    print("\n--- 6-gene score (IL1A, IL1B, TNF, IL6, OSM, LTA) ---")
    p6 = run_km(merged, 'score_6gene',
                'Cytokine-high', 'Cytokine-low',
                'TCGA-STAD: survival by NF-\u03baB-activating\ncytokine score (6-gene)',
                'survival_nfkb_cytokine_6gene.png')

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    sig4 = "***" if p4 < 0.001 else "**" if p4 < 0.01 else "*" if p4 < 0.05 else "ns"
    sig6 = "***" if p6 < 0.001 else "**" if p6 < 0.01 else "*" if p6 < 0.05 else "ns"
    print(f"  4-gene (IL1A/IL1B/TNF/IL6):  P = {p4:.4f}  {sig4}")
    print(f"  6-gene (+OSM/LTA):           P = {p6:.4f}  {sig6}")


if __name__ == '__main__':
    main()
