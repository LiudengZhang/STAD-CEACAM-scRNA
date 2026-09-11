#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Step 3: Survival analysis using BayesPrism-deconvolved cell fractions.

Stratifies TCGA-STAD patients by:
  1. Fibroblast fraction (proxy for C6_Fib)
  2. Monocytes/Macrophages fraction (proxy for C3_Mac)
  3. Combined Fib + MoMac fraction
  4. Fibroblast-specific IL6 expression (deconvolved)
  5. MoMac-specific IL1B expression (deconvolved)

Each produces a KM survival plot with log-rank test.
"""
import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

SCALE = 4

# Migrated 2026-09-03 on the author's ruling. Input and output used to be one
# directory, submission-tree/temp-workspace/02162026_BayesPrism,
# a temp workspace that can be cleared at any time - and the whole TCGA
# BayesPrism chain went with it when it was. The tables now come from the
# pipeline's declared home, paths.py TCGA_BAYESPRISM_DIR, which already held two
# of them. Not a byte moved: every md5 in ../VERDICT.md was re-checked on both
# sides of the copy, and the temp workspace was copied from, never emptied.
# Outputs are written beside this script rather than back into the deposited
# directory, so a re-run cannot overwrite what it is being compared against.
INPUT_DIR = '/path/to/submission-tree/02_Preparation_for_Panels/BayesPrism_TCGA'
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

COLOR_HIGH = '#B2182B'
COLOR_LOW = '#2166AC'


def load_data():
    """Load BayesPrism fractions and clinical data."""
    fractions = pd.read_csv(f'{INPUT_DIR}/tcga_bayesprism_fractions.tsv', sep='\t', index_col=0)
    clinical = pd.read_csv(f'{INPUT_DIR}/tcga_clinical.tsv', sep='\t')

    # Merge
    fractions['submitter_id'] = fractions.index
    merged = fractions.merge(clinical, on='submitter_id', how='inner')

    # Survival time
    merged['time'] = merged.apply(
        lambda r: r['days_to_death'] if pd.notna(r['days_to_death']) and r['vital_status'] == 'Dead'
        else r['days_to_last_follow_up'], axis=1)
    merged['event'] = (merged['vital_status'] == 'Dead').astype(int)
    merged = merged.dropna(subset=['time'])
    merged = merged[merged['time'] > 0]
    merged['time_months'] = merged['time'] / 30.44

    print(f"Samples with survival data: {len(merged)}")
    return merged, fractions


def km_plot(merged, score_col, label_high, label_low, title, filename):
    """Create KM survival plot stratified by median of score_col."""
    median_val = merged[score_col].median()
    merged = merged.copy()
    merged['group'] = merged[score_col].apply(lambda x: label_high if x >= median_val else label_low)

    n_high = (merged['group'] == label_high).sum()
    n_low = (merged['group'] == label_low).sum()

    fig, ax = plt.subplots(figsize=(7 * SCALE / 2.54, 6 * SCALE / 2.54))
    kmf = KaplanMeierFitter()

    colors = {label_high: COLOR_HIGH, label_low: COLOR_LOW}
    for grp in [label_high, label_low]:
        mask = merged['group'] == grp
        n = mask.sum()
        kmf.fit(merged.loc[mask, 'time_months'], merged.loc[mask, 'event'],
                label=f'{grp} (n={n})')
        kmf.plot_survival_function(ax=ax, color=colors[grp], linewidth=0.8,
                                    ci_show=True, ci_alpha=0.15)

    # Log-rank
    high = merged[merged['group'] == label_high]
    low = merged[merged['group'] == label_low]
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
    print(f"  {title}: {p_text} (n_high={n_high}, n_low={n_low})")
    print(f"  Saved: {out_path}")
    return lr.p_value


def main():
    print("=" * 60)
    print("Step 3: Survival analysis by deconvolved fractions")
    print("=" * 60)

    merged, fractions = load_data()

    # Print available cell types
    ct_cols = [c for c in fractions.columns if c != 'submitter_id']
    print(f"\nCell types: {ct_cols}")
    for ct in ct_cols:
        if ct in merged.columns:
            print(f"  {ct}: mean={merged[ct].mean():.4f} ({merged[ct].mean()*100:.1f}%)")

    results = {}

    # 1. Fibroblast fraction
    if 'Fibroblast' in merged.columns:
        p = km_plot(merged, 'Fibroblast',
                    'Fib-high', 'Fib-low',
                    'TCGA-STAD: survival by\nfibroblast fraction (BayesPrism)',
                    'survival_fibroblast_fraction.png')
        results['Fibroblast fraction'] = p

    # 2. MoMac fraction
    momac_col = None
    for col in merged.columns:
        if 'Monocytes' in col or 'Macrophages' in col:
            momac_col = col
            break
    if momac_col:
        p = km_plot(merged, momac_col,
                    'MoMac-high', 'MoMac-low',
                    'TCGA-STAD: survival by\nMoMac fraction (BayesPrism)',
                    'survival_momac_fraction.png')
        results['MoMac fraction'] = p

    # 3. Combined Fib + MoMac
    if 'Fibroblast' in merged.columns and momac_col:
        merged['fib_momac_combined'] = merged['Fibroblast'] + merged[momac_col]
        p = km_plot(merged, 'fib_momac_combined',
                    'Fib+MoMac high', 'Fib+MoMac low',
                    'TCGA-STAD: survival by\ncombined Fib + MoMac fraction',
                    'survival_fib_momac_combined.png')
        results['Fib+MoMac combined'] = p

    # 4. Fibroblast-specific IL6 (if available)
    try:
        fib_expr = pd.read_csv(f'{INPUT_DIR}/tcga_bayesprism_fibroblast_expression.tsv',
                                sep='\t', index_col=0)
        if 'IL6' in fib_expr.columns:
            merged_fib = merged.copy()
            fib_il6 = fib_expr['IL6']
            fib_il6.name = 'fib_IL6'
            # Repaired 2026-09-03 on the author's ruling. As deposited this was
            # merged_fib.join(fib_il6, how='inner'), which joins on merged_fib's
            # OWN index - and merged_fib gets a fresh index from the .merge()
            # in load_data(), not the sample ids fib_il6 is indexed by. Nothing
            # matched, the frame came back empty, and lifelines raised
            # "Values must be numeric" - so sections 4 and 5 have never run.
            # submitter_id is already a column here, so the join is given it.
            # Index alignment only: same table, same gene, same median split.
            merged_fib = merged_fib.join(fib_il6, on='submitter_id', how='inner')
            merged_fib['fib_IL6_log'] = np.log2(merged_fib['fib_IL6'] + 1)
            p = km_plot(merged_fib, 'fib_IL6_log',
                        'Fib-IL6 high', 'Fib-IL6 low',
                        'TCGA-STAD: survival by\nfibroblast-specific $\\it{IL6}$ (BayesPrism)',
                        'survival_fib_il6.png')
            results['Fib-specific IL6'] = p
    except FileNotFoundError:
        print("  Fibroblast expression file not found, skipping IL6 analysis")

    # 5. MoMac-specific IL1B (if available)
    try:
        momac_expr = pd.read_csv(f'{INPUT_DIR}/tcga_bayesprism_momac_expression.tsv',
                                  sep='\t', index_col=0)
        if 'IL1B' in momac_expr.columns:
            merged_mac = merged.copy()
            mac_il1b = momac_expr['IL1B']
            mac_il1b.name = 'mac_IL1B'
            # Same repair as section 4 above, same reason.
            merged_mac = merged_mac.join(mac_il1b, on='submitter_id', how='inner')
            merged_mac['mac_IL1B_log'] = np.log2(merged_mac['mac_IL1B'] + 1)
            p = km_plot(merged_mac, 'mac_IL1B_log',
                        'Mac-IL1B high', 'Mac-IL1B low',
                        'TCGA-STAD: survival by\nMoMac-specific $\\it{IL1B}$ (BayesPrism)',
                        'survival_momac_il1b.png')
            results['MoMac-specific IL1B'] = p
    except FileNotFoundError:
        print("  MoMac expression file not found, skipping IL1B analysis")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, pval in results.items():
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
        print(f"  {name:30s}  P = {pval:.4f}  {sig}")


if __name__ == '__main__':
    main()
