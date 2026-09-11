#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
CEACAM vs Immune Infiltration - Mixed Effects Model Analysis

Controls for total epithelial proportion (tumor density) to identify
immune cell recruitment patterns specific to CEACAM-high vs CEACAM-low cancer cells.

Model: Immune ~ CEACAM_ratio + Total_Epi + (1|Sample)
- CEACAM_ratio = CEACAM-high / (CEACAM-high + CEACAM-low)
- Total_Epi = CEACAM-high + CEACAM-low (tumor density proxy)

Author: Generated for Round 4 Analysis
Date: 2026-01-04
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'


def load_proportion_data(data_dir):
    """Load all sample proportion CSVs."""
    print("Loading proportion data...")
    csv_files = sorted(Path(data_dir).glob('*_proportions_14types.csv'))

    all_data = []
    for csv_path in csv_files:
        sample_name = csv_path.stem.replace('_proportions_14types', '')
        df = pd.read_csv(csv_path, index_col=0)
        df['sample'] = sample_name
        df['spot_id'] = df.index
        all_data.append(df)
        print(f"  {sample_name}: {len(df)} spots")

    df_all = pd.concat(all_data, ignore_index=True)
    print(f"Total spots: {len(df_all):,}")
    return df_all


def main():
    print("="*70)
    print("CEACAM vs Immune - Mixed Effects Model (Controlling for Tumor Density)")
    print("="*70)

    # Paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'ceacam_immune_mixed_model'
    output_dir.mkdir(parents=True, exist_ok=True)

    data_dir = script_dir.parent.parent / '03_Extract_Results' / 'results' / 'deconvolution_matrices'

    # Load data
    df_all = load_proportion_data(data_dir)

    # Calculate derived variables
    df_all['Total_Epi'] = df_all['Epi CEACAM-high'] + df_all['Epi CEACAM-low']

    # CEACAM ratio (avoid division by zero)
    # Only calculate for spots with sufficient epithelial content
    epi_threshold = 0.05  # Minimum 5% epithelial to be considered
    df_epi = df_all[df_all['Total_Epi'] >= epi_threshold].copy()
    df_epi['CEACAM_ratio'] = df_epi['Epi CEACAM-high'] / df_epi['Total_Epi']

    print(f"\nFiltered to spots with ≥{epi_threshold*100:.0f}% epithelial: {len(df_epi):,} spots")
    print(f"CEACAM_ratio range: {df_epi['CEACAM_ratio'].min():.3f} - {df_epi['CEACAM_ratio'].max():.3f}")
    print(f"Total_Epi range: {df_epi['Total_Epi'].min():.3f} - {df_epi['Total_Epi'].max():.3f}")

    # Define immune cell types
    immune_cells = [
        'B cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells',
        'DC cells', 'Mast cells', 'Neutrophils', 'Plasma cells',
        'Monocytes/Macrophages'
    ]

    # ========== MIXED EFFECTS MODELS ==========
    print("\n" + "="*70)
    print("Fitting Mixed Effects Models...")
    print("Model: Immune ~ CEACAM_ratio + Total_Epi + (1|Sample)")
    print("="*70)

    results = []

    for immune_cell in immune_cells:
        print(f"\nProcessing: {immune_cell}")

        # Safe column name for formula
        safe_col = immune_cell.replace(' ', '_').replace('/', '_').replace('+', 'pos')
        df_epi[safe_col] = df_epi[immune_cell]

        # Fit mixed effects model
        formula = f"{safe_col} ~ CEACAM_ratio + Total_Epi"

        try:
            model = smf.mixedlm(formula, df_epi, groups=df_epi['sample'])
            fit = model.fit(method='lbfgs', maxiter=1000)

            # Extract CEACAM_ratio coefficient (the effect after controlling for Total_Epi)
            coef = fit.params['CEACAM_ratio']
            se = fit.bse['CEACAM_ratio']
            pval = fit.pvalues['CEACAM_ratio']
            ci_low = coef - 1.96 * se
            ci_high = coef + 1.96 * se

            # Also get Total_Epi coefficient for reference
            coef_epi = fit.params['Total_Epi']
            pval_epi = fit.pvalues['Total_Epi']

            results.append({
                'cell_type': immune_cell,
                'CEACAM_coef': coef,
                'CEACAM_se': se,
                'CEACAM_pval': pval,
                'CEACAM_ci_low': ci_low,
                'CEACAM_ci_high': ci_high,
                'TotalEpi_coef': coef_epi,
                'TotalEpi_pval': pval_epi,
                'n_spots': len(df_epi),
                'n_samples': df_epi['sample'].nunique()
            })

            print(f"  CEACAM_ratio: coef={coef:.4f}, p={pval:.4f}")
            print(f"  Total_Epi:    coef={coef_epi:.4f}, p={pval_epi:.4f}")

        except Exception as e:
            print(f"  Error fitting model: {e}")
            continue

    # Create results dataframe
    df_results = pd.DataFrame(results)
    df_results = df_results.sort_values('CEACAM_coef', ascending=False)

    # ========== VISUALIZATION 1: Forest Plot ==========
    print("\n" + "="*70)
    print("Creating forest plot...")
    print("="*70)

    fig, ax = plt.subplots(figsize=(10, 8))

    y_pos = np.arange(len(df_results))
    cell_types = df_results['cell_type'].values
    coefs = df_results['CEACAM_coef'].values
    ci_lows = df_results['CEACAM_ci_low'].values
    ci_highs = df_results['CEACAM_ci_high'].values
    pvals = df_results['CEACAM_pval'].values

    # Color by direction
    colors = ['#d62728' if c > 0 else '#1f77b4' for c in coefs]

    # Plot points with error bars
    ax.errorbar(coefs, y_pos, xerr=[coefs - ci_lows, ci_highs - coefs],
                fmt='o', color='black', capsize=4, capthick=1.5, markersize=8,
                elinewidth=1.5, zorder=3)

    # Color the points
    for i, (coef, color, pval) in enumerate(zip(coefs, colors, pvals)):
        alpha = 1.0 if pval < 0.05 else 0.4
        ax.scatter(coef, i, color=color, s=100, zorder=4, alpha=alpha)

    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1.5, zorder=1)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(cell_types)
    ax.set_xlabel('CEACAM_ratio Coefficient\n(Effect of CEACAM-high vs CEACAM-low on Immune Infiltration,\nControlling for Total Epithelial Proportion)', fontsize=11)
    ax.set_title('Immune Recruitment by CEACAM Status\n(Mixed Effects Model: Immune ~ CEACAM_ratio + Total_Epi + (1|Sample))',
                 fontsize=12, fontweight='bold')

    # Add significance markers
    for i, pval in enumerate(pvals):
        sig = ''
        if pval < 0.001:
            sig = '***'
        elif pval < 0.01:
            sig = '**'
        elif pval < 0.05:
            sig = '*'

        x_pos = ci_highs[i] + 0.002
        ax.text(x_pos, i, sig, ha='left', va='center', fontsize=12, fontweight='bold')

    # Add legend
    ax.text(0.02, 0.98, 'Positive coef = More immune in CEACAM-high spots\nNegative coef = More immune in CEACAM-low spots',
            transform=ax.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    output_path = output_dir / 'ceacam_immune_forest_plot.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

    # ========== VISUALIZATION 2: Comparison Bar Chart ==========
    print("\nCreating comparison bar chart...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # Panel A: CEACAM_ratio effect (controlled)
    ax = axes[0]
    df_sorted = df_results.sort_values('CEACAM_coef', ascending=True)
    colors = ['#d62728' if c > 0 else '#1f77b4' for c in df_sorted['CEACAM_coef']]
    bars = ax.barh(range(len(df_sorted)), df_sorted['CEACAM_coef'], color=colors, edgecolor='black', alpha=0.8)

    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted['cell_type'])
    ax.set_xlabel('CEACAM_ratio Coefficient\n(Controlling for Tumor Density)', fontsize=11)
    ax.set_title('A. CEACAM Effect\n(After Controlling for Total Epithelial)', fontsize=12, fontweight='bold')

    # Add significance
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        pval = row['CEACAM_pval']
        sig = ''
        if pval < 0.001:
            sig = '***'
        elif pval < 0.01:
            sig = '**'
        elif pval < 0.05:
            sig = '*'
        x_pos = row['CEACAM_coef']
        offset = 0.002 if x_pos >= 0 else -0.002
        ha = 'left' if x_pos >= 0 else 'right'
        ax.text(x_pos + offset, i, sig, ha=ha, va='center', fontsize=11, fontweight='bold')

    # Panel B: Total_Epi effect (confounder)
    ax = axes[1]
    df_sorted2 = df_results.sort_values('TotalEpi_coef', ascending=True)
    colors2 = ['#2ca02c' if c > 0 else '#9467bd' for c in df_sorted2['TotalEpi_coef']]
    ax.barh(range(len(df_sorted2)), df_sorted2['TotalEpi_coef'], color=colors2, edgecolor='black', alpha=0.8)

    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(range(len(df_sorted2)))
    ax.set_yticklabels(df_sorted2['cell_type'])
    ax.set_xlabel('Total_Epi Coefficient\n(Tumor Density Effect)', fontsize=11)
    ax.set_title('B. Tumor Density Effect\n(Confounding Variable)', fontsize=12, fontweight='bold')

    # Add significance
    for i, (_, row) in enumerate(df_sorted2.iterrows()):
        pval = row['TotalEpi_pval']
        sig = ''
        if pval < 0.001:
            sig = '***'
        elif pval < 0.01:
            sig = '**'
        elif pval < 0.05:
            sig = '*'
        x_pos = row['TotalEpi_coef']
        offset = 0.002 if x_pos >= 0 else -0.002
        ha = 'left' if x_pos >= 0 else 'right'
        ax.text(x_pos + offset, i, sig, ha=ha, va='center', fontsize=11, fontweight='bold')

    plt.tight_layout()

    output_path = output_dir / 'ceacam_immune_comparison.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

    # ========== VISUALIZATION 3: Scatter plot for top findings ==========
    print("\nCreating scatter plots for significant findings...")

    # Get significant results
    sig_results = df_results[df_results['CEACAM_pval'] < 0.05].sort_values('CEACAM_pval')

    if len(sig_results) > 0:
        n_plots = min(4, len(sig_results))
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()

        for idx, (_, row) in enumerate(sig_results.head(4).iterrows()):
            if idx >= 4:
                break
            ax = axes[idx]
            cell_type = row['cell_type']

            # Sample data for visualization (too many points otherwise)
            df_sample = df_epi.sample(n=min(5000, len(df_epi)), random_state=42)

            scatter = ax.scatter(df_sample['CEACAM_ratio'], df_sample[cell_type],
                               c=df_sample['Total_Epi'], cmap='viridis',
                               alpha=0.3, s=10)

            # Add trend line
            z = np.polyfit(df_sample['CEACAM_ratio'], df_sample[cell_type], 1)
            p = np.poly1d(z)
            x_line = np.linspace(0, 1, 100)
            ax.plot(x_line, p(x_line), 'r-', linewidth=2, label='Trend')

            ax.set_xlabel('CEACAM_ratio (0=CEACAM-low, 1=CEACAM-high)', fontsize=10)
            ax.set_ylabel(f'{cell_type} Proportion', fontsize=10)
            ax.set_title(f'{cell_type}\nCoef={row["CEACAM_coef"]:.4f}, p={row["CEACAM_pval"]:.2e}',
                        fontsize=11, fontweight='bold')

            plt.colorbar(scatter, ax=ax, label='Total Epi')

        # Hide unused subplots
        for idx in range(len(sig_results), 4):
            axes[idx].set_visible(False)

        plt.tight_layout()

        output_path = output_dir / 'ceacam_immune_scatter_top4.png'
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()
        print(f"Saved: {output_path}")

    # ========== SAVE RESULTS ==========
    csv_path = output_dir / 'ceacam_immune_mixed_model_results.csv'
    df_results.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # ========== SUMMARY ==========
    print("\n" + "="*70)
    print("SUMMARY - Mixed Effects Model Results")
    print("="*70)
    print(f"\nModel: Immune ~ CEACAM_ratio + Total_Epi + (1|Sample)")
    print(f"Spots analyzed: {len(df_epi):,} (with ≥5% epithelial)")
    print(f"Samples: {df_epi['sample'].nunique()}")

    print(f"\n{'Cell Type':<25} {'CEACAM Coef':>12} {'p-value':>12} {'Interpretation':>25}")
    print("-"*75)

    for _, row in df_results.iterrows():
        sig = '*' if row['CEACAM_pval'] < 0.05 else ''
        if row['CEACAM_coef'] > 0:
            interp = "Recruited by CEACAM-high"
        else:
            interp = "Recruited by CEACAM-low"
        print(f"{row['cell_type']:<25} {row['CEACAM_coef']:>12.4f} {row['CEACAM_pval']:>12.4f} {interp:>25} {sig}")

    print("\n" + "="*70)
    print("Key Finding:")
    print("="*70)

    # Summarize key findings
    positive_sig = df_results[(df_results['CEACAM_coef'] > 0) & (df_results['CEACAM_pval'] < 0.05)]
    negative_sig = df_results[(df_results['CEACAM_coef'] < 0) & (df_results['CEACAM_pval'] < 0.05)]

    if len(positive_sig) > 0:
        print(f"\nCEACAM-HIGH recruits: {', '.join(positive_sig['cell_type'].tolist())}")
    else:
        print("\nNo immune cells significantly recruited by CEACAM-high")

    if len(negative_sig) > 0:
        print(f"CEACAM-LOW recruits: {', '.join(negative_sig['cell_type'].tolist())}")
    else:
        print("No immune cells significantly recruited by CEACAM-low")

    print("\n" + "="*70)
    print("Done!")
    print("="*70)


if __name__ == "__main__":
    main()
