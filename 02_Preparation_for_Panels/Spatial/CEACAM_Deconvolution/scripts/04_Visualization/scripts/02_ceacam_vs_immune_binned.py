#!/usr/bin/env python3
"""
CEACAM Score vs Immune Infiltration - Binned Analysis

Compare CEACAM-high/low score (delta) against immune cell infiltration
using binned boxplots (sample-level paired analysis).

CEACAM_delta = Epi CEACAM-high - Epi CEACAM-low
- Positive: CEACAM-high dominant spots
- Negative: CEACAM-low dominant spots

Author: Generated for Round 4 Analysis
Date: 2026-01-04
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import friedmanchisquare, wilcoxon, spearmanr
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
    print("CEACAM Score vs Immune Infiltration - Binned Analysis")
    print("="*70)

    # Paths
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'ceacam_immune_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)

    data_dir = script_dir.parent.parent / '03_Extract_Results' / 'results' / 'deconvolution_matrices'

    # Load data
    df_all = load_proportion_data(data_dir)

    # Calculate CEACAM delta
    df_all['CEACAM_delta'] = df_all['Epi CEACAM-high'] - df_all['Epi CEACAM-low']

    # Define immune cell types
    immune_cells = [
        'B cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells',
        'DC cells', 'Mast cells', 'Neutrophils', 'Plasma cells',
        'Monocytes/Macrophages'
    ]

    samples = sorted(df_all['sample'].unique())
    n_bins = 5  # Use 5 bins (quintiles)

    # Calculate global quantile thresholds
    quantiles = [i/n_bins for i in range(n_bins + 1)]
    thresholds = [df_all['CEACAM_delta'].quantile(q) for q in quantiles]

    print(f"\nCEACAM_delta quintile thresholds:")
    for i, (q, t) in enumerate(zip(quantiles, thresholds)):
        print(f"  Q{i}: {t:.4f}")

    # ========== BINNED ANALYSIS FOR EACH IMMUNE CELL TYPE ==========
    print("\n" + "="*70)
    print(f"Generating binned boxplots ({n_bins} bins) for each immune cell type...")
    print("="*70)

    results_summary = []

    for immune_cell in immune_cells:
        print(f"\nProcessing: {immune_cell}")

        # Calculate sample-level means per bin
        sample_means = {i: [] for i in range(n_bins)}
        valid_samples = []

        for sample in samples:
            df_sample = df_all[df_all['sample'] == sample]
            bin_means = []
            valid = True

            for i in range(n_bins):
                if i == n_bins - 1:
                    mask = df_sample['CEACAM_delta'] >= thresholds[i]
                else:
                    mask = (df_sample['CEACAM_delta'] >= thresholds[i]) & \
                           (df_sample['CEACAM_delta'] < thresholds[i+1])
                bin_data = df_sample[mask]

                if len(bin_data) < 5:  # Require at least 5 spots per bin
                    valid = False
                    break
                bin_means.append(bin_data[immune_cell].mean())

            if valid:
                for i, m in enumerate(bin_means):
                    sample_means[i].append(m)
                valid_samples.append(sample)

        sample_means = {k: np.array(v) for k, v in sample_means.items()}

        if len(valid_samples) < 3:
            print(f"  Skipping - too few valid samples ({len(valid_samples)})")
            continue

        # Statistical test (Friedman for paired samples across bins)
        stat, p_val = friedmanchisquare(*[sample_means[i] for i in range(n_bins)])

        # Spearman correlation (bin index vs mean proportion)
        all_bin_indices = []
        all_proportions = []
        for i in range(n_bins):
            all_bin_indices.extend([i] * len(sample_means[i]))
            all_proportions.extend(sample_means[i])
        rho, p_corr = spearmanr(all_bin_indices, all_proportions)

        results_summary.append({
            'cell_type': immune_cell,
            'friedman_p': p_val,
            'spearman_rho': rho,
            'spearman_p': p_corr,
            'n_samples': len(valid_samples),
            'trend': 'positive' if rho > 0 else 'negative'
        })

        # ========== CREATE PLOT ==========
        fig, ax = plt.subplots(figsize=(10, 6))

        # Boxplot
        data = [sample_means[i] for i in range(n_bins)]
        labels = [f'Q{i+1}' for i in range(n_bins)]

        # Color gradient from blue (CEACAM-low) to red (CEACAM-high)
        import matplotlib.cm as cm
        cmap = cm.get_cmap('coolwarm')
        colors = [cmap(i / (n_bins - 1)) for i in range(n_bins)]

        bp = ax.boxplot(data, patch_artist=True, labels=labels, widths=0.5)

        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # Add paired lines for each sample
        x_pos = list(range(1, n_bins + 1))
        for i in range(len(valid_samples)):
            ax.plot(x_pos, [sample_means[j][i] for j in range(n_bins)],
                    'o-', color='gray', alpha=0.4, markersize=4, linewidth=1)

        ax.set_ylabel(f'{immune_cell} Proportion', fontsize=12)
        ax.set_xlabel('CEACAM Score Bins\n(Q1=CEACAM-low dominant, Q5=CEACAM-high dominant)', fontsize=11)

        # Title with stats
        sig_marker = ''
        if p_val < 0.001:
            sig_marker = '***'
        elif p_val < 0.01:
            sig_marker = '**'
        elif p_val < 0.05:
            sig_marker = '*'

        ax.set_title(f'{immune_cell} vs CEACAM Score\n'
                     f'Friedman p={p_val:.3f}{sig_marker}, Spearman ρ={rho:.3f}',
                     fontsize=12, fontweight='bold')

        # Add significance bracket
        y_max = max([sample_means[i].max() for i in range(n_bins)])
        y_min = min([sample_means[i].min() for i in range(n_bins)])
        y_range = y_max - y_min

        bracket_y = y_max + 0.05 * y_range
        ax.plot([1, 1, n_bins, n_bins],
                [bracket_y, bracket_y + 0.02*y_range, bracket_y + 0.02*y_range, bracket_y],
                color='black', linewidth=1.5)

        trend_text = '↑' if rho > 0 else '↓'
        ax.text((1 + n_bins)/2, bracket_y + 0.04*y_range,
                f'p={p_val:.3f} {trend_text}', ha='center', fontsize=11, fontweight='bold')

        plt.tight_layout()

        # Safe filename
        safe_name = immune_cell.replace('/', '_').replace(' ', '_').replace('+', 'pos')
        output_path = output_dir / f'ceacam_vs_{safe_name}_5bins.png'
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_path}")
        print(f"  Friedman p={p_val:.4f}, Spearman ρ={rho:.3f}")

    # ========== SUMMARY HEATMAP ==========
    print("\n" + "="*70)
    print("Generating summary heatmap...")
    print("="*70)

    df_summary = pd.DataFrame(results_summary)
    df_summary = df_summary.sort_values('spearman_rho', ascending=False)

    fig, ax = plt.subplots(figsize=(10, 8))

    # Create heatmap data
    cell_types = df_summary['cell_type'].values
    rho_values = df_summary['spearman_rho'].values
    p_values = df_summary['friedman_p'].values

    # Plot horizontal bars
    colors = ['#d62728' if r > 0 else '#1f77b4' for r in rho_values]
    y_pos = np.arange(len(cell_types))
    bars = ax.barh(y_pos, rho_values, color=colors, edgecolor='black', alpha=0.8)

    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(cell_types)
    ax.set_xlabel('Spearman ρ (CEACAM Score vs Immune Proportion)', fontsize=12)
    ax.set_title('CEACAM Score Correlation with Immune Infiltration\n'
                 '(Red = positive, Blue = negative)', fontsize=12, fontweight='bold')

    # Add significance markers
    for i, (bar, p) in enumerate(zip(bars, p_values)):
        sig = ''
        if p < 0.001:
            sig = '***'
        elif p < 0.01:
            sig = '**'
        elif p < 0.05:
            sig = '*'

        x_pos = bar.get_width()
        offset = 0.02 if x_pos >= 0 else -0.02
        ha = 'left' if x_pos >= 0 else 'right'
        ax.text(x_pos + offset, bar.get_y() + bar.get_height()/2,
                sig, ha=ha, va='center', fontsize=12, fontweight='bold')

    ax.set_xlim(-0.5, 0.5)
    plt.tight_layout()

    output_path = output_dir / 'ceacam_immune_correlation_summary.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

    # Save summary CSV
    csv_path = output_dir / 'ceacam_immune_correlation_summary.csv'
    df_summary.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\n{'Cell Type':<25} {'Spearman ρ':>12} {'Friedman p':>12} {'Trend':>10}")
    print("-"*65)
    for _, row in df_summary.iterrows():
        sig = '*' if row['friedman_p'] < 0.05 else ''
        print(f"{row['cell_type']:<25} {row['spearman_rho']:>12.3f} {row['friedman_p']:>12.4f} {row['trend']:>10} {sig}")

    print("\n" + "="*70)
    print("Done!")
    print("="*70)


if __name__ == "__main__":
    main()
