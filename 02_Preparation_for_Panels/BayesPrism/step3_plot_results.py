#!/usr/bin/env python3
"""
Step 3: Plot epithelial-specific CEACAM5/6 from BayesPrism deconvolution.
Violin + boxplot, R vs NR, TIGER cohort (PRJEB25780).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent

# Nature Cancer style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 6,
    'axes.labelsize': 6,
    'xtick.labelsize': 5,
    'ytick.labelsize': 5,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 2,
    'ytick.major.size': 2,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# Standard R vs NR colors
COLOR_R = '#2166AC'
COLOR_NR = '#B2182B'


def main():
    print("=" * 60)
    print("Step 3: Plot Deconvolved Epithelial CEACAM Expression")
    print("=" * 60)

    # Load deconvolved expression
    epi_expr = pd.read_csv(OUTPUT_DIR / "bayesprism_epithelial_expression.tsv",
                           sep='\t', index_col=0)
    meta = pd.read_csv(OUTPUT_DIR / "bulk_metadata.tsv", sep='\t')
    fractions = pd.read_csv(OUTPUT_DIR / "bayesprism_fractions.tsv",
                            sep='\t', index_col=0)

    print(f"Epithelial expression: {epi_expr.shape}")
    print(f"Fractions: {fractions.shape}")

    # Merge with metadata
    meta = meta.set_index('sample_id')
    common = epi_expr.index.intersection(meta.index)
    epi_expr = epi_expr.loc[common]
    meta = meta.loc[common]

    print(f"\nSamples: {len(common)}")
    print(f"R: {(meta['response_NR'] == 'R').sum()}, NR: {(meta['response_NR'] == 'N').sum()}")

    # Print epithelial fraction summary
    if 'Epithelial cells' in fractions.columns:
        epi_frac = fractions.loc[common, 'Epithelial cells']
        print(f"\nEpithelial fraction: mean={epi_frac.mean():.3f}, median={epi_frac.median():.3f}")

    # Genes to plot
    genes = ['CEACAM5', 'CEACAM6']
    available = [g for g in genes if g in epi_expr.columns]
    if not available:
        print("ERROR: No CEACAM genes found in deconvolved expression!")
        print(f"Available CEACAM genes: {[c for c in epi_expr.columns if 'CEACAM' in c]}")
        return

    print(f"\nPlotting: {available}")

    # Panel target size (same as Figure 2 stacked panels: 27 × 24.5mm)
    MM_TO_INCH = 1 / 25.4
    panel_w = 40 * MM_TO_INCH  # slightly wider for two panels side by side
    panel_h = 35 * MM_TO_INCH

    fig, axes = plt.subplots(1, len(available),
                             figsize=(panel_w * len(available), panel_h))
    if len(available) == 1:
        axes = [axes]

    for ax, gene in zip(axes, available):
        r_vals = epi_expr.loc[meta['response_NR'] == 'R', gene].values
        nr_vals = epi_expr.loc[meta['response_NR'] == 'N', gene].values

        # Log2 transform
        r_log = np.log2(r_vals + 1)
        nr_log = np.log2(nr_vals + 1)

        data = [r_log, nr_log]
        positions = [1, 2]
        colors = [COLOR_R, COLOR_NR]

        # Violin
        vp = ax.violinplot(data, positions=positions, widths=0.6,
                          showmeans=False, showmedians=False, showextrema=False)
        for i, body in enumerate(vp['bodies']):
            body.set_facecolor(colors[i])
            body.set_edgecolor('black')
            body.set_linewidth(0.5)
            body.set_alpha(0.7)

        # Boxplot
        bp = ax.boxplot(data, positions=positions, widths=0.15, patch_artist=True,
                       showfliers=False,
                       boxprops=dict(facecolor='white', linewidth=0.5),
                       whiskerprops=dict(color='black', linewidth=0.5),
                       capprops=dict(color='black', linewidth=0.5),
                       medianprops=dict(color='black', linewidth=0.8))

        # Individual points
        np.random.seed(42)
        for i, (pos, vals, color) in enumerate(zip(positions, data, colors)):
            jitter = np.random.uniform(-0.08, 0.08, len(vals))
            ax.scatter(pos + jitter, vals, c=color, s=6, alpha=0.8,
                      edgecolors='white', linewidths=0.3, zorder=3)

        # Stats
        stat, pval = stats.mannwhitneyu(nr_log, r_log, alternative='greater')
        y_max = max(np.max(r_log), np.max(nr_log))
        y_bracket = y_max * 1.1
        ax.plot([1, 1, 2, 2],
               [y_bracket, y_bracket * 1.03, y_bracket * 1.03, y_bracket],
               'k-', linewidth=0.5)
        p_text = f'P = {pval:.3f}' if pval >= 0.001 else f'P = {pval:.1e}'
        ax.text(1.5, y_bracket * 1.05, p_text, ha='center', va='bottom', fontsize=5)

        ax.set_title(f'{gene}\n(Epithelial, deconvolved)', fontsize=5, pad=2)
        ax.set_ylabel('Expression (log2)', fontsize=5)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f'R\n(n={len(r_log)})', f'NR\n(n={len(nr_log)})'], fontsize=5)

        print(f"\n{gene} (epithelial-specific, log2):")
        print(f"  R:  mean={np.mean(r_log):.3f}, median={np.median(r_log):.3f}, n={len(r_log)}")
        print(f"  NR: mean={np.mean(nr_log):.3f}, median={np.median(nr_log):.3f}, n={len(nr_log)}")
        print(f"  Mann-Whitney (NR > R): U={stat:.0f}, P={pval:.4f}")

    plt.tight_layout(pad=0.5)

    # Save
    for ext in ['png', 'svg', 'pdf']:
        out = OUTPUT_DIR / f"ceacam_epithelial_deconvolved.{ext}"
        fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: ceacam_epithelial_deconvolved.png/svg/pdf")
    plt.close()

    # Also save raw data for inspection
    result_df = pd.DataFrame({
        'sample_id': common,
        'response': meta['response_NR'].values,
    })
    for gene in available:
        result_df[f'{gene}_epi_raw'] = epi_expr.loc[common, gene].values
        result_df[f'{gene}_epi_log2'] = np.log2(epi_expr.loc[common, gene].values + 1)
    if 'Epithelial cells' in fractions.columns:
        result_df['epithelial_fraction'] = fractions.loc[common, 'Epithelial cells'].values
    result_df.to_csv(OUTPUT_DIR / "ceacam_deconvolved_data.tsv", sep='\t', index=False)
    print("Saved: ceacam_deconvolved_data.tsv")

    print("\nStep 3 complete!")


if __name__ == '__main__':
    main()
