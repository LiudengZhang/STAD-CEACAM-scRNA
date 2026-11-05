#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
CEACAM vs Immune Infiltration - Mixed Effects Model Analysis
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
    csv_files = sorted(Path(data_dir).glob('*_proportions_14types.csv'))

    all_data = []
    for csv_path in csv_files:
        sample_name = csv_path.stem.replace('_proportions_14types', '')
        df = pd.read_csv(csv_path, index_col=0)
        df['sample'] = sample_name
        df['spot_id'] = df.index
        all_data.append(df)

    df_all = pd.concat(all_data, ignore_index=True)
    return df_all


def main():

    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / 'results' / 'ceacam_immune_mixed_model'
    output_dir.mkdir(parents=True, exist_ok=True)

    data_dir = script_dir.parent.parent / '03_Extract_Results' / 'results' / 'deconvolution_matrices'

    df_all = load_proportion_data(data_dir)

    df_all['Total_Epi'] = df_all['Epi CEACAM-high'] + df_all['Epi CEACAM-low']

    # CEACAM ratio (avoid division by zero)
    # Only calculate for spots with sufficient epithelial content
    epi_threshold = 0.05  # Minimum 5% epithelial to be considered
    df_epi = df_all[df_all['Total_Epi'] >= epi_threshold].copy()
    df_epi['CEACAM_ratio'] = df_epi['Epi CEACAM-high'] / df_epi['Total_Epi']

    immune_cells = [
        'B cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells',
        'DC cells', 'Mast cells', 'Neutrophils', 'Plasma cells',
        'Monocytes/Macrophages'
    ]

    results = []

    for immune_cell in immune_cells:

        safe_col = immune_cell.replace(' ', '_').replace('/', '_').replace('+', 'pos')
        df_epi[safe_col] = df_epi[immune_cell]

        formula = f"{safe_col} ~ CEACAM_ratio + Total_Epi"

        try:
            model = smf.mixedlm(formula, df_epi, groups=df_epi['sample'])
            fit = model.fit(method='lbfgs', maxiter=1000)

            coef = fit.params['CEACAM_ratio']
            se = fit.bse['CEACAM_ratio']
            pval = fit.pvalues['CEACAM_ratio']
            ci_low = coef - 1.96 * se
            ci_high = coef + 1.96 * se

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

        except Exception as e:
            continue

    df_results = pd.DataFrame(results)
    df_results = df_results.sort_values('CEACAM_coef', ascending=False)


    csv_path = output_dir / 'ceacam_immune_mixed_model_results.csv'
    df_results.to_csv(csv_path, index=False)

    for _, row in df_results.iterrows():
        sig = '*' if row['CEACAM_pval'] < 0.05 else ''
        if row['CEACAM_coef'] > 0:
            interp = "Recruited by CEACAM-high"
        else:
            interp = "Recruited by CEACAM-low"

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


if __name__ == "__main__":
    main()
