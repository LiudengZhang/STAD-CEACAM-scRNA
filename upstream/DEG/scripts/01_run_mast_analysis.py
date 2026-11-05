#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Treg MAST Analysis Script
Uses the unified deg_pathway_analysis_library for MAST analysis of Treg cells
Handles single cell type mode (Treg only, not iterating through multiple cell types)
"""

import sys
import argparse
import logging
from pathlib import Path
from datetime import datetime

library_path = Path(__file__).resolve().parents[4] / '00_Script_Library'
sys.path.insert(0, str(library_path))

from deg_pathway_analysis_library import (
    DEGConfig,
    DataLoader,
    MASTStrategy,
)

import pandas as pd
import numpy as np


class TregMASTAnalyzer:

    def __init__(self, config_path: str):
        self.config = DEGConfig(config_path)
        self.setup_logging()

        self.data_loader = DataLoader(self.config, self.logger)

        self.strategies = {
            'mast': MASTStrategy(),
        }

    def setup_logging(self):
        output_dir = Path(self.config.get('data', 'output_dir'))
        log_dir = output_dir.parent / '04_Logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_file = log_dir / f"treg_mast_analysis_{timestamp}.log"

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)

    def run_treg_mast_analysis(self):

        methods = self.config.get('analysis', 'methods', default=['wilcoxon'])
        if not methods:
            methods = ['wilcoxon']

        primary_method = methods[0]

        if primary_method not in self.strategies:
            self.logger.error(f"Unknown method: {primary_method}")
            return False

        strategy = self.strategies[primary_method]

        # Check method dependencies - NO FALLBACKS, fail cleanly
        if not strategy.check_dependencies():
            self.logger.error(f"{primary_method} dependencies not available!")
            self.logger.error(f"Description: {strategy.get_description()}")
            self.logger.error(f"Please ensure the correct conda environment is activated")
            self.logger.error(f"For MAST: conda activate r_bayesprism")
            return False

        adata = self.data_loader.load_h5ad_data()
        self.data_loader.validate_data()
        self.data_loader.prepare_for_deg()

        cell_type_name = self.config.get('data', 'cell_type_name', default='Treg')
        comparisons_to_run = self.config.get('analysis', 'comparisons', default=[])
        min_cells = self.config.get('analysis', 'min_cells_per_group', default=10)

        all_comparison_configs = self.config.get_comparison_configs()

        output_base = Path(self.config.get('data', 'output_dir')) / '01_DEG_Analysis'
        output_base.mkdir(parents=True, exist_ok=True)

        method_params = self.config.get('analysis', f'{primary_method}_parameters',
                                        default=self.config.get('analysis', 'mast_parameters', default={}))
        max_cells = method_params.get('max_cells_mast', 8000) if 'mast' in primary_method else None

        if max_cells:
            self.logger.info(f"  Max cells: {max_cells}")
        if 'model_formula' in method_params:
            self.logger.info(f"  Model: {method_params.get('model_formula', 'default')}")

        all_results = {}

        for comp_name in comparisons_to_run:
            if comp_name not in all_comparison_configs:
                continue

            comp_config = all_comparison_configs[comp_name]

            adata_comp = self.data_loader.filter_by_comparison(comp_config)

            if adata_comp.n_obs < min_cells * 2:
                continue

            group_col = comp_config['column']
            group_counts = adata_comp.obs[group_col].value_counts()
            for group, count in group_counts.items():
                self.logger.info(f"  {group}: {count} cells")

            try:

                strategy_kwargs = {}
                if primary_method == 'mast' and max_cells:
                    strategy_kwargs['max_cells_mast'] = max_cells
                    if 'min_pct_expressed' in method_params:
                        strategy_kwargs['min_pct_expressed'] = method_params['min_pct_expressed']

                deg_results = strategy.run(adata_comp, comp_config, **strategy_kwargs)

                if deg_results is not None and len(deg_results) > 0:

                    sig_results = self._filter_significant_genes(deg_results)

                    if len(sig_results) > 0:
                        n_up = len(sig_results[sig_results['logfoldchanges'] > 0])
                        n_down = len(sig_results[sig_results['logfoldchanges'] < 0])

                    self._save_results(
                        deg_results,
                        sig_results,
                        comp_name,
                        cell_type_name,
                        output_base,
                        group_counts
                    )

                    all_results[comp_name] = {
                        'all_results': deg_results,
                        'significant': sig_results,
                        'n_cells': adata_comp.n_obs,
                        'group_counts': group_counts.to_dict()
                    }
                else:
                    self.logger.warning(f"No results from {primary_method} for {comp_name}")

            except Exception as e:
                self.logger.error(f"Error in {primary_method} analysis for {comp_name}: {str(e)}")
                import traceback
                self.logger.error(traceback.format_exc())
                continue

        for comp_name, results in all_results.items():
            self.logger.info(f"  {comp_name}:")
            self.logger.info(f"    Total genes: {len(results['all_results'])}")
            self.logger.info(f"    Significant: {len(results['significant'])}")
            self.logger.info(f"    Cells: {results['n_cells']}")

        return len(all_results) > 0

    def _filter_significant_genes(self, deg_results: pd.DataFrame) -> pd.DataFrame:
        """Filter significant genes based on thresholds"""
        logfc_thresh = self.config.get('analysis', 'logfc_threshold', default=0.25)
        pval_thresh = self.config.get('analysis', 'pval_threshold', default=0.05)

        mask = (np.abs(deg_results['logfoldchanges']) >= logfc_thresh) & \
               (deg_results['pvals_adj'] < pval_thresh)

        return deg_results[mask].copy()

    def _save_results(self, all_results, sig_results, comp_name,
                     cell_type_name, output_base, group_counts):
        comp_dir = output_base / comp_name / f"{cell_type_name}_cell_type"
        comp_dir.mkdir(parents=True, exist_ok=True)

        methods = self.config.get('analysis', 'methods', default=['wilcoxon'])
        method_name = methods[0] if methods else 'wilcoxon'

        all_file = comp_dir / f"{cell_type_name}_{method_name}_deg.csv"
        all_results.to_csv(all_file, index=True)

        sig_file = comp_dir / f"{cell_type_name}_{method_name}_deg_significant.csv"
        sig_results.to_csv(sig_file, index=True)

        summary_data = {
            'cell_type': cell_type_name,
            'comparison': comp_name,
            'method': method_name,
            'n_genes_tested': len(all_results),
            'n_significant': len(sig_results),
            'n_upregulated': len(sig_results[sig_results['logfoldchanges'] > 0]),
            'n_downregulated': len(sig_results[sig_results['logfoldchanges'] < 0]),
            'group_counts': str(group_counts.to_dict())
        }

        summary_df = pd.DataFrame([summary_data])
        summary_file = comp_dir / f"{cell_type_name}_{method_name}_summary.csv"
        summary_df.to_csv(summary_file, index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Treg MAST Analysis using unified library",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        """
    )

    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to YAML configuration file'
    )

    args = parser.parse_args()

    if not Path(args.config).exists():
        sys.exit(1)

    analyzer = TregMASTAnalyzer(args.config)
    success = analyzer.run_treg_mast_analysis()

    if success:
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
