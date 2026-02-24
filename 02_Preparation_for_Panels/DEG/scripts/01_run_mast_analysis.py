#!/usr/bin/env python3
"""
Treg MAST Analysis Script
Uses the unified deg_pathway_analysis_library for MAST analysis of Treg cells
Handles single cell type mode (Treg only, not iterating through multiple cell types)
"""

import sys
import os
import argparse
import logging
from pathlib import Path
from datetime import datetime

# Import library components (co-located)
from deg_pathway_analysis_library import (
    DEGConfig,
    DataLoader,
    DEGAnalyzer,
    MASTStrategy,
    WilcoxonStrategy,
    DESeq2Strategy,
    TTestStrategy
)

import pandas as pd
import numpy as np


class TregMASTAnalyzer:
    """Specialized MAST analyzer for single cell type (Treg) analysis"""

    def __init__(self, config_path: str):
        """Initialize Treg MAST analyzer"""
        self.config = DEGConfig(config_path)
        self.setup_logging()

        # Initialize components
        self.data_loader = DataLoader(self.config, self.logger)

        # Initialize available strategies
        self.strategies = {
            'mast': MASTStrategy(),
            'wilcoxon': WilcoxonStrategy(),
            'deseq2': DESeq2Strategy(),
            'ttest': TTestStrategy()
        }

    def setup_logging(self):
        """Setup logging configuration"""
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
        self.logger.info(f"Logging to {log_file}")

    def run_treg_mast_analysis(self):
        """Run DEG analysis for Treg cells (method determined by config)"""
        self.logger.info("=" * 70)
        self.logger.info("TREG DIFFERENTIAL EXPRESSION ANALYSIS")
        self.logger.info("=" * 70)

        # Get method from config
        methods = self.config.get('analysis', 'methods', default=['wilcoxon'])
        if not methods:
            methods = ['wilcoxon']

        primary_method = methods[0]
        self.logger.info(f"Primary method: {primary_method}")

        # Check if method is available
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

        self.logger.info(f"✓ {primary_method} method available")
        self.logger.info(f"  {strategy.get_description()}")

        # Load data
        self.logger.info("\nLoading Treg data...")
        adata = self.data_loader.load_h5ad_data()
        self.data_loader.validate_data()
        self.data_loader.prepare_for_deg()

        # Get configuration
        cell_type_name = self.config.get('data', 'cell_type_name', default='Treg')
        comparisons_to_run = self.config.get('analysis', 'comparisons', default=[])
        min_cells = self.config.get('analysis', 'min_cells_per_group', default=10)

        self.logger.info(f"Cell type: {cell_type_name}")
        self.logger.info(f"Comparisons to run: {comparisons_to_run}")

        # Get all comparison configurations
        all_comparison_configs = self.config.get_comparison_configs()

        # Create output base directory
        output_base = Path(self.config.get('data', 'output_dir')) / '01_DEG_Analysis'
        output_base.mkdir(parents=True, exist_ok=True)

        # Method parameters
        method_params = self.config.get('analysis', f'{primary_method}_parameters',
                                        default=self.config.get('analysis', 'mast_parameters', default={}))
        max_cells = method_params.get('max_cells_mast', 8000) if 'mast' in primary_method else None

        self.logger.info(f"\n{primary_method.upper()} Parameters:")
        if max_cells:
            self.logger.info(f"  Max cells: {max_cells}")
        if 'model_formula' in method_params:
            self.logger.info(f"  Model: {method_params.get('model_formula', 'default')}")

        # Process each comparison
        all_results = {}

        for comp_name in comparisons_to_run:
            if comp_name not in all_comparison_configs:
                self.logger.warning(f"Unknown comparison: {comp_name}")
                continue

            comp_config = all_comparison_configs[comp_name]

            self.logger.info("\n" + "=" * 70)
            self.logger.info(f"COMPARISON: {comp_name}")
            self.logger.info(f"Description: {comp_config['description']}")
            self.logger.info("=" * 70)

            # Filter data for this comparison
            adata_comp = self.data_loader.filter_by_comparison(comp_config)

            if adata_comp.n_obs < min_cells * 2:
                self.logger.warning(f"Insufficient cells for {comp_name}: {adata_comp.n_obs}")
                continue

            # Check group sizes
            group_col = comp_config['column']
            group_counts = adata_comp.obs[group_col].value_counts()
            self.logger.info(f"\nGroup sizes:")
            for group, count in group_counts.items():
                self.logger.info(f"  {group}: {count} cells")

            # Run DEG analysis
            try:
                self.logger.info(f"\nRunning {primary_method} analysis...")
                self.logger.info(f"  Total cells: {adata_comp.n_obs}")

                # Build kwargs for strategy
                strategy_kwargs = {}
                if primary_method == 'mast' and max_cells:
                    strategy_kwargs['max_cells_mast'] = max_cells
                    if 'min_pct_expressed' in method_params:
                        strategy_kwargs['min_pct_expressed'] = method_params['min_pct_expressed']

                deg_results = strategy.run(adata_comp, comp_config, **strategy_kwargs)

                if deg_results is not None and len(deg_results) > 0:
                    self.logger.info(f"✓ {primary_method.upper()} analysis completed: {len(deg_results)} genes tested")

                    # Filter significant results
                    sig_results = self._filter_significant_genes(deg_results)
                    self.logger.info(f"  Significant DEGs: {len(sig_results)}")

                    if len(sig_results) > 0:
                        n_up = len(sig_results[sig_results['logfoldchanges'] > 0])
                        n_down = len(sig_results[sig_results['logfoldchanges'] < 0])
                        self.logger.info(f"    Upregulated: {n_up}")
                        self.logger.info(f"    Downregulated: {n_down}")

                    # Save results
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

        # Final summary
        self.logger.info("\n" + "=" * 70)
        self.logger.info(f"{primary_method.upper()} ANALYSIS SUMMARY")
        self.logger.info("=" * 70)
        self.logger.info(f"Comparisons processed: {len(all_results)}")
        for comp_name, results in all_results.items():
            self.logger.info(f"  {comp_name}:")
            self.logger.info(f"    Total genes: {len(results['all_results'])}")
            self.logger.info(f"    Significant: {len(results['significant'])}")
            self.logger.info(f"    Cells: {results['n_cells']}")

        self.logger.info(f"\nResults saved to: {output_base}")

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
        """Save DEG results to files"""
        # Create comparison directory
        comp_dir = output_base / comp_name / f"{cell_type_name}_cell_type"
        comp_dir.mkdir(parents=True, exist_ok=True)

        # Determine method name for file naming
        methods = self.config.get('analysis', 'methods', default=['wilcoxon'])
        method_name = methods[0] if methods else 'wilcoxon'

        # Save all results
        all_file = comp_dir / f"{cell_type_name}_{method_name}_deg.csv"
        all_results.to_csv(all_file, index=True)
        self.logger.info(f"  Saved all results: {all_file}")

        # Save significant results
        sig_file = comp_dir / f"{cell_type_name}_{method_name}_deg_significant.csv"
        sig_results.to_csv(sig_file, index=True)
        self.logger.info(f"  Saved significant results: {sig_file}")

        # Save summary
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
        self.logger.info(f"  Saved summary: {summary_file}")


def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(
        description="Treg MAST Analysis using unified library",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run pre-treatment analysis
  python 01_run_mast_analysis.py --config ../01_Config/treg_mast_pre_response.yaml

  # Run post-treatment analysis
  python 01_run_mast_analysis.py --config ../01_Config/treg_mast_post_response.yaml

  # Run complete analysis (both comparisons)
  python 01_run_mast_analysis.py --config ../01_Config/treg_mast_complete.yaml
        """
    )

    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to YAML configuration file'
    )

    args = parser.parse_args()

    # Verify config file exists
    if not Path(args.config).exists():
        print(f"Error: Configuration file not found: {args.config}")
        sys.exit(1)

    # Run analysis
    analyzer = TregMASTAnalyzer(args.config)
    success = analyzer.run_treg_mast_analysis()

    if success:
        print("\n✓ Treg MAST analysis completed successfully!")
        sys.exit(0)
    else:
        print("\n✗ Treg MAST analysis failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
