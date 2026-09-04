#!/usr/bin/env python3
"""
Visualization Module for NicheNet Results

Generate figures from NicheNet ligand-target prediction results.

Author: Standardized Pipeline Team
Date: October 2025
"""

import os
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class NicheNetVisualizer:
    """Generate visualizations for NicheNet results."""

    def __init__(self, config):
        self.config = config
        self.nichenet_dir = os.path.join(
            config['output']['base_dir'],
            config['output'].get('nichenet_dir', '03_Results/nichenet_output')
        )
        self.output_dir = os.path.join(
            config['output']['base_dir'],
            config['output'].get('visualization_dir', '03_Results/figures')
        )

        os.makedirs(self.output_dir, exist_ok=True)

        # Visualization parameters
        self.dpi = config.get('visualization', {}).get('dpi', 300)
        self.top_n = config.get('nichenet', {}).get('top_n_ligands', 20)

    def plot_ligand_activities(self, comparison_name):
        """
        Plot ligand activity scores.

        Parameters:
        -----------
        comparison_name : str
            Name of comparison
        """
        logging.info(f"  Creating ligand activity plot for: {comparison_name}")

        nichenet_dir = os.path.join(self.nichenet_dir, comparison_name)

        if not os.path.exists(nichenet_dir):
            logging.warning(f"  Results not found for {comparison_name}, skipping")
            return

        # Find all result files
        result_files = [f for f in os.listdir(nichenet_dir)
                       if f.endswith('_ligand_activities.csv')]

        if not result_files:
            logging.warning(f"  No ligand activity files found for {comparison_name}")
            return

        # Plot each sender-receiver pair
        for result_file in result_files:
            # Load results
            results = pd.read_csv(os.path.join(nichenet_dir, result_file))

            if len(results) == 0:
                continue

            # Sort by activity score
            results = results.sort_values('pearson', ascending=False).head(self.top_n)

            # Plot
            fig, ax = plt.subplots(figsize=(10, 8))

            ax.barh(range(len(results)), results['pearson'].values)
            ax.set_yticks(range(len(results)))
            ax.set_yticklabels(results['test_ligand'].values)
            ax.set_xlabel('Ligand Activity Score (Pearson)', fontsize=12)
            ax.set_ylabel('Ligand', fontsize=12)
            ax.set_title(f'Top {self.top_n} Active Ligands\n{result_file.replace("_ligand_activities.csv", "")}',
                        fontsize=14, fontweight='bold')
            ax.invert_yaxis()

            plt.tight_layout()

            output_file = os.path.join(self.output_dir,
                                      result_file.replace('.csv', '.png'))
            plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
            plt.close()

            logging.info(f"  Saved: {output_file}")

    def generate_all_visualizations(self):
        """Generate all visualizations for all comparisons."""
        logging.info("Generating NicheNet visualizations for all comparisons")

        figures = {}

        for comparison in self.config['comparisons']:
            comparison_name = comparison['name']
            logging.info(f"\nGenerating visualizations for: {comparison_name}")

            self.plot_ligand_activities(comparison_name)

            figures[comparison_name] = {
                'ligand_activities': self.output_dir
            }

        logging.info(f"\nGenerated visualizations for {len(figures)} comparisons")

        return figures
