#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
NicheNet Runner Module

Executes NicheNet R analysis using prepared data.

Author: Standardized Pipeline Team
Date: October 2025
"""

import os
import logging
import subprocess
from pathlib import Path


class NicheNetRunner:
    """Run NicheNet analysis via R script."""

    def __init__(self, config):
        self.config = config
        self.r_script = os.path.join(
            os.path.dirname(__file__),
            'run_nichenet_analysis.R'
        )
        self.output_dir = os.path.join(
            config['output']['base_dir'],
            config['output'].get('nichenet_dir', '03_Results/nichenet_output')
        )
        self.prepared_data_dir = os.path.join(
            config['output']['base_dir'],
            config['output'].get('prepared_data_dir', '03_Results/prepared_data')
        )

        os.makedirs(self.output_dir, exist_ok=True)

    def run_nichenet_for_comparison(self, comparison_name):
        """
        Run NicheNet analysis for a single comparison.

        Parameters:
        -----------
        comparison_name : str
            Name of comparison

        Returns:
        --------
        dict : Analysis results information
        """
        logging.info(f"  Running NicheNet analysis for: {comparison_name}")

        # Data directory
        data_dir = os.path.join(self.prepared_data_dir, comparison_name)

        # Output directory
        output_dir = os.path.join(self.output_dir, comparison_name)
        os.makedirs(output_dir, exist_ok=True)

        # R executable
        r_executable = self.config.get('computational', {}).get('r_executable', 'Rscript')

        # Build R command
        cmd = [
            r_executable,
            self.r_script,
            data_dir,
            output_dir,
            self.config['cell_types']['sender'],
            ','.join(self.config['cell_types']['receivers'])
        ]

        logging.info(f"  R command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            logging.info(f"  NicheNet analysis completed for {comparison_name}")
            logging.info(f"  Results saved to: {output_dir}")

            return {
                'status': 'success',
                'output_dir': output_dir
            }

        except subprocess.CalledProcessError as e:
            logging.error(f"  NicheNet analysis failed for {comparison_name}")
            logging.error(f"  Error: {e.stderr}")

            return {
                'status': 'failed',
                'error': str(e)
            }

        except FileNotFoundError:
            logging.error(f"  R executable not found: {r_executable}")
            logging.info("  Skipping NicheNet analysis")

            return {
                'status': 'skipped',
                'reason': 'R not available'
            }

    def run_all_analyses(self):
        """Run NicheNet analysis for all comparisons."""
        logging.info("Starting NicheNet analysis for all comparisons")

        results = {}

        for comparison in self.config['comparisons']:
            comparison_name = comparison['name']
            logging.info(f"\nProcessing comparison: {comparison_name}")

            result = self.run_nichenet_for_comparison(comparison_name)
            results[comparison_name] = result

        return results
