#!/usr/bin/env python3
"""
Core NicheNet Analysis Pipeline
YAML-Driven Ligand-Target Prediction Analysis

This script orchestrates the complete NicheNet analysis workflow:
1. Load and filter single-cell data based on YAML config
2. Prepare data for NicheNet analysis (expression matrices, cell types)
3. Run NicheNet R analysis for ligand-target predictions
4. Parse and visualize results

Usage:
    python core_nichenet_pipeline.py --config path/to/config.yaml

Author: Standardized Pipeline Team
Date: October 2025
"""

import os
import sys
import yaml
import argparse
import logging
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Add prototype directory to path
sys.path.insert(0, os.path.dirname(__file__))

# Import pipeline modules
from data_preparation import NicheNetDataPreparator
from nichenet_runner import NicheNetRunner
from visualization import NicheNetVisualizer
from utils import setup_logging, validate_config, create_output_dirs


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Run NicheNet analysis pipeline from YAML configuration'
    )
    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to YAML configuration file'
    )
    parser.add_argument(
        '--steps',
        type=str,
        default='all',
        help='Steps to run: all, prepare, analyze, visualize (comma-separated)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    return parser.parse_args()


def load_config(config_path):
    """Load and validate configuration file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Validate configuration
    validate_config(config, pipeline_type='nichenet')

    return config


def run_pipeline(config, steps='all'):
    """
    Run the complete NicheNet analysis pipeline.

    Parameters:
    -----------
    config : dict
        Configuration dictionary from YAML file
    steps : str
        Which steps to run (default: 'all')
    """

    logging.info("="*80)
    logging.info("NicheNet Analysis Pipeline Started")
    logging.info("="*80)

    # Create output directories
    output_dirs = create_output_dirs(config['output'], pipeline_type='nichenet')

    # Parse steps
    if steps == 'all':
        run_steps = ['prepare', 'analyze', 'visualize']
    else:
        run_steps = [s.strip() for s in steps.split(',')]

    logging.info(f"Pipeline steps to run: {', '.join(run_steps)}")

    # Step 1: Data Preparation
    if 'prepare' in run_steps:
        logging.info("\n" + "="*80)
        logging.info("STEP 1: Data Preparation for NicheNet")
        logging.info("="*80)

        preparator = NicheNetDataPreparator(config)
        prepared_data = preparator.prepare_all_data()

        logging.info("Data preparation completed successfully")
        logging.info(f"Prepared datasets saved to: {output_dirs['prepared_data']}")

    # Step 2: NicheNet Analysis (R)
    if 'analyze' in run_steps:
        logging.info("\n" + "="*80)
        logging.info("STEP 2: NicheNet Ligand-Target Prediction")
        logging.info("="*80)

        runner = NicheNetRunner(config)
        results = runner.run_all_analyses()

        logging.info("NicheNet analysis completed successfully")
        logging.info(f"Results saved to: {output_dirs['nichenet']}")

    # Step 3: Visualization
    if 'visualize' in run_steps:
        logging.info("\n" + "="*80)
        logging.info("STEP 3: Visualization Generation")
        logging.info("="*80)

        visualizer = NicheNetVisualizer(config)
        figures = visualizer.generate_all_visualizations()

        logging.info("Visualization generation completed successfully")
        logging.info(f"Figures saved to: {output_dirs['visualization']}")

    # Pipeline complete
    logging.info("\n" + "="*80)
    logging.info("PIPELINE COMPLETED SUCCESSFULLY")
    logging.info("="*80)
    logging.info(f"\nAll results saved to: {config['output']['base_dir']}")


def main():
    """Main execution function."""
    # Parse arguments
    args = parse_arguments()

    # Load configuration
    config = load_config(args.config)

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    log_file = setup_logging(
        config['output'].get('log_dir', 'logs'),
        'nichenet_pipeline',
        log_level
    )

    logging.info(f"Configuration loaded from: {args.config}")
    logging.info(f"Log file: {log_file}")

    try:
        # Run pipeline
        run_pipeline(config, args.steps)

        return 0

    except Exception as e:
        logging.error(f"Pipeline failed with error: {str(e)}", exc_info=True)
        return 1


if __name__ == '__main__':
    sys.exit(main())
