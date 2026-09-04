#!/usr/bin/env python3
"""
Utility Functions for NicheNet Pipeline

Helper functions for logging, configuration validation, and directory management.

Author: Standardized Pipeline Team
Date: October 2025
"""

import os
import sys
import logging
from pathlib import Path
from datetime import datetime


def setup_logging(log_dir, log_name, log_level=logging.INFO):
    """
    Setup logging configuration.

    Parameters:
    -----------
    log_dir : str
        Directory for log files
    log_name : str
        Base name for log file
    log_level : int
        Logging level (default: logging.INFO)

    Returns:
    --------
    str : Path to log file
    """
    os.makedirs(log_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"{log_name}_{timestamp}.log")

    # Remove existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

    return log_file


def validate_config(config, pipeline_type='cellphonedb'):
    """
    Validate configuration structure and required fields.

    Parameters:
    -----------
    config : dict
        Configuration dictionary
    pipeline_type : str
        Type of pipeline ('cellphonedb' or 'nichenet')

    Raises:
    -------
    ValueError : If required fields are missing
    """
    errors = []

    # Check required sections
    if pipeline_type == 'nichenet':
        required_sections = ['data', 'cell_types', 'comparisons', 'nichenet', 'output']
    else:
        required_sections = ['data', 'comparisons', 'cellphonedb', 'output']

    for section in required_sections:
        if section not in config:
            errors.append(f"Missing required section: {section}")

    # Validate data section
    if 'data' in config:
        if 'h5ad_file' not in config['data']:
            errors.append("Missing required field in data section: h5ad_file")

    # Validate output section
    if 'output' in config:
        if 'base_dir' not in config['output']:
            errors.append("Missing required field in output section: base_dir")

    # NicheNet-specific validation
    if pipeline_type == 'nichenet':
        if 'cell_types' in config:
            if 'sender' not in config['cell_types']:
                errors.append("Missing required field in cell_types: sender")
            if 'receivers' not in config['cell_types']:
                errors.append("Missing required field in cell_types: receivers")

    if errors:
        raise ValueError("Configuration validation failed:\n" + "\n".join(f"  - {e}" for e in errors))

    logging.info("Configuration validated successfully")


def create_output_dirs(output_config, pipeline_type='cellphonedb'):
    """
    Create all output directories from configuration.

    Parameters:
    -----------
    output_config : dict
        Output configuration section
    pipeline_type : str
        Type of pipeline

    Returns:
    --------
    dict : Dictionary of created directory paths
    """
    base_dir = output_config['base_dir']
    os.makedirs(base_dir, exist_ok=True)

    # Standard directories
    if pipeline_type == 'nichenet':
        dirs = {
            'prepared_data': os.path.join(base_dir, output_config.get('prepared_data_dir', '03_Results/prepared_data')),
            'nichenet': os.path.join(base_dir, output_config.get('nichenet_dir', '03_Results/nichenet_output')),
            'visualization': os.path.join(base_dir, output_config.get('visualization_dir', '03_Results/figures')),
            'log': os.path.join(base_dir, output_config.get('log_dir', '04_Logs'))
        }
    else:
        dirs = {
            'prepared_data': os.path.join(base_dir, output_config.get('prepared_data_dir', '03_Results/prepared_data')),
            'cellphonedb': os.path.join(base_dir, output_config.get('cellphonedb_dir', '03_Results/cellphonedb_output')),
            'network': os.path.join(base_dir, output_config.get('network_dir', '03_Results/networks')),
            'visualization': os.path.join(base_dir, output_config.get('visualization_dir', '03_Results/figures')),
            'log': os.path.join(base_dir, output_config.get('log_dir', '04_Logs'))
        }

    # Create all directories
    for dir_name, dir_path in dirs.items():
        os.makedirs(dir_path, exist_ok=True)

    logging.info("Output directories created:")
    for dir_name, dir_path in dirs.items():
        logging.info(f"  {dir_name}: {dir_path}")

    return dirs


def format_cell_type_name(cell_type):
    """
    Format cell type name for file naming (remove special characters).

    Parameters:
    -----------
    cell_type : str
        Cell type name

    Returns:
    --------
    str : Formatted cell type name
    """
    formatted = cell_type.replace(' ', '_')
    formatted = formatted.replace('/', '_')
    formatted = formatted.replace('+', 'pos')
    formatted = formatted.replace('-', 'neg')

    return formatted
