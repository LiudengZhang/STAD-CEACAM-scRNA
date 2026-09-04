#!/usr/bin/env python3
# The interpreter, R and library locations below were machine-specific
# paths under the author's home directory. They are replaced by a
# /path/to/home placeholder for deposition; point them at your own conda
# environment before running this script.
"""
Batch Process All Samples with GraphST
15 Cell Types (12 Non-Epithelial + 3 Epithelial) - Korean Gut Dataset

This script runs GraphST on all 10 Visium samples sequentially.

Author: Generated for Round 4 Analysis
Date: 2025-11-21
"""

import os
import sys
import yaml
import subprocess
import time
from datetime import datetime

# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)

def log(msg):
    """Print with timestamp and flush."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def load_config(config_path):
    """Load configuration."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def run_sample(config_path, sample_name, log_dir, status_file):
    """Run GraphST on a single sample with real-time output."""
    log(f"{'='*60}")
    log(f"STARTING: {sample_name}")
    log(f"{'='*60}")

    # Update status file
    with open(status_file, 'a') as sf:
        sf.write(f"{datetime.now().strftime('%H:%M:%S')} | {sample_name} | STARTED\n")
        sf.flush()

    # Prepare command - use full path to conda environment's Python
    script_path = os.path.join(os.path.dirname(__file__), '01_run_graphst_per_sample.py')
    python_path = '/path/to/home/miniforge3/envs/Liudeng_Python_310/bin/python'
    cmd = [
        python_path, '-u',  # -u for unbuffered output
        script_path,
        config_path,
        sample_name
    ]

    # Run with REAL-TIME logging
    log_file = os.path.join(log_dir, f'{sample_name}_graphst.log')
    start = time.time()

    with open(log_file, 'w') as f:
        # Use Popen for real-time streaming
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1  # Line buffered
        )

        # Stream output in real-time
        for line in process.stdout:
            print(line, end='', flush=True)  # Print to console
            f.write(line)  # Write to log
            f.flush()

        process.wait()
        elapsed = time.time() - start

    if process.returncode == 0:
        log(f"SUCCESS: {sample_name} completed in {elapsed/60:.1f} minutes")
        with open(status_file, 'a') as sf:
            sf.write(f"{datetime.now().strftime('%H:%M:%S')} | {sample_name} | COMPLETED ({elapsed/60:.1f} min)\n")
            sf.flush()
        return True
    else:
        log(f"FAILED: {sample_name} (check {log_file})")
        with open(status_file, 'a') as sf:
            sf.write(f"{datetime.now().strftime('%H:%M:%S')} | {sample_name} | FAILED\n")
            sf.flush()
        return False


def main(config_path):
    """Main execution."""
    log("="*60)
    log("GraphST Batch Processing - All Samples")
    log("="*60)

    # Load config
    config = load_config(config_path)

    # Get all samples
    samples = list(config['samples'].keys())
    log(f"Samples to process: {len(samples)}")
    for s in samples:
        print(f"  - {s}", flush=True)

    # Create log directory
    log_dir = os.path.join(os.path.dirname(config_path), '../logs')
    os.makedirs(log_dir, exist_ok=True)

    # Create status file for quick monitoring
    status_file = os.path.join(log_dir, 'status.txt')
    with open(status_file, 'w') as sf:
        sf.write(f"GraphST Batch Run - Started {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        sf.write("="*50 + "\n")
        sf.flush()
    log(f"Status file: {status_file}")

    # Process each sample
    results = {}
    overall_start = time.time()

    for i, sample_name in enumerate(samples, 1):
        log(f"")
        log(f"### Sample {i}/{len(samples)}: {sample_name} ###")

        success = run_sample(config_path, sample_name, log_dir, status_file)
        results[sample_name] = success

        # Estimated time remaining
        if i < len(samples):
            elapsed = time.time() - overall_start
            avg_time = elapsed / i
            remaining = avg_time * (len(samples) - i)
            log(f"Estimated time remaining: {remaining/60:.1f} minutes")

    # Summary
    overall_elapsed = time.time() - overall_start

    log("")
    log("="*60)
    log("BATCH PROCESSING COMPLETE")
    log("="*60)
    log(f"Total time: {overall_elapsed/60:.1f} minutes ({overall_elapsed/3600:.1f} hours)")

    successful = [s for s, success in results.items() if success]
    failed = [s for s, success in results.items() if not success]

    log(f"Successful: {len(successful)}/{len(samples)}")
    for s in successful:
        print(f"    OK {s}", flush=True)

    if failed:
        log(f"Failed: {len(failed)}/{len(samples)}")
        for s in failed:
            print(f"    XX {s}", flush=True)

    # Write final status
    with open(status_file, 'a') as sf:
        sf.write("="*50 + "\n")
        sf.write(f"COMPLETED: {len(successful)}/{len(samples)} successful\n")
        sf.write(f"Total time: {overall_elapsed/60:.1f} minutes\n")
        sf.flush()

    log(f"Logs: {log_dir}")
    log(f"Status: {status_file}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python 02_batch_process_all_samples.py <config_file>")
        print("\nExample:")
        print("  python 02_batch_process_all_samples.py config/graphst_config.yaml")
        sys.exit(1)

    config_path = sys.argv[1]

    if not os.path.exists(config_path):
        print(f"Error: Configuration file not found: {config_path}")
        sys.exit(1)

    main(config_path)
