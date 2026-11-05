#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
TCGA-STAD Data Download Script (Python Version)
Downloads FPKM gene expression and clinical data from GDC using Python API

This is the script that fetched the TCGA-STAD inputs on 2025-10-01 and wrote
three of the tables the release's analyses read: Clinical/
TCGA_STAD_clinical_data.tsv, Gene_Expression/FPKM/file_metadata.tsv and
Gene_Expression/FPKM/TCGA_STAD_gene_info.tsv. Only the first is a download; the
other two it DERIVES - file_metadata.tsv from the GDC file query's own answer
and gene_info.tsv from the gene_id/gene_name/gene_type columns of the first
STAR file it reads. See PROVENANCE.md beside this directory: re-running it today
will NOT reproduce the 2025-10-01 result, because five of the 448 files the
query returned have since been withdrawn from GDC.

BASE_DIR was an absolute path inside the author's upstream-pipeline tree, which no longer
exists in this layout. It now reads the STAD_TCGA_DOWNLOAD_DIR environment
variable and refuses to run without it, rather than carrying a path that
resolves on nobody's machine - including ours.
"""

import os
import sys
import json
import requests
import pandas as pd
import gzip
import shutil
from pathlib import Path
from datetime import datetime
import time

# Configuration
# Where to write. No default: a wrong default here would quietly download 1.9 GB
# into whatever directory the script happened to be run from, and an absent one
# would have the script create the tree and report success having written
# nowhere anybody reads. The deposit's own layout for these files is
# 02_External/Bulk/TCGA_STAD/, so that is what to point this at.
BASE_DIR = os.environ.get("STAD_TCGA_DOWNLOAD_DIR")
if not BASE_DIR:
    raise SystemExit(
        "Set STAD_TCGA_DOWNLOAD_DIR to the directory to download into.\n"
        "In the deposit's layout that is .../02_External/Bulk/TCGA_STAD, and\n"
        "this script writes 02_Raw_Data/ and 04_Logs/ beneath it.")
GDC_API_URL = "https://api.gdc.cancer.gov"
PROJECT_ID = "TCGA-STAD"

# Set up paths
RAW_DATA_DIR = os.path.join(BASE_DIR, "02_Raw_Data")
GENE_EXP_DIR = os.path.join(RAW_DATA_DIR, "Gene_Expression", "FPKM")
CLINICAL_DIR = os.path.join(RAW_DATA_DIR, "Clinical")
LOG_DIR = os.path.join(BASE_DIR, "04_Logs")

# Create log file
log_file = os.path.join(LOG_DIR, "download_python.log")
log_fh = open(log_file, 'w')

def log(message):
    """Write message to both console and log file"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_message = f"[{timestamp}] {message}"
    print(full_message)
    log_fh.write(full_message + "\n")
    log_fh.flush()

log("=" * 60)
log("TCGA-STAD Data Download Script (Python)")
log("=" * 60)
log("")

# ==============================================================================
# 1. Download Clinical Data
# ==============================================================================

log("=" * 60)
log("Step 1: Downloading Clinical Data")
log("=" * 60)

# Query clinical data
log("Querying clinical data from GDC...")
clinical_endpoint = f"{GDC_API_URL}/cases"
clinical_params = {
    "filters": json.dumps({
        "op": "in",
        "content": {
            "field": "project.project_id",
            "value": [PROJECT_ID]
        }
    }),
    "expand": "diagnoses,demographic,exposures",
    "size": 10000,
    "format": "JSON"
}

clinical_response = requests.get(clinical_endpoint, params=clinical_params)
clinical_json = clinical_response.json()

log(f"Found {len(clinical_json['data']['hits'])} cases")

# Parse clinical data
clinical_data = []
for case in clinical_json['data']['hits']:
    case_id = case['case_id']
    submitter_id = case['submitter_id']

    # Demographics
    demographic = case.get('demographic', {})

    # Diagnoses
    diagnoses_list = case.get('diagnoses', [])
    diagnosis = diagnoses_list[0] if diagnoses_list else {}

    # Exposures
    exposures_list = case.get('exposures', [])
    exposure = exposures_list[0] if exposures_list else {}

    clinical_data.append({
        'case_id': case_id,
        'submitter_id': submitter_id,
        'age_at_diagnosis': diagnosis.get('age_at_diagnosis'),
        'vital_status': demographic.get('vital_status'),
        'gender': demographic.get('gender'),
        'race': demographic.get('race'),
        'ethnicity': demographic.get('ethnicity'),
        'tumor_stage': diagnosis.get('tumor_stage'),
        'primary_diagnosis': diagnosis.get('primary_diagnosis'),
        'days_to_death': demographic.get('days_to_death'),
        'days_to_last_follow_up': diagnosis.get('days_to_last_follow_up'),
        'year_of_diagnosis': diagnosis.get('year_of_diagnosis')
    })

clinical_df = pd.DataFrame(clinical_data)
clinical_file = os.path.join(CLINICAL_DIR, "TCGA_STAD_clinical_data.tsv")
clinical_df.to_csv(clinical_file, sep='\t', index=False)
log(f"Saved clinical data to: {clinical_file}")
log(f"Number of patients: {len(clinical_df)}")
log("")

# ==============================================================================
# 2. Query and Download FPKM Gene Expression Data
# ==============================================================================

log("=" * 60)
log("Step 2: Downloading FPKM Gene Expression Data")
log("=" * 60)

# Query for FPKM files
log("Querying FPKM gene expression files from GDC...")
files_endpoint = f"{GDC_API_URL}/files"
files_params = {
    "filters": json.dumps({
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [PROJECT_ID]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_type",
                    "value": ["Gene Expression Quantification"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "analysis.workflow_type",
                    "value": ["STAR - Counts"]
                }
            }
        ]
    }),
    "fields": "file_id,file_name,cases.submitter_id,cases.samples.sample_type",
    "size": 10000,
    "format": "JSON"
}

files_response = requests.get(files_endpoint, params=files_params)
files_json = files_response.json()

log(f"Found {len(files_json['data']['hits'])} gene expression files")

# Download files
log("Downloading gene expression files...")
download_dir = os.path.join(GENE_EXP_DIR, "raw_files")
os.makedirs(download_dir, exist_ok=True)

file_metadata = []
for i, file_info in enumerate(files_json['data']['hits']):
    file_id = file_info['file_id']
    file_name = file_info['file_name']

    # Get case information
    cases = file_info.get('cases', [])
    if cases:
        submitter_id = cases[0].get('submitter_id', 'Unknown')
        samples = cases[0].get('samples', [])
        sample_type = samples[0].get('sample_type', 'Unknown') if samples else 'Unknown'
    else:
        submitter_id = 'Unknown'
        sample_type = 'Unknown'

    file_metadata.append({
        'file_id': file_id,
        'file_name': file_name,
        'submitter_id': submitter_id,
        'sample_type': sample_type
    })

    # Download file with retry logic
    download_url = f"{GDC_API_URL}/data/{file_id}"
    output_path = os.path.join(download_dir, file_name)

    if i % 10 == 0:
        log(f"Progress: {i+1}/{len(files_json['data']['hits'])} files...")

    if not os.path.exists(output_path):
        max_retries = 3
        retry_count = 0
        downloaded = False

        while retry_count < max_retries and not downloaded:
            try:
                response = requests.get(download_url, stream=True, timeout=60)
                response.raise_for_status()
                with open(output_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                downloaded = True
                time.sleep(0.5)  # Rate limiting
            except Exception as e:
                retry_count += 1
                if retry_count < max_retries:
                    log(f"Retry {retry_count}/{max_retries} for {file_name}")
                    time.sleep(2 * retry_count)  # Exponential backoff
                else:
                    log(f"Error downloading {file_name} after {max_retries} attempts: {str(e)}")
                continue

log(f"Downloaded {len(file_metadata)} files")

# Save file metadata
metadata_file = os.path.join(GENE_EXP_DIR, "file_metadata.tsv")
metadata_df = pd.DataFrame(file_metadata)
metadata_df.to_csv(metadata_file, sep='\t', index=False)
log(f"Saved file metadata to: {metadata_file}")

# ==============================================================================
# 3. Parse FPKM Data from Downloaded Files
# ==============================================================================

log("=" * 60)
log("Step 3: Parsing FPKM Data")
log("=" * 60)

log("Extracting FPKM values from downloaded files...")

# Initialize dataframe for FPKM values
fpkm_dict = {}
gene_info = None

for idx, row in metadata_df.iterrows():
    file_name = row['file_name']
    submitter_id = row['submitter_id']
    file_path = os.path.join(download_dir, file_name)

    if not os.path.exists(file_path):
        continue

    try:
        # Read the file (gzipped)
        with gzip.open(file_path, 'rt') as f:
            df = pd.read_csv(f, sep='\t', comment='#', header=0)

        # The file should have columns including gene_id and fpkm_unstrand
        if 'fpkm_unstrand' in df.columns and 'gene_id' in df.columns:
            if gene_info is None:
                gene_info = df[['gene_id', 'gene_name', 'gene_type']].copy() if 'gene_name' in df.columns else df[['gene_id']].copy()

            fpkm_dict[submitter_id] = df.set_index('gene_id')['fpkm_unstrand']

        if (idx + 1) % 50 == 0:
            log(f"Processed {idx + 1}/{len(metadata_df)} files...")

    except Exception as e:
        log(f"Error processing {file_name}: {str(e)}")
        continue

# Create FPKM matrix
log("Creating FPKM matrix...")
fpkm_df = pd.DataFrame(fpkm_dict)

# Save FPKM matrix
fpkm_matrix_file = os.path.join(GENE_EXP_DIR, "TCGA_STAD_FPKM_matrix.tsv")
fpkm_df.to_csv(fpkm_matrix_file, sep='\t')
log(f"Saved FPKM matrix to: {fpkm_matrix_file}")
log(f"Matrix dimensions: {fpkm_df.shape[0]} genes x {fpkm_df.shape[1]} samples")

# Save gene information
if gene_info is not None:
    gene_info_file = os.path.join(GENE_EXP_DIR, "TCGA_STAD_gene_info.tsv")
    gene_info.to_csv(gene_info_file, sep='\t', index=False)
    log(f"Saved gene information to: {gene_info_file}")

# ==============================================================================
# 4. Create Summary Report
# ==============================================================================

log("")
log("=" * 60)
log("Step 4: Creating Summary Report")
log("=" * 60)

summary_file = os.path.join(LOG_DIR, "download_summary.txt")
with open(summary_file, 'w') as f:
    f.write("=" * 60 + "\n")
    f.write("TCGA-STAD Data Download Summary\n")
    f.write("=" * 60 + "\n")
    f.write(f"Download Date: {datetime.now()}\n\n")
    f.write("Gene Expression Data:\n")
    f.write(f"  - Number of genes: {fpkm_df.shape[0]}\n")
    f.write(f"  - Number of samples: {fpkm_df.shape[1]}\n")
    f.write(f"  - Data file: {fpkm_matrix_file}\n\n")
    f.write("Clinical Data:\n")
    f.write(f"  - Number of patients: {len(clinical_df)}\n")
    f.write(f"  - Number of variables: {len(clinical_df.columns)}\n")
    f.write(f"  - Data file: {clinical_file}\n\n")
    f.write("All files saved to:\n")
    f.write(f"  {BASE_DIR}\n")
    f.write("=" * 60 + "\n")

log(f"Summary report saved to: {summary_file}")

# ==============================================================================
# Completion
# ==============================================================================

log("")
log("=" * 60)
log("Download Complete!")
log("=" * 60)
log(f"All data downloaded and saved successfully.")
log(f"Check logs at: {log_file}")

log_fh.close()
