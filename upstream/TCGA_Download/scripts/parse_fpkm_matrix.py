#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Parse FPKM Matrix from Downloaded TCGA Files
Creates a unified gene expression matrix from individual sample files

This writes TCGA_STAD_FPKM_matrix.tsv, which is the file the revision's TCGA
analyses read. It is 410 columns wide against 443 patients and 448 aliquot
files, and the reason is at line ~68 below: the matrix is assembled as
`fpkm_dict[submitter_id] = ...`, keyed by PATIENT. A patient with more than one
aliquot therefore keeps whichever file the loop reads last, and the assignment
leaves no record of which that was.

That is not a defect to repair here - the published analyses were run on this
matrix - but it is the reason
04_Revision_Analyses/05_R1.6_Spatial_Confounders/scripts/tcga_aliquot_provenance.py
exists: it recovers each column's tissue by comparing the column against the
patient's own STAR files rather than by reading a label, because the label was
overwritten here. Do not "fix" this keying without reading that analysis first.

BASE_DIR was an absolute path inside the author's upstream-pipeline tree; it now reads
STAD_TCGA_DOWNLOAD_DIR, like the two scripts beside it. See PROVENANCE.md.
"""

import os
import pandas as pd
from datetime import datetime

# Configuration
BASE_DIR = os.environ.get("STAD_TCGA_DOWNLOAD_DIR")
if not BASE_DIR:
    raise SystemExit(
        "Set STAD_TCGA_DOWNLOAD_DIR to the directory the download wrote into.\n"
        "In the deposit's layout that is .../02_External/Bulk/TCGA_STAD.")
RAW_FILES_DIR = os.path.join(BASE_DIR, "02_Raw_Data/Gene_Expression/FPKM/raw_files")
OUTPUT_DIR = os.path.join(BASE_DIR, "02_Raw_Data/Gene_Expression/FPKM")
METADATA_FILE = os.path.join(OUTPUT_DIR, "file_metadata.tsv")
LOG_FILE = os.path.join(BASE_DIR, "04_Logs/fpkm_parsing.log")

def log(message):
    """Write message to both console and log file"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_message = f"[{timestamp}] {message}"
    print(full_message)
    with open(LOG_FILE, 'a') as f:
        f.write(full_message + "\n")

log("=" * 60)
log("FPKM Matrix Parsing")
log("=" * 60)
log("")

# Load metadata to map file IDs to patient IDs
log("Loading file metadata...")
metadata_df = pd.read_csv(METADATA_FILE, sep='\t')
log(f"Found metadata for {len(metadata_df)} files")

# Initialize containers
fpkm_dict = {}
gene_info = None
files_processed = 0
files_skipped = 0

log("")
log("Parsing FPKM data from individual files...")

for idx, row in metadata_df.iterrows():
    file_name = row['file_name']
    submitter_id = row['submitter_id']
    file_path = os.path.join(RAW_FILES_DIR, file_name)

    if not os.path.exists(file_path):
        files_skipped += 1
        continue

    try:
        # Read the file (NOT gzipped - plain TSV)
        df = pd.read_csv(file_path, sep='\t', comment='#', header=0)

        # Column name is 'fpkm_unstranded' (not 'fpkm_unstrand')
        if 'fpkm_unstranded' in df.columns and 'gene_id' in df.columns:
            # Store gene information from first file
            if gene_info is None:
                if 'gene_name' in df.columns and 'gene_type' in df.columns:
                    gene_info = df[['gene_id', 'gene_name', 'gene_type']].copy()
                else:
                    gene_info = df[['gene_id']].copy()

            # Store FPKM values for this sample
            fpkm_dict[submitter_id] = df.set_index('gene_id')['fpkm_unstranded']
            files_processed += 1

        if (idx + 1) % 50 == 0:
            log(f"  Processed {idx + 1}/{len(metadata_df)} files ({files_processed} successful)...")

    except Exception as e:
        files_skipped += 1
        log(f"  Error processing {file_name}: {str(e)}")
        continue

log("")
log(f"Files processed successfully: {files_processed}")
log(f"Files skipped: {files_skipped}")

# Create FPKM matrix
log("")
log("Creating FPKM matrix...")
fpkm_df = pd.DataFrame(fpkm_dict)

# Save FPKM matrix
fpkm_matrix_file = os.path.join(OUTPUT_DIR, "TCGA_STAD_FPKM_matrix.tsv")
fpkm_df.to_csv(fpkm_matrix_file, sep='\t')
log(f"Saved FPKM matrix to: {fpkm_matrix_file}")
log(f"Matrix dimensions: {fpkm_df.shape[0]} genes × {fpkm_df.shape[1]} samples")

# Save gene information
if gene_info is not None:
    gene_info_file = os.path.join(OUTPUT_DIR, "TCGA_STAD_gene_info.tsv")
    gene_info.to_csv(gene_info_file, sep='\t', index=False)
    log(f"Saved gene information to: {gene_info_file}")
    log(f"Gene info contains {len(gene_info)} genes")

log("")
log("=" * 60)
log("FPKM Matrix Parsing Complete!")
log("=" * 60)
