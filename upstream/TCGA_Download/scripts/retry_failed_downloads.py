#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Retry Failed TCGA-STAD Downloads
Attempts to download files that failed in the main download process

Run on 2025-10-01 against the four files FAILED_FILES names. All four failed
again with a permanent `404 NOT FOUND` from the GDC data endpoint, five attempts
each; the log of that run is what records the gap between the 448 files the
query returned and the 443 on disk. The fifth absent file
(79946f8f-c71f-477a-985b-0fb923937225) is not in this list because it failed on
a dropped connection rather than a 404 and was never added to the retry set.

Kept as the record of that attempt, not as something to run again: a 404 from
GDC means the file has been withdrawn from the release, and no number of retries
changes that. See PROVENANCE.md.

BASE_DIR was an absolute path inside the author's upstream-pipeline tree; it now reads
STAD_TCGA_DOWNLOAD_DIR, like the download script beside it.
"""

import os
import requests
import time
from datetime import datetime

# Configuration
BASE_DIR = os.environ.get("STAD_TCGA_DOWNLOAD_DIR")
if not BASE_DIR:
    raise SystemExit(
        "Set STAD_TCGA_DOWNLOAD_DIR to the directory the download wrote into.\n"
        "In the deposit's layout that is .../02_External/Bulk/TCGA_STAD.")
GDC_API_URL = "https://api.gdc.cancer.gov"
DOWNLOAD_DIR = os.path.join(BASE_DIR, "02_Raw_Data/Gene_Expression/FPKM/raw_files")
LOG_FILE = os.path.join(BASE_DIR, "04_Logs/retry_downloads.log")

# Failed file IDs
FAILED_FILES = [
    "1b06673a-c03b-48e4-a220-e9195c515c9e",
    "691e4c86-7ad9-40b7-9312-2d625472bf97",
    "1455889c-26f3-4d83-aa90-271d6dbdea65",
    "71106ad8-a2c4-4337-9682-c407356d39ef"
]

def log(message):
    """Write message to both console and log file"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_message = f"[{timestamp}] {message}"
    print(full_message)
    with open(LOG_FILE, 'a') as f:
        f.write(full_message + "\n")

log("=" * 60)
log("Retry Failed Downloads Script")
log("=" * 60)
log(f"Attempting to download {len(FAILED_FILES)} failed files")
log("")

# Enhanced retry settings
MAX_RETRIES = 5
TIMEOUT = 120
BASE_DELAY = 5

successful = 0
failed = 0

for file_id in FAILED_FILES:
    file_name = f"{file_id}.rna_seq.augmented_star_gene_counts.tsv"
    output_path = os.path.join(DOWNLOAD_DIR, file_name)

    if os.path.exists(output_path):
        log(f"File {file_name} already exists, skipping")
        successful += 1
        continue

    log(f"Attempting to download {file_id}...")
    download_url = f"{GDC_API_URL}/data/{file_id}"

    retry_count = 0
    downloaded = False

    while retry_count < MAX_RETRIES and not downloaded:
        try:
            log(f"  Attempt {retry_count + 1}/{MAX_RETRIES}...")
            response = requests.get(download_url, stream=True, timeout=TIMEOUT)
            response.raise_for_status()

            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            downloaded = True
            successful += 1
            log(f"  ✓ Successfully downloaded {file_name}")
            time.sleep(2)  # Rate limiting

        except Exception as e:
            retry_count += 1
            if retry_count < MAX_RETRIES:
                delay = BASE_DELAY * retry_count
                log(f"  ✗ Failed: {str(e)}")
                log(f"  Waiting {delay} seconds before retry...")
                time.sleep(delay)
            else:
                failed += 1
                log(f"  ✗ Failed to download {file_name} after {MAX_RETRIES} attempts")
                log(f"  Final error: {str(e)}")

log("")
log("=" * 60)
log("Retry Summary")
log("=" * 60)
log(f"Successfully downloaded: {successful}/{len(FAILED_FILES)}")
log(f"Failed: {failed}/{len(FAILED_FILES)}")
log("=" * 60)
