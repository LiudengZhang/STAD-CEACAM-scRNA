> **Upstream record.** Paths in this document refer to the upstream
> processing pipeline that produced the deposited intermediates; it is not
> part of this release's run path. The files under `upstream/` are included
> as a record of how the inputs were produced and are not executed by
> `./run` or `_run_all_panels.sh`.

# TCGA-STAD — how the bulk inputs were obtained, and why re-running will not reproduce them

Added 2026-10-04. Three shipped analyses read tables that nothing in the release
explained the origin of. These are the scripts that produced them, moved here
from `submission-tree/01_Raw_Inputs/02_External/Bulk/TCGA_STAD/01_Scripts/` with one
change each: `BASE_DIR` was an absolute path inside the author's **upstream-pipeline**
tree, which does not exist in this layout, and now reads
`STAD_TCGA_DOWNLOAD_DIR` and refuses to run without it.

    scripts/download_tcga_python.py     the GDC query and fetch, 2025-10-01
    scripts/retry_failed_downloads.py   the retry of the four that failed
    scripts/parse_fpkm_matrix.py        assembles TCGA_STAD_FPKM_matrix.tsv

## What each file is, and which of them is actually a download

Only one of the four tables is downloaded as such. The other three are derived
by these scripts from what the GDC answered, which is why they belong beside the
code and not in a list of external inputs:

| file | md5 | bytes | how it is produced |
|---|---|---|---|
| `Clinical/TCGA_STAD_clinical_data.tsv` | `56a280aac796b4a69312e7a81058f018` | 60,392 | **downloaded**: the GDC `/cases` endpoint for project `TCGA-STAD`, 443 cases |
| `…/FPKM/file_metadata.tsv` | `00f7f58738dbbd024c8dd882806282c4` | 62,979 | **derived**: the GDC file query's own answer, 448 rows (`file_id`, `file_name`, `submitter_id`, `sample_type`) |
| `…/FPKM/TCGA_STAD_gene_info.tsv` | `4bc63bc8139a72c045f8e60f223074fd` | 2,457,448 | **derived**: the `gene_id`/`gene_name`/`gene_type` columns of the first STAR file read |
| `…/FPKM/TCGA_STAD_FPKM_matrix.tsv` | `9b3f83d7a795e2a16b8c92500732e461` | 148,357,653 | **derived** by `parse_fpkm_matrix.py` from the 443 STAR files |

The 443 `*.rna_seq.augmented_star_gene_counts.tsv` files are 1.8 GB and are not
deposited. `file_metadata.tsv` carries the `file_id` of every one of them, so it
is the manifest: `https://api.gdc.cancer.gov/data/<file_id>` fetches each.

## Re-running this will NOT reproduce the 2025-10-01 result

The GDC file query returned **448** files on 2025-10-01 and **443** are on disk.
The five are not a partial download of ours — they are gone from GDC:

    1b06673a-c03b-48e4-a220-e9195c515c9e    404 NOT FOUND, 5 attempts
    691e4c86-7ad9-40b7-9312-2d625472bf97    404 NOT FOUND, 5 attempts
    71106ad8-a2c4-4337-9682-c407356d39ef    404 NOT FOUND, 5 attempts
    1455889c-26f3-4d83-aa90-271d6dbdea65    404 NOT FOUND, 5 attempts
    79946f8f-c71f-477a-985b-0fb923937225    connection dropped, 3 attempts

The first four are what `retry_failed_downloads.py` was written for, and all four
failed again; the fifth failed differently and was never added to its list. A
404 from the GDC data endpoint means the file has been withdrawn from the data
release, so retrying cannot close the gap and a query run today will return a
different set. The four tables above are therefore a **frozen record of one
query**, and their md5s are the only way to tell whether a reader is holding
what the paper was computed on.

Measured 2026-10-04: no analysis result depends on which five went missing. Each
of the five is its patient's only aliquot, all five are `Primary Tumor`, and none
of the five patients is a column of the FPKM matrix.
`tcga_aliquot_provenance.py` refuses to run if that ever stops being true — if
an absent file's patient keeps aliquots of a *different* `sample_type`, the
tissue assigned to that column would be decided by the absence rather than by
the comparison, and it raises instead.

## The matrix is 410 columns wide, and that is deliberate to know about

443 patients, 448 aliquot files, 410 columns. `parse_fpkm_matrix.py` assembles
the matrix as `fpkm_dict[submitter_id] = …`, keyed by **patient**, so a patient
with more than one aliquot keeps whichever file the loop read last and the
assignment leaves no record of which.

That is why
`04_Revision_Analyses/05_R1.6_Spatial_Confounders/scripts/tcga_aliquot_provenance.py`
exists. It recovers each column's tissue by comparing the column against that
patient's own STAR files and keeping the closest match — a measurement, not a
label — because the label was overwritten here. Its answer: 393 of the 410
columns carry a Primary Tumor aliquot and **17 carry Solid Tissue Normal**, and
the sensitivity analysis beside it re-runs six published comparisons on the
tumour aliquots only. Its verdict, verbatim:

> every association keeps its direction and gains effect size, no P value loses
> an order of magnitude, and both survival comparisons stay null

It is careful to say what it does not do, and so is this: it does **not** re-run
the deconvolution. Feeding the cleaned matrix to BayesPrism would move printed
values in figures that are frozen against the submitted pages, which is the
author's ruling to make; the measurement exists to tell him which direction they
would move.

**Do not re-key this matrix by aliquot without reading that analysis.** The
published numbers were computed on the matrix as it stands.

## Who reads these in the release

    02_External/Bulk/TCGA_STAD/02_Raw_Data/Clinical/TCGA_STAD_clinical_data.tsv
      03_Final_Panels/03_Figure_3/03_C/create_tcga_ceacam_scatter.py
      03_Final_Panels/03_Figure_3/03_D/create_tcga_ceacam_scatter.py
      03_Final_Panels/Supplementary_New/S3_CEACAM_Metaprogram_Validation/S3_H/create_S3_H_tcga_survival.py

    …/FPKM/file_metadata.tsv
      04_Revision_Analyses/05_R1.6_Spatial_Confounders/scripts/tcga_aliquot_provenance.py

    …/FPKM/TCGA_STAD_gene_info.tsv
      upstream/BayesPrism_TCGA/scripts/step1_prepare_tcga_bulk.py
      upstream/BayesPrism_TCGA/scripts/step1_prepare_tcga_bulk_pinned.py

The ABSOLUTE purity calls these analyses also use are a separate download, from
the TCGA PanCanAtlas publication page rather than from the project endpoint, and
are documented in `01_Reviewer_Materials/public_data/PROVENANCE.txt`.
