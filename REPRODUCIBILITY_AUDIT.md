# Reproducibility audit

Audit date: 2026-09-17

This release was checked from a prospective clean checkout against the inputs
described in the README and deposited at Zenodo (DOI
`10.5281/zenodo.18737073`). The audit covered import and path resolution,
manuscript-number verification, panel provenance, figure rendering, release
synchronization, archive integrity, and publication-safety checks.

## Results

| Check | Result |
|---|---:|
| Release self-check | 302 Python files passed |
| Clean-checkout dependency and path audit | 216 runnable Python files passed |
| Default main workflow | 98 scripts passed; 0 failed |
| Manuscript numerical checks | 183 checks passed |
| Manuscript/table text checks | 38 checks passed |
| DOCX revision and embedded-figure integrity | Passed |
| Supplementary Figure S1–S10 rebuild | 37 scripts passed; 0 failed |
| Supplementary page comparison | 10 of 10 pixel-identical |
| Panel geometry and content sweep | 62 passed; 1 declared Figure 4A exemption; 0 failed |
| Main and supplementary page sweep | 5 of 5 main and 10 of 10 supplementary pages passed |
| Frozen source baseline | 882 analysis, 4 panel-written, 69 panel-data, and 348 drawing files matched |
| Submission archive | All members matched their sources; reviewer commands passed |
| Publication-safety scan | 2,696 text files scanned; no fatal findings |
| Data-deposit accession scan | 1,383 files scanned; no pathology accessions found |

The supplementary rebuild left cached data hashes unchanged. All ten rebuilt
pages matched the approved PDFs in pixels, text geometry, fonts, page size, and
bounding boxes. Figure S9 retained panels A–E; panels A and B contain no section
counts, panel C has clear annotation spacing, and panels D and E match their
source analysis axes.

Figure S10A uses the versioned lineage-score table. Figure S10B installs its
versioned prepared panel because its scaled expression values require the
feature-selected source matrix used for the analysis; the public deposited
MoMac object is feature-complete. The standard build therefore reproduces the
approved S10 page exactly.

The fresh default workflow regenerated 63 individual main panels: 62 were
pixel-identical, and Figure 4B differed only by one color-channel unit. Its 56
common CSV outputs were hash-identical.

## Workflow scope

Run the final revision analyses and Supplementary Figures S1–S10 with:

```bash
bash 03_Final_Panels/_run_all_panels.sh revision
```

Run the main panel workflow with:

```bash
bash 03_Final_Panels/_run_all_panels.sh
```

The revised main-figure pages use the published slot geometry recorded in the
repository; this release reproduces their individual panels. Supplementary
Figures S1–S10 are assembled in full.

Clinical response labels are treated as fixed. The public cohort cache contains
only the six fields required by released consumers. Private clinical-source
recomputation and internal method-development audits are outside the default
driver. CellTypist and CellPhoneDB supporting analyses remain available through
`STAD_RUN_OPTIONAL=1`; their deposited results are not required to assemble
S1–S10.

Archived upstream pipelines document how deposited intermediates were produced.
The default workflows consume the versioned deposited intermediates, as stated
in the README.
