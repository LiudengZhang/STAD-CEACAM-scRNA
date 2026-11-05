> **Upstream record.** Paths in this document refer to the upstream
> processing pipeline that produced the deposited intermediates; it is not
> part of this release's run path. The files under `upstream/` are included
> as a record of how the inputs were produced and are not executed by
> `./run` or `_run_all_panels.sh`.

# GraphST - where this code came from

Source: `upstream-pipeline/01_standardized-pipeline/11_Spatial_Analysis/02_Applications/App_GraphST_Deconvolution_14Types_CEACAM`

Produces: spot-level deconvolution matrices and spot_data.csv

Consumed by: Figure 3F, 3G, 3K, 3L, 3M

| file | md5 | bytes |
|---|---|---|
| `01_Reference_Preparation/scripts/01_prepare_stomach_reference_14types.py` | `3dfb5c1f1f1e94e68083988c9b5d557b` | 9857 |
| `02_GraphST_Analysis/scripts/01_run_graphst_per_sample.py` | `ea6f30b1169a19e75b84489c72c4ade3` | 16725 |
| `03_Extract_Results/scripts/01_extract_proportions_and_domains.py` | `9909ed4d2b581164eb48a80d2b7cb879` | 6684 |
| `04_Visualization/scripts/03_ceacam_immune_mixed_model.py` | `2d6fc060852bf340dcc6ace5c971b5ad` | 13591 |
| `04_Visualization/scripts/04_ceacam_immune_region_analysis.py` | `05ec1b69d2658191b219a165e940ae65` | 14974 |
| `04_Visualization/scripts/08_ceacam_tumor_centered_analysis.py` | `1610a6338d4c64f92dd04b64f799f904` | 20639 |

Copied verbatim; any change made for determinism is recorded in `run.sh`, never by editing these.
