# GraphST - where this code came from

Source: `Round_4/01_Round_4.2_Standardized_Pipeline/11_Spatial_Analysis/02_Applications/App_GraphST_Deconvolution_14Types_CEACAM`

Produces: spot-level deconvolution matrices and spot_data.csv

Consumed by: Figure 3F, 3G, 3K, 3L, 3M

| file | md5 | bytes |
|---|---|---|
| `01_Reference_Preparation/scripts/01_prepare_stomach_reference_14types.py` | `3dfb5c1f1f1e94e68083988c9b5d557b` | 9857 |
| `01_Reference_Preparation/scripts/01_prepare_stomach_reference_14types_metacell.py` | `55a5eb95fa4920be4b2b692d72577242` | 12393 |
| `02_GraphST_Analysis/scripts/01_run_graphst_per_sample.py` | `ea6f30b1169a19e75b84489c72c4ade3` | 16725 |
| `02_GraphST_Analysis/scripts/02_batch_process_all_samples.py` | `db919d4671dd265acb9cfa614dc5e59f` | 5507 |
| `03_Extract_Results/scripts/01_extract_proportions_and_domains.py` | `9909ed4d2b581164eb48a80d2b7cb879` | 6684 |
| `04_Visualization/scripts/01_visualize_graphst_results.py` | `4855e855e860ce83c8fd970bd32c01a0` | 9748 |
| `04_Visualization/scripts/02_ceacam_vs_immune_binned.py` | `39ec6a72f2f2a36b81b5f027a4204983` | 9800 |
| `04_Visualization/scripts/02_plot_epithelial_deltas.py` | `15647dc94bfd0c110855ab8800c901d1` | 2846 |
| `04_Visualization/scripts/03_ceacam_immune_mixed_model.py` | `2d6fc060852bf340dcc6ace5c971b5ad` | 13591 |
| `04_Visualization/scripts/03_plot_epithelial_deltas_spatial.py` | `aab9f96d3519d6790a2105c3b242812d` | 7306 |
| `04_Visualization/scripts/04_ceacam_immune_region_analysis.py` | `05ec1b69d2658191b219a165e940ae65` | 14974 |
| `04_Visualization/scripts/04_delta_correlation.py` | `c80d4a390461f53ddcd175473a348e72` | 5144 |
| `04_Visualization/scripts/05_delta_correlation_by_sample.py` | `ce1a8cf835fd488cebcd269d3fa4bb03` | 4542 |
| `04_Visualization/scripts/06_delta_correlation_high_epi.py` | `d9ac6e4e5f3f6fad90fe52d0d1597ebc` | 5944 |
| `04_Visualization/scripts/07_domain_annotation.py` | `19ad003ea9f19c57cdfef258b4d0a63a` | 6605 |
| `04_Visualization/scripts/08_ceacam_tumor_centered_analysis.py` | `1610a6338d4c64f92dd04b64f799f904` | 20639 |
| `04_Visualization/scripts/10_tumor_score_histogram.py` | `726bb208794bf3bb796d941162f72a93` | 3247 |
| `_Archived/17type_reference/01_prepare_stomach_reference_17types.py` | `88ccff60941839bb4838e44b07f2b658` | 4263 |

Copied verbatim; any change made for determinism is recorded in `run.sh`, never by editing these.
