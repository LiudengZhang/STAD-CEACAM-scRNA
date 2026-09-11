> **Upstream record.** Paths in this document refer to the upstream
> processing pipeline that produced the deposited intermediates; it is not
> part of this release's run path. The files under `upstream/` are included
> as a record of how the inputs were produced and are not executed by
> `./run` or `_run_all_panels.sh`.

# NMF - where this code came from

Source: `upstream-pipeline/01_standardized-pipeline/01.1_meta_program_epithelial`

Produces: per-sample NMF programs and the epithelial metaprograms

Consumed by: Figure 2G, 2H

| file | md5 | bytes |
|---|---|---|
| `02_Scripts/01_Per_Sample_NMF/01_prepare_per_sample_data.py` | `4c895602c81767b636fd5fb0e8baa74d` | 3016 |
| `02_Scripts/01_Per_Sample_NMF/03_collect_nmf_programs.py` | `ca60245276ba9d71081a9662ef719e84` | 2526 |
| `02_Scripts/02_Robust_Filtering/01_within_sample_recurrence.py` | `617999689ae80c1b3cbef9d59b6398e1` | 2019 |
| `02_Scripts/02_Robust_Filtering/02_cross_sample_recurrence.py` | `3c76aa0cb33234379fedd529d9a23903` | 1921 |
| `02_Scripts/02_Robust_Filtering/03_remove_redundancy.py` | `cd9eb0bccc248c0f86d0ee34dea416e8` | 2748 |
| `02_Scripts/03_MetaProgram_Clustering/01_compute_jaccard_similarity.py` | `2ded4a53147013e7faba99e152e48eca` | 1744 |
| `02_Scripts/03_MetaProgram_Clustering/02_cluster_metaprograms.py` | `6232d7266a4226ea5d0d4c09b1c63bef` | 2789 |
| `02_Scripts/03_MetaProgram_Clustering/03_extract_gene_signatures.py` | `22b0c246d6187203946f7d618545b3af` | 3343 |
| `02_Scripts/04_Activity_Scoring/01_score_metaprogram_activity.py` | `18e2a604312b77932e53b764d9337def` | 2447 |
| `02_Scripts/04_Activity_Scoring/02_aggregate_by_sample.py` | `e48d87f30bd0807eb18a596788c64330` | 2061 |
| `02_Scripts/05_Visualization/01_k_selection_robustness.py` | `e807641c964d7d284b346386a6c306b9` | 2669 |
| `02_Scripts/05_Visualization/02_gene_heatmap.py` | `f06cea78a66ab18fab3b35c119701b89` | 2517 |
| `02_Scripts/05_Visualization/03_jaccard_heatmap.py` | `c83d5935c3cf65a3cfa4e6d0ccd9cfba` | 2469 |
| `02_Scripts/05_Visualization/04_dendrogram.py` | `39589d81ddaa01518bab51a3a66a2286` | 2029 |
| `02_Scripts/05_Visualization/05_umap_activity.py` | `3f9c3ada48b1c24c440b1165ac03edb8` | 2380 |
| `02_Scripts/05_Visualization/06_boxplot_pre_responder.py` | `cfe5f8e0cf136624bba518a20f462ebb` | 2795 |
| `02_Scripts/05_Visualization/07_boxplot_post_responder.py` | `46641953b26616a68ea36736374f2034` | 2806 |
| `02_Scripts/05_Visualization/08_pathway_dotplot.py` | `7c3924988fdaf0fb79fabc5d7e3fe3ce` | 4197 |
| `02_Scripts/01_Per_Sample_NMF/02_run_nmf_per_sample.R` | `b20cbc524f609ae912cf4384ed18e717` | 3030 |
| `00_Config/config.yaml` | `c22f398d062036782d4bcf3c31709571` | 3065 |

Copied verbatim; any change made for determinism is recorded in `run.sh`, never by editing these.
