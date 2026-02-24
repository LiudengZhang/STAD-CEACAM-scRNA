#!/bin/bash
# =============================================================================
# Run all panel-creation and assembly scripts for main figures & supplementaries
# Runs sequentially to avoid memory issues from parallel h5ad loading
#
# Usage: bash _run_all_panels.sh [figure_number]
#   e.g. bash _run_all_panels.sh 4      # only Figure 4
#        bash _run_all_panels.sh supp    # only supplementaries
#        bash _run_all_panels.sh         # all figures
#
# Conda env: stad_ceacam
# Updated: 2026-02-23 (rebuilt from live directory audit)
# =============================================================================

set -euo pipefail

CONDA_CMD="conda run -n stad_ceacam python"
BASE="$(cd "$(dirname "$0")" && pwd)"
PASS=0
FAIL=0
FAILED_SCRIPTS=""

run_script() {
    local script="$1"
    if [ -f "$script" ]; then
        echo ">>> Running: $script"
        if (cd "$(dirname "$script")" && $CONDA_CMD "$(basename "$script")" 2>&1); then
            PASS=$((PASS + 1))
        else
            FAIL=$((FAIL + 1))
            FAILED_SCRIPTS="$FAILED_SCRIPTS\n  $script"
        fi
    else
        echo "WARNING: $script not found, skipping"
    fi
}

echo "========================================"
echo "  Running ALL panel scripts"
echo "  Started: $(date)"
echo "========================================"


# ── Preparation Pipeline (uncomment to regenerate from scratch) ─────────────
# $CONDA_CMD "$BASE/../02_Preparation_for_Panels/DEG/scripts/01_run_mast_analysis.py"
# $CONDA_CMD "$BASE/../02_Preparation_for_Panels/GSEA/run_gsea_from_mast.py"
# $CONDA_CMD "$BASE/../02_Preparation_for_Panels/IHC/quantify_ceacam_ihc.py"
# $CONDA_CMD "$BASE/../02_Preparation_for_Panels/BayesPrism/step1_prepare_reference.py"
# $CONDA_CMD "$BASE/../02_Preparation_for_Panels/BayesPrism/step1b_reduce_genes.py"
# conda run -n r_bayesprism Rscript "$BASE/../02_Preparation_for_Panels/BayesPrism/step2_run_bayesprism.R"
# $CONDA_CMD "$BASE/../02_Preparation_for_Panels/BayesPrism/step3_plot_results.py"
# $CONDA_CMD "$BASE/../02_Preparation_for_Panels/Metaprogram_Permutation/mp4_pre_r_permutation_analysis.py"


# =====================================================================
# FIGURE 1 (Panels A=schematic, B, C)
# =====================================================================
run_fig1() {
    echo -e "\n=== FIGURE 1 (A-C) ==="
    run_script "$BASE/01_Figure_1/01_B/generate_stomach_umap.py"
    run_script "$BASE/01_Figure_1/01_C/generate_stomach_stacked_bar.py"
    run_script "$BASE/01_Figure_1/assemble_figure_1.py"
}

# =====================================================================
# FIGURE 2 (16 panels: A-Q, skipping E/F/L)
# =====================================================================
run_fig2() {
    echo -e "\n=== FIGURE 2 (A-Q) ==="
    run_script "$BASE/02_Figure_2/02_A/create_epithelial_umap.py"
    run_script "$BASE/02_Figure_2/02_B/create_ceacam5_umap.py"
    run_script "$BASE/02_Figure_2/02_C/create_ceacam6_umap.py"
    run_script "$BASE/02_Figure_2/02_C1/create_tumor_score_umap.py"
    run_script "$BASE/02_Figure_2/02_C2/create_cnv_score_umap.py"
    run_script "$BASE/02_Figure_2/02_D/create_ceacam_correlation.py"
    run_script "$BASE/02_Figure_2/02_E/create_c2_proportion_boxplot.py"
    run_script "$BASE/02_Figure_2/02_F/create_epithelial_milo_pre_rvsnr.py"
    run_script "$BASE/02_Figure_2/02_G/create_mp45_horizontal_boxplot.py"
    run_script "$BASE/02_Figure_2/02_H/create_mp45_gene_dotplot.py"
    run_script "$BASE/02_Figure_2/02_I/create_checkpoint_dotplot_flipped.py"
    run_script "$BASE/02_Figure_2/02_J/create_ceacam6_boxplot.py"
    run_script "$BASE/02_Figure_2/02_K/create_ceacam5_boxplot.py"
    run_script "$BASE/02_Figure_2/02_L/create_ihc_representative_2x2.py"
    run_script "$BASE/02_Figure_2/02_M/create_ihc_combined_boxplot.py"
    run_script "$BASE/02_Figure_2/02_N/create_ceacam6_prjeb25780_boxplot.py"
    run_script "$BASE/02_Figure_2/02_O/create_ceacam5_prjeb25780_boxplot.py"
    run_script "$BASE/02_Figure_2/02_P/create_stomach_jaccard_heatmap.py"
    run_script "$BASE/02_Figure_2/assemble_figure_2.py"
}

# =====================================================================
# FIGURE 3 (14 panels: A-N)
# =====================================================================
run_fig3() {
    echo -e "\n=== FIGURE 3 (A-N) ==="
    run_script "$BASE/03_Figure_3/03_A/create_ceacam_cd274_scatter.py"
    run_script "$BASE/03_Figure_3/03_B/create_ceacam_tex_scatter.py"
    run_script "$BASE/03_Figure_3/03_C/create_tcga_ceacam_scatter.py"
    run_script "$BASE/03_Figure_3/03_D/create_tcga_ceacam_scatter.py"
    run_script "$BASE/03_Figure_3/03_E/create_cd8_umap.py"
    run_script "$BASE/03_Figure_3/03_F/create_spatial_ceacam_ratio.py"
    run_script "$BASE/03_Figure_3/03_G/create_spatial_epi_density.py"
    run_script "$BASE/03_Figure_3/03_H/create_spatial_boxplots.py"
    run_script "$BASE/03_Figure_3/03_I/create_spatial_boxplot_stroma.py"
    run_script "$BASE/03_Figure_3/03_J/create_spatial_boxplot_immune.py"
    run_script "$BASE/03_Figure_3/03_K/create_spatial_dist_stroma.py"
    run_script "$BASE/03_Figure_3/03_L/create_spatial_dist_immune.py"
    run_script "$BASE/03_Figure_3/03_M/create_spatial_momac_fraction.py"
    run_script "$BASE/03_Figure_3/03_N/create_immune_recruitment_barplot.py"
    run_script "$BASE/03_Figure_3/assemble_figure_3.py"
}

# =====================================================================
# FIGURE 4 (8 panels: A-H)
# =====================================================================
run_fig4() {
    echo -e "\n=== FIGURE 4 (A-H) ==="
    run_script "$BASE/04_Figure_4/04_A/create_correlation_heatmap.py"
    run_script "$BASE/04_Figure_4/04_B/create_merged_networks.py"
    run_script "$BASE/04_Figure_4/04_C/create_module2_boxplot.py"
    run_script "$BASE/04_Figure_4/04_E/create_momac_umap.py"
    run_script "$BASE/04_Figure_4/04_E/create_momac_marker_dotplot.py"
    run_script "$BASE/04_Figure_4/04_F/create_density_panel.py"
    run_script "$BASE/04_Figure_4/04_G/create_forest_plot.py"
    run_script "$BASE/04_Figure_4/04_H/create_external_validation_boxplots.py"
    run_script "$BASE/04_Figure_4/assemble_figure_4.py"
}

# =====================================================================
# FIGURE 5 (22 panels: A-Q incl. GSEA Q1-Q4)
# =====================================================================
run_fig5() {
    echo -e "\n=== FIGURE 5 (A-Q) ==="
    run_script "$BASE/05_Figure_5/05_A/create_panel_a_momac_enrichment.py"
    run_script "$BASE/05_Figure_5/05_B/create_4pathway_umap_matched.py"
    run_script "$BASE/05_Figure_5/05_C/create_panel_c_ligand_violin.py"
    run_script "$BASE/05_Figure_5/05_D/create_gene_umaps.py"
    run_script "$BASE/05_Figure_5/05_TF/create_tf_4group_panels.py"
    run_script "$BASE/05_Figure_5/05_H/create_cytokine_dotplot.py"
    run_script "$BASE/05_Figure_5/05_F/create_panel_f_radar.py"
    run_script "$BASE/05_Figure_5/05_G/create_panel_g_gsea_2types.py"
    run_script "$BASE/05_Figure_5/05_I/create_cd274_momac_boxplot.py"
    run_script "$BASE/05_Figure_5/05_J/create_cd274_epithelial_boxplot.py"
    run_script "$BASE/05_Figure_5/05_K/create_cd274_fibroblast_boxplot.py"
    run_script "$BASE/05_Figure_5/05_DC_CD274/create_cd274_dc_boxplot.py"
    run_script "$BASE/05_Figure_5/05_O/create_tex_nfkb_scatter_cd8.py"
    run_script "$BASE/05_Figure_5/05_IL6_CD4/create_il6_stat3_cd4t_boxplot.py"
    run_script "$BASE/05_Figure_5/05_Q/create_th17_nfkb_scatter_cd4.py"
    # GSEA summary (compute then plot)
    run_script "$BASE/05_Figure_5/05_GSEA_Summary/run_gsea_momac.py"
    run_script "$BASE/05_Figure_5/05_GSEA_Summary/run_gsea_5celltypes.py"
    run_script "$BASE/05_Figure_5/05_GSEA_Summary/create_gsea_4panels.py"
    run_script "$BASE/05_Figure_5/assemble_figure_5.py"
}

# =====================================================================
# SUPPLEMENTARY FIGURES (S1-S6)
# =====================================================================
run_supp() {
    echo -e "\n=== SUPPLEMENTARY FIGURES ==="

    # S1 — QC & Annotation (8 panels A-H)
    echo -e "\n>>> S1: QC & Annotation <<<"
    run_script "$BASE/10_Supplementaries/S1_QC_Annotation/S1_A/create_S1_A_qc_umaps.py"
    run_script "$BASE/10_Supplementaries/S1_QC_Annotation/S1_B/create_S1_B_qc_metrics.py"
    run_script "$BASE/10_Supplementaries/S1_QC_Annotation/S1_C/create_S1_C_qc_non_stomach.py"
    run_script "$BASE/10_Supplementaries/S1_QC_Annotation/S1_D/create_S1_D_batch_integration.py"
    run_script "$BASE/10_Supplementaries/S1_QC_Annotation/S1_E/create_S1_E_doublet_score.py"
    run_script "$BASE/10_Supplementaries/S1_QC_Annotation/S1_F/create_S1_F_marker_dotplot_12types.py"
    run_script "$BASE/10_Supplementaries/S1_QC_Annotation/S1_G/create_S1_G_tcell_subtype_dotplot.py"
    run_script "$BASE/10_Supplementaries/S1_QC_Annotation/S1_H/create_S1_H_non_stomach_stacked_bar.py"
    run_script "$BASE/10_Supplementaries/S1_QC_Annotation/assemble_S1.py"

    # S2 — CEACAM & Metaprogram Validation (9 panels A-I + data_prep)
    echo -e "\n>>> S2: CEACAM & Metaprogram Validation <<<"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/data_prep/generate_liver_stomach_ceacam.py"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/S2_A/create_S2_A_epithelial_dotplot.py"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/S2_B/create_S2_B_tumor_score_umap.py"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/S2_C/create_S2_C_cnv_score_umap.py"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/S2_D/create_S2_D_survival.py"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/S2_E/create_S2_E_coexpression_prjeb.py"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/S2_F/create_S2_F_coexpression_tcga.py"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/S2_G/create_S2_G_metaprogram_heatmap.py"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/S2_H/create_S2_H_liver_ceacam_boxplot.py"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/S2_I/create_S2_I_liver_c2ceacam_proportion.py"
    run_script "$BASE/10_Supplementaries/S2_CEACAM_Metaprogram_Validation/assemble_S2.py"

    # S3 — CD8+ T Cells (4 panels A-D)
    echo -e "\n>>> S3: CD8+ T Cells <<<"
    run_script "$BASE/10_Supplementaries/S3_CD8_TCells/S3_A/create_S3_A_cd8_dotplot.py"
    run_script "$BASE/10_Supplementaries/S3_CD8_TCells/S3_B/create_S3_B_cd8_pdcd1_umap.py"
    run_script "$BASE/10_Supplementaries/S3_CD8_TCells/S3_C/create_S3_C_cd8_havcr2_umap.py"
    run_script "$BASE/10_Supplementaries/S3_CD8_TCells/S3_D/create_S3_D_cd8_tex_score_umap.py"
    run_script "$BASE/10_Supplementaries/S3_CD8_TCells/assemble_S3.py"

    # S4 — Spatial Validation (6 panels A-F)
    echo -e "\n>>> S4: Spatial Validation <<<"
    run_script "$BASE/10_Supplementaries/S4_Spatial_Validation/S4_A/create_S4_A_ceacam_ratio_all_samples.py"
    run_script "$BASE/10_Supplementaries/S4_Spatial_Validation/S4_B/create_S4_B_momac_fraction_all_samples.py"
    run_script "$BASE/10_Supplementaries/S4_Spatial_Validation/S4_C/create_S4_C_epithelial_density_all_samples.py"
    run_script "$BASE/10_Supplementaries/S4_Spatial_Validation/S4_D/create_S4_D_distance_to_stroma.py"
    run_script "$BASE/10_Supplementaries/S4_Spatial_Validation/S4_E/create_S4_E_distance_to_immune.py"
    run_script "$BASE/10_Supplementaries/S4_Spatial_Validation/S4_F/create_S4_F_stacked_deconvolution.py"
    run_script "$BASE/10_Supplementaries/S4_Spatial_Validation/assemble_S4.py"

    # S5 — Immune Modules (5 panels A-E, single generator)
    echo -e "\n>>> S5: Immune Modules <<<"
    run_script "$BASE/10_Supplementaries/S5_Immune_Modules/create_S5_module_boxplots.py"
    run_script "$BASE/10_Supplementaries/S5_Immune_Modules/assemble_S5.py"

    # S6 — CD274 Remaining Cell Types (9 panels A-I, single generator)
    echo -e "\n>>> S6: CD274 Remaining Cell Types <<<"
    run_script "$BASE/10_Supplementaries/S6_CD274_Remaining/create_S6_cd274_all_celltypes.py"
    run_script "$BASE/10_Supplementaries/S6_CD274_Remaining/assemble_S6.py"
}


# =====================================================================
# MAIN
# =====================================================================
if [ -n "${1:-}" ]; then
    case "$1" in
        1) run_fig1 ;;
        2) run_fig2 ;;
        3) run_fig3 ;;
        4) run_fig4 ;;
        5) run_fig5 ;;
        supp|s|S) run_supp ;;
        *) echo "Usage: $0 [1|2|3|4|5|supp]"; exit 1 ;;
    esac
else
    run_fig1
    run_fig2
    run_fig3
    run_fig4
    run_fig5
    run_supp
fi

echo ""
echo "========================================"
echo "  COMPLETE: $(date)"
echo "  PASS: $PASS"
echo "  FAIL: $FAIL"
if [ -n "$FAILED_SCRIPTS" ]; then
    echo "  Failed scripts:"
    echo -e "$FAILED_SCRIPTS"
fi
echo "========================================"
