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
#
# THIS IS THE SOURCE OF THE DEPOSITED DRIVER.
# The driver is maintained here and nowhere else. Keeping the only corrected
# copy inside 05_Code_Release/github_repo/ would make a generated tree the
# master of a hand-maintained file: an edit here would be reverted by the next
# build, and an edit there would be invisible to every check that reads the
# working tree. update_release.py copies this file into the release and then
# applies extend_driver(), prune_driver(), parameterize_driver_env(),
# fix_milo_env() and pin_driver_hash_seed() on top of it - all five are no-ops
# against this text, because this text already carries what they add.
# =============================================================================

set -euo pipefail

# Set so that set iteration order, and with it Figure 4B, repeats.
export PYTHONHASHSEED=0

CONDA_CMD="conda run -n ${STAD_CONDA_ENV:-stad_ceacam} python"
BASE="$(cd "$(dirname "$0")" && pwd)"
PASS=0
FAIL=0
FAILED_SCRIPTS=""

run_script() {
    local script="$1"
    if [ -f "$script" ]; then
        echo ">>> Running: $script"
        local env="${STAD_CONDA_ENV:-stad_ceacam}"
        if (cd "$(dirname "$script")" && conda run -n "$env" python "$(basename "$script")" 2>&1); then
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
# $CONDA_CMD "$BASE/../upstream/DEG/scripts/01_run_mast_analysis.py"
# $CONDA_CMD "$BASE/../upstream/GSEA/run_gsea_from_mast.py"
# $CONDA_CMD "$BASE/../upstream/IHC/quantify_ceacam_ihc.py"
# $CONDA_CMD "$BASE/../upstream/BayesPrism/step1_prepare_reference.py"
# $CONDA_CMD "$BASE/../upstream/BayesPrism/step1b_reduce_genes.py"
# conda run -n r_bayesprism Rscript "$BASE/../upstream/BayesPrism/step2_run_bayesprism.R"
# $CONDA_CMD "$BASE/../upstream/BayesPrism/step3_plot_results.py"
# $CONDA_CMD "$BASE/../upstream/Metaprogram_Permutation/mp4_pre_r_permutation_analysis.py"


# =====================================================================
# FIGURE 1 (Panels A=schematic, B, C)
# =====================================================================
run_fig1() {
    echo -e "\n=== FIGURE 1 (A-C) ==="
    # B and C are drawn at 1:1 by the create_*.py scripts (renamed from
    # generate_*.py on 2026-09-15 when the panels were redrawn; the old
    # names no longer exist). A is the study-overview artwork under
    # _panel_1A/. The shipped page - A across the top, B and C on one row
    # at one height, on the 171.10 x 229.31 mm page every main figure uses
    # since 2026-09-15 - is assembled from the slot tables
    # (panel_rects_v2.csv) by 03_Final_Panels/assemble_slotted.py, which
    # is not part of this release (RELEASE_GAPS.csv; the author's ruling).
    # assemble_figure_1.py beside these scripts is the pre-submission
    # assembler of the submitted 254 mm landscape page and is NOT run: it
    # would print a page the revised article does not.
    run_script "$BASE/01_Figure_1/01_B/create_stomach_umap.py"
    run_script "$BASE/01_Figure_1/01_C/create_stomach_stacked_bar.py"
}

# =====================================================================
# FIGURE 2 (16 panels: A-Q, skipping E/F/L)
# =====================================================================
run_fig2() {
    echo -e "\n=== FIGURE 2 (A-Q) ==="
    run_script "$BASE/02_Figure_2/02_A/create_epithelial_umap.py"
    run_script "$BASE/02_Figure_2/02_B/create_ceacam5_umap.py"
    run_script "$BASE/02_Figure_2/02_C/create_ceacam6_umap.py"
    run_script "$BASE/02_Figure_2/02_D/create_ceacam_correlation.py"
    run_script "$BASE/02_Figure_2/02_E/create_c2_proportion_boxplot.py"
    # Milo needs its own environment; see the Dockerfile.
    STAD_CONDA_ENV="${STAD_MILO_ENV:-pertpy_milo}" \
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


# --- revision target ---
# Added for CIR-26-0753-ET: the analyses and supplementary figures produced in
# response to the reviewers. Scripts whose name starts with an underscore are
# build tools rather than reproduction steps and are skipped.
run_revision() {
    for s in "$BASE"/../04_Revision_Analyses/*/scripts/*.py; do
        [ -e "$s" ] || continue
        case "$(basename "$s")" in _*) continue ;; esac
        case "$s" in
            */01_R1.3_Cohort_Pairing/*)
                echo "SKIPPED private clinical-audit/reference module 01_R1.3_Cohort_Pairing (de-identified outputs are deposited)"; continue ;;
            */14_MAST_Specification/*|*/16_GSEA_Metric_Sensitivity/*|*/17_NFkB_Claim_Ledger/*)
                echo "SKIPPED internal audit module $(basename "$(dirname "$(dirname "$s")")")"; continue ;;
            */15_Pseudobulk_Sample_Level/*)
                echo "SKIPPED full-workflow reference module 15_Pseudobulk_Sample_Level (use its run_all.sh)"; continue ;;
        esac
        # CellTypist and pyDESeq2 require numpy>=2; the main environment is
        # pinned to numpy 1.23.5, so they run in their own. See the Dockerfile.
        case "$(basename "$s")" in
            celltypist_annotation.py)
                if [ "${STAD_RUN_OPTIONAL:-0}" = "1" ]; then
                    STAD_CONDA_ENV="${STAD_NUMPY2_ENV:-stad_numpy2}" run_script "$s"
                else
                    echo "SKIPPED optional $(basename "$s") (set STAD_RUN_OPTIONAL=1 to run)"
                fi ;;
            ceacam_family_and_interactions.py)
                if [ "${STAD_RUN_OPTIONAL:-0}" = "1" ]; then
                    run_script "$s"
                else
                    echo "SKIPPED optional $(basename "$s") (set STAD_RUN_OPTIONAL=1 to run)"
                fi ;;
            pin_gene_sets.py)
                echo "SKIPPED preparation-only pin_gene_sets.py (the pinned GMT is deposited)" ;;
            rescore_mast_gsea.py)
                echo "SKIPPED audit/preparation-only rescore_mast_gsea.py (the adopted sound13 outputs are deposited)" ;;
            momac_lineage.py)
                if [ "${STAD_RUN_PREPARATION:-0}" = "1" ]; then
                    run_script "$s"
                else
                    echo "SKIPPED preparation/reference-only momac_lineage.py (the approved lineage-score tables are deposited; set STAD_RUN_PREPARATION=1 to recompute from a matching source object)"
                fi ;;
            spatial_positive_evidence.py)
                echo "SKIPPED audit/reference-only spatial_positive_evidence.py (its internal S11 output is not submitted)" ;;
            *)  run_script "$s" ;;
        esac
    done
    # The S7-S9 panels are drawn by the drivers under Supplementary_New/
    # _drivers/, not by create_*.py scripts in the panel directories: those
    # directories hold panel outputs (.svg/.pdf/.png) and no code at all, so
    # copy_tree(), which ships code only, never carries them into the release
    # and the old `*/*/create_*.py` glob matched nothing there.
    #
    # Measured in the reviewer simulation of 2026-09-10: the glob ran zero
    # scripts, so assemble_new_supplementaries.py found no panel SVG and died
    # with "no panel SVG in Supplementary_New/S7_CEACAM5_vs_CEACAM6/S7_A".
    # Three of the paper's nine supplementary figures could not be rebuilt from
    # the capsule.
    #
    # draw_*.py writes to Supplementary_New/<FIG>/<panel>/ through
    # _driver_base.save(), where <FIG> is exactly the key
    # assemble_new_supplementaries.FIGURES uses, so the drawing and the
    # assembly meet where the assembler already looks.
    for s in "$BASE"/Supplementary_New/_drivers/draw_*.py; do
        [ -e "$s" ] || continue
        run_script "$s"
    done
    # S1-S7, since 2026-09-15 (evening): every panel of the six submitted
    # pages is redrawn at print size by a create_*.py in its own panel
    # directory (the immune-module and PD-L1 figures by one script at the
    # figure level), reading the deposited h5ads and tables. They must run
    # before the assembler too. Since 2026-09-16 the directories carry the
    # printed numbers S1-S10 (S1 split into S1 and S2, S2-S9 renumbered
    # S3-S10), so the glob takes every S*_ directory rather than a fixed
    # range; the drivers above draw the S8-S10 panels.
    for s in $(ls -d "$BASE"/Supplementary_New/S*_*/ \
               | sed 's#\(.*/S\([0-9]*\)_[^/]*/\)$#\2 \1#' | sort -n | cut -d" " -f2); do
        for c in "$s"S*_?/create_*.py "$s"create_*.py; do
            [ -e "$c" ] || continue
            run_script "$c"
        done
    done
    run_script "$BASE/Supplementary_New/assemble_new_supplementaries.py"
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
        revision|rev) run_revision ;;
        *) echo "Usage: $0 [1|2|3|4|5|supp|revision]"; exit 1 ;;
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

# run_script() deliberately does not abort on a failing script - one broken
# panel must not take the other hundred and fifty down with it - but that
# tolerance must stop here. If the driver printed its FAIL count and still
# exited 0, ./run would exit 0 too, and a capsule in which twenty-six scripts
# had died would report success with nothing downstream able to see it. The
# count is the exit status, and code/run propagates it.
if [ "$FAIL" -gt 0 ]; then
    echo "  $FAIL script(s) failed; exiting non-zero." >&2
    exit 1
fi
