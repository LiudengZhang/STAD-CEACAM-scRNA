#!/bin/bash
# Run the panel scripts for the main and supplementary figures, one at a
# time (each loads large h5ad files).
#
# Usage: bash _run_all_panels.sh [figure_number|supp|revision]
#   e.g. bash _run_all_panels.sh 4             # only Figure 4
#        bash _run_all_panels.sh supp          # Supplementary Figures S1-S10
#        bash _run_all_panels.sh revision      # revision analyses, then S1-S10
#        bash _run_all_panels.sh               # Figures 1-5 and S1-S10
#
#
# Panels only; page assembly is not part of this repository.
#
#
# Conda env: stad_ceacam (override with STAD_CONDA_ENV)
#

set -euo pipefail

# Set so that set iteration order, and with it Figure 4B, repeats.
export PYTHONHASHSEED=0

CONDA_CMD="conda run -n ${STAD_CONDA_ENV:-stad_ceacam} python"
BASE="$(cd "$(dirname "$0")" && pwd)"

PASS=0
FAIL=0
MISSING=0
FAILED_SCRIPTS=""
MISSING_SCRIPTS=""

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
        echo "MISSING: $script is named by this driver and is not in the tree"
        MISSING=$((MISSING + 1))
        MISSING_SCRIPTS="$MISSING_SCRIPTS\n  $script"
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


# FIGURE 1 (B, C; panel A is artwork and is not generated here)
run_fig1() {
    echo -e "\n=== FIGURE 1 (A-C) ==="
    # Figure 1A is artwork and is not part of this repository.
    echo "SKIPPED _panel_1A/build_figure.py - Figure 1A is the study-overview artwork and is not included in this repository"
    run_script "$BASE/01_Figure_1/01_B/create_stomach_umap.py"
    run_script "$BASE/01_Figure_1/01_C/create_stomach_stacked_bar.py"
}

# FIGURE 2 (panels A-N)
run_fig2() {
    echo -e "\n=== FIGURE 2 (A-Q) ==="
    run_script "$BASE/02_Figure_2/02_A/create_epithelial_umap.py"
    run_script "$BASE/02_Figure_2/02_B/create_ceacam5_umap.py"
    run_script "$BASE/02_Figure_2/02_C/create_ceacam6_umap.py"
    run_script "$BASE/02_Figure_2/02_D/create_ceacam_metacell_correlation.py"
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
    # 02_P is not run; that heatmap is Supplementary Figure S3 G, drawn by
    # the S3 scripts.
    echo "SKIPPED 02_Figure_2/02_P/create_stomach_jaccard_heatmap.py - this heatmap is not part of Figure 2 (panels A-N)"
}

# FIGURE 3 (14 panels: A-N)
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
}

# FIGURE 4 (8 panels: A-H)
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
}

# FIGURE 5 (panels A-N)
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
}



# --- revision target ---
# DISPATCH.csv lists the revision analyses in dependency order, one per row,
# with interpreter, environment, an optional enabling variable, and the reason
# for any unit that does not run.
DISPATCH="$BASE/../04_Revision_Analyses/DISPATCH.csv"

# A unit is a Python script or a shell script.
run_unit() {
    local interp="$1" script="$2"
    if [ "$interp" != bash ]; then
        run_script "$script"
        return
    fi
    if [ -f "$script" ]; then
        echo ">>> Running: $script"
        if (cd "$(dirname "$script")" && bash "$(basename "$script")" 2>&1); then
            PASS=$((PASS + 1))
        else
            FAIL=$((FAIL + 1))
            FAILED_SCRIPTS="$FAILED_SCRIPTS\n  $script"
        fi
    else
        echo "MISSING: $script is named by DISPATCH.csv and is not in the tree"
        MISSING=$((MISSING + 1))
        MISSING_SCRIPTS="$MISSING_SCRIPTS\n  $script"
    fi
}

run_revision() {
    if [ ! -f "$DISPATCH" ]; then
        echo "MISSING: $DISPATCH is the dispatch manifest and is not in the tree"
        MISSING=$((MISSING + 1))
        MISSING_SCRIPTS="$MISSING_SCRIPTS\n  $DISPATCH"
        return
    fi
    # Columns: order,module,unit,interpreter,env,condition,status,reason.
    # `reason` takes the rest of the line, so it may not contain a comma.
    local order module unit interp env cond status reason script
    while IFS=, read -r order module unit interp env cond status reason; do
        [ -n "$order" ] || continue
        script="$BASE/../04_Revision_Analyses/$module/$unit"
        if [ "$status" != run ]; then
            echo "SKIPPED $module/$unit - $reason"
            continue
        fi
        if [ -n "$cond" ] && [ "${!cond:-0}" != "1" ]; then
            echo "SKIPPED $module/$unit - $reason (set $cond=1 to run)"
            continue
        fi
        # Apply the env column only where it differs from the default, so an
        # environment chosen with STAD_CONDA_ENV for the whole run still applies.
        case "$env" in
            stad_numpy2) STAD_CONDA_ENV="${STAD_NUMPY2_ENV:-stad_numpy2}" \
                             run_unit "$interp" "$script" ;;
            ""|stad_ceacam) run_unit "$interp" "$script" ;;
            *) STAD_CONDA_ENV="$env" run_unit "$interp" "$script" ;;
        esac
    done < <(tail -n +2 "$DISPATCH")
    run_shipped_supp
}


# SUPPLEMENTARY FIGURES S1-S10
run_shipped_supp() {
    echo -e "\n=== SUPPLEMENTARY FIGURES S1-S10 (shipped) ==="
    # Panels drawn from the revision analyses: Supplementary_New/_drivers/.
    for s in "$BASE"/Supplementary_New/_drivers/draw_*.py; do
        [ -e "$s" ] || continue
        run_script "$s"
    done
    # Then every create_*.py under Supplementary_New/S<n>_*/, in figure order.
    for s in $(ls -d "$BASE"/Supplementary_New/S*_*/ \
               | sed 's#\(.*/S\([0-9]*\)_[^/]*/\)$#\2 \1#' | sort -n | cut -d" " -f2); do
        for c in "$s"S*_?/create_*.py "$s"create_*.py; do
            [ -e "$c" ] || continue
            run_script "$c"
        done
    done
}



# MAIN
if [ -n "${1:-}" ]; then
    case "$1" in
        1) run_fig1 ;;
        2) run_fig2 ;;
        3) run_fig3 ;;
        4) run_fig4 ;;
        5) run_fig5 ;;
        supp|s|S) run_shipped_supp ;;
        revision|rev) run_revision ;;
        *) echo "Usage: $0 [1|2|3|4|5|supp|revision]"; exit 1 ;;
    esac
else
    run_fig1
    run_fig2
    run_fig3
    run_fig4
    run_fig5
    run_shipped_supp
fi

echo ""
echo "========================================"
echo "  COMPLETE: $(date)"
echo "  PASS: $PASS"
echo "  FAIL: $FAIL"
echo "  MISSING: $MISSING"
if [ -n "$FAILED_SCRIPTS" ]; then
    echo "  Failed scripts:"
    echo -e "$FAILED_SCRIPTS"
fi
if [ -n "$MISSING_SCRIPTS" ]; then
    echo "  Missing scripts:"
    echo -e "$MISSING_SCRIPTS"
fi
echo "========================================"

# run_script() does not stop on a failing script, so that one broken panel
# does not stop the rest; the exit status reports any failed or missing script.
#
if [ "$FAIL" -gt 0 ] || [ "$MISSING" -gt 0 ]; then
    [ "$FAIL" -gt 0 ] && echo "  $FAIL script(s) failed; exiting non-zero." >&2
    [ "$MISSING" -gt 0 ] && echo "  $MISSING script(s) named by this driver are not in the tree; exiting non-zero." >&2
    exit 1
fi
