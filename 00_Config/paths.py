"""
Central path configuration for Round_5 project.
All scripts should import paths from this module.

Usage:
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent / "00_Config"))
    from paths import *
"""

import os
from pathlib import Path

# =============================================================================
# Project Root
# =============================================================================
PROJECT_ROOT = Path(__file__).parent.parent
# Input data are not bundled with the code. Set STAD_RAW_INPUTS to
# point at the deposited record; the default is the in-repository
# directory, which holds placeholders only.
RAW_INPUTS = Path(os.environ.get(
    "STAD_RAW_INPUTS", PROJECT_ROOT / "01_Raw_Inputs"))
# Intermediate results are deposited with the data, not with the code.
# STAD_PREPARED_INPUTS points at them; the default is the in-repository
# directory, which holds the scripts that produce them.
PREPARATION = Path(os.environ.get(
    "STAD_PREPARED_INPUTS", PROJECT_ROOT / "02_Preparation_for_Panels"))
FINAL_PANELS = PROJECT_ROOT / "03_Final_Panels"
# Defined below, in the revision block, as 05_Manuscript.

# =============================================================================
# Raw Inputs - H5AD Files
# =============================================================================
H5AD_DIR = RAW_INPUTS / "01_H5AD"

# Major cell types
B_CELLS_H5AD = H5AD_DIR / "B_cells.h5ad"
DC_CELLS_H5AD = H5AD_DIR / "DC_cells.h5ad"
ENDOTHELIAL_H5AD = H5AD_DIR / "Endothelial.h5ad"
EPITHELIAL_H5AD = H5AD_DIR / "Epithelial.h5ad"
FIBROBLAST_H5AD = H5AD_DIR / "Fibroblast.h5ad"
MOMAC_H5AD = H5AD_DIR / "MoMac.h5ad"
NEUTROPHILS_H5AD = H5AD_DIR / "Neutrophils.h5ad"
NK_CELLS_H5AD = H5AD_DIR / "NK_cells.h5ad"
PERICYTE_H5AD = H5AD_DIR / "Pericyte.h5ad"
TCD4_H5AD = H5AD_DIR / "TCD4.h5ad"
TCD8_H5AD = H5AD_DIR / "TCD8.h5ad"

# Full dataset (all cells)
FULL_DATASET_H5AD = H5AD_DIR / "full_dataset.h5ad"

# Specialized h5ad variants
EPITHELIAL_TUMOR_SCORED_H5AD = H5AD_DIR / "Epithelial_tumor_scored.h5ad"
EPITHELIAL_RAW_COUNTS_H5AD = H5AD_DIR / "epithelial_raw_counts_full.h5ad"

# =============================================================================
# Raw Inputs - External Data
# =============================================================================
EXTERNAL_DIR = RAW_INPUTS / "02_External"

# Bulk RNA-seq (PRJEB25780 - TIGER)
BULK_DIR = EXTERNAL_DIR / "Bulk"
TIGER_DIR = BULK_DIR / "PRJEB25780_Tiger"
TIGER_EXPR = TIGER_DIR / "expression_matrix.txt"
TIGER_META = TIGER_DIR / "metadata_with_purity.tsv"

# Bulk RNA-seq (TCGA-STAD)
TCGA_DIR = BULK_DIR / "TCGA_STAD"
TCGA_CLINICAL = TCGA_DIR / "02_Raw_Data" / "Clinical" / "TCGA_STAD_clinical_data.tsv"
TCGA_RAW_DIR = TCGA_DIR / "02_Raw_Data" / "Gene_Expression" / "FPKM" / "raw_files"

# External single-cell datasets
SINGLE_CELL_DIR = EXTERNAL_DIR / "Single_Cell"
KUMAR_DIR = SINGLE_CELL_DIR / "GSE183904_Kumar"
LINGHUA_DIR = SINGLE_CELL_DIR / "GSE239676_Linghua"

# External spatial datasets
SPATIAL_EXT_DIR = EXTERNAL_DIR / "Spatial"
KOREAN_GUT_DIR = SPATIAL_EXT_DIR / "GSE251950_Korean_Gut"
GSE246011_DIR = SPATIAL_EXT_DIR / "GSE246011"

# Raw Inputs - IHC Images
IHC_DIR = RAW_INPUTS / "02_IHC_Raw"
IHC_THUMBNAILS = IHC_DIR / "Thumbnails"

# =============================================================================
# Preparation - Analysis Results
# =============================================================================

# DEG results (MAST)
DEG_DIR = PREPARATION / "DEG"
DEG_PRE_DIR = DEG_DIR / "pre"
DEG_POST_DIR = DEG_DIR / "post"

# GSEA results
GSEA_DIR = PREPARATION / "GSEA"
GSEA_PRE_DIR = GSEA_DIR / "pre"
GSEA_POST_DIR = GSEA_DIR / "post"
# Figures 5F/5G were produced from the Welch t-test analysis. The
# MAST table beside it is the independent check added in revision;
# see NFKB_RANKINGS_README.txt in that directory.
NFKB_RANKINGS_CSV = GSEA_DIR / "nfkb_rankings_13types_ttest_reproduce.csv"

# SCENIC results (MoMac — pySCENIC: GRNBoost2 + cisTarget + AUCell)
SCENIC_DIR = PREPARATION / "SCENIC"
SCENIC_AUCELL = SCENIC_DIR / "aucell_matrix.h5ad"
SCENIC_REGULONS = SCENIC_DIR / "regulons.csv"

# NicheNet results (used by Panel 5C)
NICHENET_DIR = PREPARATION / "NicheNet"
LIGAND_ACTIVITIES_RAW = NICHENET_DIR / "c3mac_ligand_activities_raw.csv"
LIGAND_ACTIVITIES_AGG = NICHENET_DIR / "c3mac_ligand_activities_aggregated.csv"

# Spatial analysis results
SPATIAL_DIR = PREPARATION / "Spatial"
SPATIAL_CEACAM_DIR = SPATIAL_DIR / "CEACAM_Deconvolution"
SPATIAL_SPOT_DATA = SPATIAL_CEACAM_DIR / "spot_data.csv"
SPATIAL_MIXED_MODEL = SPATIAL_CEACAM_DIR / "mixed_model_results.csv"
SPATIAL_REGION_COMPARISON = SPATIAL_CEACAM_DIR / "region_comparison.csv"
SPATIAL_DECONV_MATRICES = SPATIAL_CEACAM_DIR / "deconvolution_matrices"
SPATIAL_DISTRIBUTION_RESULTS = SPATIAL_DIR / "spatial_distribution_results.csv"

# NMF results
NMF_DIR = PREPARATION / "NMF"
NMF_PER_SAMPLE = NMF_DIR / "per_sample_nmf"
NMF_INTERMEDIATE = NMF_DIR / "intermediate"

# BayesPrism results
TIGER_BAYESPRISM_DIR = PREPARATION / "BayesPrism"
TIGER_BAYESPRISM_EPI = TIGER_BAYESPRISM_DIR / "bayesprism_epithelial_expression.tsv"
TIGER_BAYESPRISM_FRAC = TIGER_BAYESPRISM_DIR / "bayesprism_fractions.tsv"

TCGA_BAYESPRISM_DIR = PREPARATION / "BayesPrism_TCGA"
TCGA_BAYESPRISM_EPI = TCGA_BAYESPRISM_DIR / "tcga_bayesprism_epithelial_expression.tsv"
TCGA_BULK_TUMOR_ONLY = TCGA_BAYESPRISM_DIR / "tcga_bulk_counts_tumor_only.tsv"

# =============================================================================
# Prepared Outputs (intermediate data ready for plotting)
# =============================================================================
PREP_OUTPUTS = PREPARATION / "outputs"

# Figure 4 prepared data
FIG4_PREP = PREP_OUTPUTS / "Figure_04"
FIG4_MODULE_PROPORTIONS = FIG4_PREP / "module_proportions.csv"
FIG4_CLINICAL_METADATA = FIG4_PREP / "clinical_metadata.csv"
FIG4_MOMAC_FOLDCHANGE = FIG4_PREP / "momac_foldchange_bootstrap_ci.csv"
FIG4_EXTERNAL_LINGHUA = FIG4_PREP / "external_validation" / "linghua_mac3_proportions.csv"
FIG4_EXTERNAL_KUMAR = FIG4_PREP / "external_validation" / "kumar_mac3_proportions.csv"

# =============================================================================
# Shared Code
# =============================================================================
SHARED_CODE = PROJECT_ROOT / "00_Config" / "shared"

# =============================================================================
# Helper function to get cell type h5ad by name
# =============================================================================
CELL_TYPE_H5AD = {
    "B_cells": B_CELLS_H5AD,
    "DC_cells": DC_CELLS_H5AD,
    "Endothelial": ENDOTHELIAL_H5AD,
    "Epithelial": EPITHELIAL_H5AD,
    "Fibroblast": FIBROBLAST_H5AD,
    "MoMac": MOMAC_H5AD,
    "Neutrophils": NEUTROPHILS_H5AD,
    "NK_cells": NK_CELLS_H5AD,
    "Pericyte": PERICYTE_H5AD,
    "TCD4": TCD4_H5AD,
    "TCD8": TCD8_H5AD,
}

def get_h5ad(cell_type: str) -> Path:
    """Get h5ad path for a cell type."""
    if cell_type not in CELL_TYPE_H5AD:
        raise ValueError(f"Unknown cell type: {cell_type}. Available: {list(CELL_TYPE_H5AD.keys())}")
    return CELL_TYPE_H5AD[cell_type]


# =============================================================================
# DEG cell type names (match MAST output filenames)
# =============================================================================
DEG_CELL_TYPES = [
    "B_cells", "DC_cells", "Endothelial_cells", "Epithelial", "Fibroblast",
    "Mast_cells", "MoMac", "Neutrophils", "NK_cells", "Pericyte",
    "Plasma_cells", "TCD4_cells", "TCD8_cells",
]

def get_deg(cell_type: str, comparison: str = "post") -> Path:
    """Get MAST DEG CSV path for a cell type and comparison."""
    d = DEG_PRE_DIR if comparison == "pre" else DEG_POST_DIR
    return d / f"{cell_type}_mast_deg.csv"

# Hallmark gene sets, shipped with the code. gseapy given the library
# NAME downloads it from Enrichr on each call, which makes the analysis
# depend on a remote service and records nothing about which release was
# used. See 00_Reference/README.txt.
HALLMARK_GMT = PROJECT_ROOT / "00_Reference" / "MSigDB_Hallmark_2020.gmt"

# =============================================================================
# Revision (CIR-26-0753-ET)
# =============================================================================
# Analyses added in response to the reviewers, one directory per reviewer point.
NEW_ANALYSES = PROJECT_ROOT / "04_Revision_Analyses"

# Panels regenerated or added for the revision. The two-sided main-figure panels
# replace their one-tailed originals in place, so this is the figure tree itself;
# the supplementary figures added in revision live under Supplementary_New/.
REVISED_PANELS = FINAL_PANELS

# Supplementary tables ST1-ST9, shipped with the code because several revision
# analyses read the cohort and signature definitions out of them.
MANUSCRIPT = PROJECT_ROOT / "05_Manuscript"

# Public reference data that are not part of the Zenodo record and are fetched
# from their own sources; see the README in that directory.
REVIEWER_MATERIALS = RAW_INPUTS

# Names the working-tree paths.py never defined although scripts in the release
# import them. Without these the SCENIC pipeline and the IHC panel fail on the
# import line.
SCENIC_RESULTS_DIR = SCENIC_DIR / "Results"
# cisTarget motif databases; download from https://resources.aertslab.org/cistarget/
SCENIC_DB_DIR = SCENIC_DIR / "databases"
IHC_COLOR_DECONV_CSV = PREPARATION / "IHC" / "ceacam_ihc_color_deconv_results.csv"

# CellPhoneDB v5 database, used by the CEACAM ligand-receptor analysis. It is
# not redistributed here; download the release zip from
# https://github.com/ventolab/cellphonedb-data and place it at this path.
CPDB_DB_ZIP = RAW_INPUTS / "02_External" / "CellPhoneDB" / "cellphonedb.zip"

# MP4/MP5 permutation results, produced by
# 02_Preparation_for_Panels/Metaprogram_Permutation/mp4_pre_r_permutation_analysis.py
# and read by the two-sided sweep and the MP direction analysis.
MP_PERMUTATION_DIR = PREPARATION / "Metaprogram_Permutation"
