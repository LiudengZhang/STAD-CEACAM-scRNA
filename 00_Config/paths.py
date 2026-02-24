"""
Central path configuration for the STAD-CEACAM project.
All scripts should import paths from this module.

Usage:
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent / "00_Config"))
    from paths import *
"""

from pathlib import Path

# =============================================================================
# Project Root
# =============================================================================
PROJECT_ROOT = Path(__file__).parent.parent
RAW_INPUTS = PROJECT_ROOT / "01_Raw_Inputs"
PREPARATION = PROJECT_ROOT / "02_Preparation_for_Panels"
FINAL_PANELS = PROJECT_ROOT / "03_Final_Panels"

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

# Raw Inputs - IHC Images
IHC_DIR = RAW_INPUTS / "03_IHC"
IHC_THUMBNAILS = IHC_DIR / "Thumbnails"

# IHC analysis results
IHC_ANALYSIS_DIR = PREPARATION / "IHC"
IHC_COLOR_DECONV_CSV = IHC_ANALYSIS_DIR / "ceacam_ihc_color_deconv_results.csv"

# Raw Inputs - Supplementary Tables
TABLES_DIR = RAW_INPUTS / "04_Tables"
ST1_CSV = TABLES_DIR / "ST1_patient_sample_characteristics.csv"

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
NFKB_RANKINGS_CSV = GSEA_DIR / "nfkb_rankings_13types_mast.csv"

# SCENIC results (MoMac — pySCENIC: GRNBoost2 + cisTarget + AUCell)
SCENIC_DIR = PREPARATION / "SCENIC"
SCENIC_RESULTS_DIR = SCENIC_DIR / "Results"
# cisTarget databases — download from https://resources.aertslab.org/cistarget/
SCENIC_DB_DIR = SCENIC_DIR / "databases"

# NicheNet results (used by Panel 5C)
NICHENET_DIR = PREPARATION / "NicheNet"
LIGAND_ACTIVITIES_RAW = NICHENET_DIR / "c3mac_ligand_activities_raw.csv"

# Spatial analysis results
SPATIAL_DIR = PREPARATION / "Spatial"
SPATIAL_CEACAM_DIR = SPATIAL_DIR / "CEACAM_Deconvolution"
SPATIAL_SPOT_DATA = SPATIAL_CEACAM_DIR / "spot_data.csv"
SPATIAL_REGION_COMPARISON = SPATIAL_CEACAM_DIR / "region_comparison.csv"

# NMF results
NMF_DIR = PREPARATION / "NMF"
NMF_PER_SAMPLE = NMF_DIR / "per_sample_nmf"
NMF_INTERMEDIATE = NMF_DIR / "intermediate"

# BayesPrism results
TIGER_BAYESPRISM_DIR = PREPARATION / "BayesPrism"
TIGER_BAYESPRISM_EPI = TIGER_BAYESPRISM_DIR / "bayesprism_epithelial_expression.tsv"

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

# Additional cell type h5ad files
MAST_CELLS_H5AD = H5AD_DIR / "Mast_cells.h5ad"
PLASMA_CELLS_H5AD = H5AD_DIR / "Plasma_cells.h5ad"

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
