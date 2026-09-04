#!/usr/bin/env python3
# Paths below refer to the upstream Round_4 processing pipeline, which is
# not part of this release. This script is included as a record of how the
# input was produced; it is not called by _run_all_panels.sh.
"""
Prepare Single-Cell Reference for GraphST Deconvolution
14 Cell Types - CEACAM-based Epithelial Split

Sources: Individual annotated cell type files from 00_Data/01_Major_Cell_Types/
- Epithelial: Split tumor cells by CEACAM5/6 expression (median threshold)
- Non-epithelial: Use as-is with standardized naming

Cell Types (14):
  Immune (8): B cells, CD4+ T cells, CD8+ T cells, NK cells, DC cells,
              Mast cells, Neutrophils, Plasma cells
  Myeloid (1): Monocytes/Macrophages
  Stromal (3): Endothelial cells, Fibroblast, Pericyte
  Epithelial (2): Epi CEACAM-high, Epi CEACAM-low

Author: Generated for Round 4 Analysis
Date: 2025-12-19
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1

# === CONFIGURATION ===
DATA_DIR = "/path/to/Project_4_05232025/Round_4/04_Final_Panels/00_Set_Ups/00_Data/01_Major_Cell_Types"
OUTPUT_DIR = "/path/to/Project_4_05232025/Round_4/01_Round_4.2_Standardized_Pipeline/11_Spatial_Analysis/02_Applications/App_GraphST_Deconvolution_14Types_CEACAM/01_Reference_Preparation/results"
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "stomach_14types_reference.h5ad")

# Source files and their cell type mapping
CELL_TYPE_FILES = {
    "B_cells_standardized_annotated.h5ad": "B cells",
    "TCD4_cells_standardized_annotated.h5ad": "CD4+ T cells",
    "TCD8_annotated.h5ad": "CD8+ T cells",
    "nk_cells_integrated_sample_harmony.h5ad": "NK cells",
    "DC_cells_standardized_annotated.h5ad": "DC cells",
    "mast_cell_integrated.h5ad": "Mast cells",
    "neutrophil_filtered_v2.1_integrated.h5ad": "Neutrophils",
    "plasma_cell_integrated.h5ad": "Plasma cells",
    "MoMac_final_annotated.h5ad": "Monocytes/Macrophages",
    "Endothelial_cells_standardized_annotated.h5ad": "Endothelial cells",
    "Fibroblast_final_annotated.h5ad": "Fibroblast",
    "Pericyte_annotated_final.h5ad": "Pericyte",
}

# Epithelial file (special handling)
EPITHELIAL_FILE = "Epithelial_functionally_annotated.h5ad"

# CEACAM genes for split
CEACAM_GENES = ["CEACAM5", "CEACAM6"]

# Expected 14 cell types
CELL_TYPES_14 = [
    "B cells", "CD4+ T cells", "CD8+ T cells", "NK cells", "DC cells",
    "Mast cells", "Neutrophils", "Plasma cells",
    "Monocytes/Macrophages",
    "Endothelial cells", "Fibroblast", "Pericyte",
    "Epi CEACAM-high", "Epi CEACAM-low"
]


def load_and_prepare_nonepithelial(data_dir, file_mapping):
    """Load non-epithelial cell types and assign cell_type_14 labels."""
    print("\n" + "="*70)
    print("Loading non-epithelial cell types...")
    print("="*70)

    adata_list = []

    for filename, cell_type in file_mapping.items():
        filepath = os.path.join(data_dir, filename)
        print(f"\n  Loading {cell_type}...")
        print(f"    File: {filename}")

        adata = sc.read_h5ad(filepath)
        print(f"    Cells: {adata.n_obs:,}")

        # Use raw data for consistent gene space
        if adata.raw is not None:
            adata_raw = sc.AnnData(
                X=adata.raw.X,
                obs=adata.obs.copy(),
                var=adata.raw.var.copy()
            )
        else:
            adata_raw = adata.copy()

        # Assign cell type
        adata_raw.obs['cell_type_14'] = cell_type

        adata_list.append(adata_raw)

    return adata_list


def load_and_prepare_epithelial(data_dir, ceacam_genes):
    """Load epithelial cells and split by CEACAM expression."""
    print("\n" + "="*70)
    print("Loading and processing epithelial cells...")
    print("="*70)

    filepath = os.path.join(data_dir, EPITHELIAL_FILE)
    print(f"\n  File: {EPITHELIAL_FILE}")

    adata = sc.read_h5ad(filepath)
    print(f"  Total epithelial cells: {adata.n_obs:,}")

    # Show existing classification
    print(f"\n  Existing major_cell_type_GraphST distribution:")
    for ct, count in adata.obs['major_cell_type_GraphST'].value_counts().items():
        print(f"    {ct}: {count:,}")

    # Filter to tumor cells only (exclude Normal and Not Selected)
    tumor_mask = adata.obs['major_cell_type_GraphST'].isin(['Tumor F10-high', 'Tumor F10-low'])
    adata_tumor = adata[tumor_mask].copy()
    n_tumor = adata_tumor.n_obs
    print(f"\n  Selected tumor cells: {n_tumor:,}")

    # Get CEACAM expression from raw
    raw = adata_tumor.raw
    if raw is None:
        raise ValueError("No raw data found in epithelial file!")

    # Check CEACAM genes
    available_genes = [g for g in ceacam_genes if g in raw.var_names]
    if len(available_genes) == 0:
        raise ValueError(f"No CEACAM genes found! Checked: {ceacam_genes}")
    print(f"\n  Using CEACAM genes: {available_genes}")

    # Extract CEACAM expression
    gene_indices = [list(raw.var_names).index(g) for g in available_genes]

    if hasattr(raw.X, 'toarray'):
        expr_matrix = raw.X[:, gene_indices].toarray()
    else:
        expr_matrix = np.array(raw.X[:, gene_indices])

    # Mean CEACAM expression per cell
    mean_ceacam = np.mean(expr_matrix, axis=1)

    print(f"\n  CEACAM expression statistics:")
    print(f"    Min: {np.min(mean_ceacam):.4f}")
    print(f"    Max: {np.max(mean_ceacam):.4f}")
    print(f"    Median: {np.median(mean_ceacam):.4f}")
    print(f"    Mean: {np.mean(mean_ceacam):.4f}")
    print(f"    % expressing (>0): {100 * np.mean(mean_ceacam > 0):.1f}%")

    # Split by median threshold
    threshold = np.median(mean_ceacam)
    ceacam_high_mask = mean_ceacam > threshold

    n_high = np.sum(ceacam_high_mask)
    n_low = n_tumor - n_high

    print(f"\n  CEACAM split (median threshold = {threshold:.4f}):")
    print(f"    Epi CEACAM-high: {n_high:,} ({100*n_high/n_tumor:.1f}%)")
    print(f"    Epi CEACAM-low: {n_low:,} ({100*n_low/n_tumor:.1f}%)")

    # Verify inverse relationship with F10
    print(f"\n  Validation - CEACAM vs F10 relationship:")
    f10_high_mask = adata_tumor.obs['major_cell_type_GraphST'] == 'Tumor F10-high'
    f10_high_ceacam_high_pct = 100 * np.mean(ceacam_high_mask[f10_high_mask])
    f10_low_ceacam_high_pct = 100 * np.mean(ceacam_high_mask[~f10_high_mask])
    print(f"    F10-high cells that are CEACAM-high: {f10_high_ceacam_high_pct:.1f}%")
    print(f"    F10-low cells that are CEACAM-high: {f10_low_ceacam_high_pct:.1f}%")

    # Create AnnData with raw expression
    adata_raw = sc.AnnData(
        X=raw.X,
        obs=adata_tumor.obs.copy(),
        var=raw.var.copy()
    )

    # Assign CEACAM-based labels
    adata_raw.obs['cell_type_14'] = 'Epi CEACAM-low'
    adata_raw.obs.loc[adata_raw.obs_names[ceacam_high_mask], 'cell_type_14'] = 'Epi CEACAM-high'

    # Store CEACAM expression for reference
    adata_raw.obs['mean_ceacam_expression'] = mean_ceacam

    return adata_raw


def combine_and_save(adata_list, output_path):
    """Combine all cell types and save reference."""
    print("\n" + "="*70)
    print("Combining all cell types...")
    print("="*70)

    # Find common genes across all datasets
    print("\n  Finding common genes...")
    gene_sets = [set(adata.var_names) for adata in adata_list]
    common_genes = gene_sets[0]
    for gs in gene_sets[1:]:
        common_genes = common_genes.intersection(gs)

    common_genes = sorted(list(common_genes))
    print(f"  Common genes: {len(common_genes):,}")

    # Subset to common genes
    adata_subset = []
    for adata in adata_list:
        adata_sub = adata[:, common_genes].copy()
        adata_subset.append(adata_sub)

    # Concatenate
    print("\n  Concatenating...")
    adata_combined = sc.concat(adata_subset, join='outer')

    print(f"  Combined shape: {adata_combined.shape}")

    # Final distribution
    print("\n" + "="*70)
    print("Final cell type distribution (14 types):")
    print("="*70)

    counts = adata_combined.obs['cell_type_14'].value_counts()
    for ct in CELL_TYPES_14:
        if ct in counts.index:
            print(f"  {ct}: {counts[ct]:,}")
        else:
            print(f"  {ct}: 0 (MISSING!)")

    print(f"\n  TOTAL: {adata_combined.n_obs:,} cells")

    # Verify all types present
    missing = [ct for ct in CELL_TYPES_14 if ct not in counts.index]
    if missing:
        print(f"\n  WARNING: Missing cell types: {missing}")

    # Clean obs columns for h5ad compatibility
    print("\n  Cleaning obs columns for h5ad compatibility...")
    cols_to_keep = ['cell_type_14']  # Only keep essential columns
    optional_cols = ['sample', 'Sample ID', 'Patient ID', 'Response', 'mean_ceacam_expression']
    for col in optional_cols:
        if col in adata_combined.obs.columns:
            cols_to_keep.append(col)

    # Keep only selected columns
    adata_combined.obs = adata_combined.obs[cols_to_keep].copy()

    # Ensure cell_type_14 is string type
    adata_combined.obs['cell_type_14'] = adata_combined.obs['cell_type_14'].astype(str)

    # Save
    print(f"\n  Saving to: {output_path}")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata_combined.write(output_path)

    file_size = os.path.getsize(output_path) / 1e9
    print(f"  Size: {file_size:.2f} GB")

    return adata_combined


def main():
    print("="*70)
    print("Reference Preparation for GraphST - 14 Cell Types (CEACAM Split)")
    print("="*70)
    print(f"\nSource directory: {DATA_DIR}")
    print(f"Output: {OUTPUT_PATH}")

    # Load non-epithelial cell types
    adata_list = load_and_prepare_nonepithelial(DATA_DIR, CELL_TYPE_FILES)

    # Load and process epithelial cells
    adata_epi = load_and_prepare_epithelial(DATA_DIR, CEACAM_GENES)
    adata_list.append(adata_epi)

    # Combine and save
    adata_combined = combine_and_save(adata_list, OUTPUT_PATH)

    print("\n" + "="*70)
    print("Reference Preparation Complete!")
    print("="*70)


if __name__ == "__main__":
    main()
