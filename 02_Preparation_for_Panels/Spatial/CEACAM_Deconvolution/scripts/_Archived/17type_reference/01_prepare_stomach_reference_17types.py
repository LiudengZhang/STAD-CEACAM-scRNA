#!/usr/bin/env python3
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
"""
Prepare Single-Cell Reference for GraphST Deconvolution
17 Cell Types - MoMac Split

Takes 15-type reference and splits MoMac into 3 subtypes:
  - Monocytes (C1, C4)
  - IL2_high_Macrophage (C3)
  - IL2_low_Macrophage (C0, C2, C5, C6)

Author: Generated for Round 4 Analysis
Date: 2025-11-26
"""

import os
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1

# === CONFIGURATION ===
REF_15TYPES = "/path/to/home/Project_4_Gastric_Cancer/01.1_GraphST_15Types/00_Data/reference/stomach_15types_reference.h5ad"
OUTPUT_PATH = "/path/to/home/Project_4_Gastric_Cancer/01.2_GraphST_17Types/00_Data/reference/stomach_17types_reference.h5ad"

# Expected 17 cell types
CELL_TYPES_17 = [
    # Non-epithelial (14)
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "NK cells",
    "DC cells",
    "Endothelial cells",
    "Fibroblast",
    "Mast cells",
    "Neutrophils",
    "Pericyte",
    "Plasma cells",
    "Monocytes",
    "IL2_high_Macrophage",
    "IL2_low_Macrophage",
    # Epithelial (3)
    "Normal Epi",
    "Cancer F10-high",
    "Cancer F10-low"
]

# MoMac cluster mapping
MOMAC_MAPPING = {
    "C1": "Monocytes",
    "C4": "Monocytes",
    "C3": "IL2_high_Macrophage",
    "C0": "IL2_low_Macrophage",
    "C2": "IL2_low_Macrophage",
    "C5": "IL2_low_Macrophage",
    "C6": "IL2_low_Macrophage"
}


def main():
    print("="*70)
    print("Reference Preparation for GraphST - 17 Cell Types (MoMac Split)")
    print("="*70)

    # 1. Load 15-type reference
    print(f"\n1. Loading 15-type reference...")
    print(f"   {REF_15TYPES}")
    adata = sc.read_h5ad(REF_15TYPES)
    print(f"   Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Show existing cell types
    print(f"\n   Existing cell types (cell_type_15):")
    for ct, count in adata.obs['cell_type_15'].value_counts().items():
        print(f"     {ct}: {count:,}")

    # 2. Create cell_type_17 column
    print(f"\n2. Creating cell_type_17 column...")

    # Start with copy of cell_type_15
    adata.obs['cell_type_17'] = adata.obs['cell_type_15'].astype(str)

    # Find MoMac cells
    momac_mask = adata.obs['cell_type_15'] == "Monocytes/Macrophages"
    momac_cells = adata.obs_names[momac_mask]
    print(f"   MoMac cells to remap: {len(momac_cells):,}")

    # Remap MoMac cells based on minor_cell_state prefix
    remap_counts = {"Monocytes": 0, "IL2_high_Macrophage": 0, "IL2_low_Macrophage": 0, "unknown": 0}

    for cell_id in momac_cells:
        minor_state = adata.obs.loc[cell_id, 'minor_cell_state']
        # Extract cluster prefix (e.g., "C0" from "C0_Mac_Classic_TREM2")
        cluster = minor_state.split("_")[0] if "_" in str(minor_state) else str(minor_state)

        if cluster in MOMAC_MAPPING:
            new_type = MOMAC_MAPPING[cluster]
            adata.obs.loc[cell_id, 'cell_type_17'] = new_type
            remap_counts[new_type] += 1
        else:
            # Unknown cluster - default to IL2_low_Macrophage
            adata.obs.loc[cell_id, 'cell_type_17'] = "IL2_low_Macrophage"
            remap_counts["unknown"] += 1

    print(f"\n   MoMac remapping results:")
    print(f"     Monocytes: {remap_counts['Monocytes']:,}")
    print(f"     IL2_high_Macrophage: {remap_counts['IL2_high_Macrophage']:,}")
    print(f"     IL2_low_Macrophage: {remap_counts['IL2_low_Macrophage']:,}")
    if remap_counts['unknown'] > 0:
        print(f"     Unknown (→ IL2_low): {remap_counts['unknown']:,}")

    # 3. Final distribution
    print("\n" + "="*70)
    print("Final cell type distribution (17 types):")
    print("="*70)
    counts = adata.obs['cell_type_17'].value_counts()
    for ct in CELL_TYPES_17:
        if ct in counts.index:
            print(f"  {ct}: {counts[ct]:,}")
    print(f"\n  TOTAL: {adata.n_obs:,} cells")

    # 4. Save
    print(f"\n4. Saving reference...")
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    adata.write(OUTPUT_PATH)

    file_size = os.path.getsize(OUTPUT_PATH) / 1e9
    print(f"   Saved: {OUTPUT_PATH}")
    print(f"   Size: {file_size:.2f} GB")

    print("\n" + "="*70)
    print("Reference Preparation Complete!")
    print("="*70)


if __name__ == "__main__":
    main()
