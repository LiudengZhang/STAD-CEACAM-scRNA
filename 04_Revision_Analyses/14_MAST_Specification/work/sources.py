"""
Cell-type sources for the MAST specification audit.

Copied from 12_R1.8_DEG_Recompute/scripts/recompute_deg.py (CELL_SOURCES,
PHASES, assert_log1p_cp10k, load_subset, filter_genes), with exactly two
changes, both of which the task brief specifies:

  1. The eleven shared-pipeline types are read from 06_Clean_Data/01_H5AD/
     (the rebuilt, singly-normalised objects) rather than submission-tree.
  2. Neutrophils come from 06_Clean_Data/02_Rebuilt/Neutrophils_sound.h5ad,
     which is what 13_R1.8_Neutrophil_Rebuilt_Recompute uses.

Mast and plasma cells have no file in 01_H5AD and keep the upstream-pipeline objects the
original names. full_dataset.h5ad is not used for anything.
"""

from pathlib import Path
import sys

import numpy as np
import scanpy as sc
from scipy import sparse

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "00_Config"))
from paths import HALLMARK_GMT  # noqa: E402,F401

CLEAN = PROJECT_ROOT / "06_Clean_Data" / "01_H5AD"
REBUILT = PROJECT_ROOT / "06_Clean_Data" / "02_Rebuilt"
ROUND_4_MAJOR = (PROJECT_ROOT.parent / "upstream-pipeline" / "04_Final_Panels"
                 / "00_Set_Ups" / "00_Data" / "01_Major_Cell_Types")

# Which slot of each file holds the sound log1p CP10K matrix.
#
# This is not a fallback and it is not a guess. build_clean_h5ad.py promoted
# `.raw.X` to `.X` when it wrote 06_Clean_Data/01_H5AD/, and dropped `.raw`
# altogether, so those eleven files have no `.raw` to read: their `.X` *is*
# the matrix the rule "read .raw, never .X" points at. Checked, not assumed -
# 06_Clean_Data/01_H5AD/Pericyte.h5ad `.X` is bit-identical to
# submission-tree/01_Raw_Inputs/01_H5AD/Pericyte.h5ad `.raw.X` (max abs difference
# 0.0, same 5,205 cells and 57,058 genes in the same order), which is the file
# the archived sound run read.
#
# The rebuilt neutrophil object and the two upstream-pipeline objects do carry `.raw`,
# and it is read, exactly as 12_ and 13_ read it.
#
# Whichever slot is named, assert_log1p_cp10k runs on it before a single
# number is used, so a mislabelled slot stops that cell type rather than
# being tested.
SLOT = {
    "B_cells": "X", "DC_cells": "X", "Endothelial_cells": "X",
    "Epithelial": "X", "Fibroblast": "X", "MoMac": "X", "NK_cells": "X",
    "Pericyte": "X", "TCD4_cells": "X", "TCD8_cells": "X",
    "Neutrophils": "raw", "Mast_cells": "raw", "Plasma_cells": "raw",
}

CELL_SOURCES = {
    "B_cells":           CLEAN / "B_cells.h5ad",
    "DC_cells":          CLEAN / "DC_cells.h5ad",
    "Endothelial_cells": CLEAN / "Endothelial.h5ad",
    "Epithelial":        CLEAN / "Epithelial.h5ad",
    "Fibroblast":        CLEAN / "Fibroblast.h5ad",
    "Mast_cells":        ROUND_4_MAJOR / "mast_cell_integrated.h5ad",
    "MoMac":             CLEAN / "MoMac.h5ad",
    "Neutrophils":       REBUILT / "Neutrophils_sound.h5ad",
    "NK_cells":          CLEAN / "NK_cells.h5ad",
    "Pericyte":          CLEAN / "Pericyte.h5ad",
    "Plasma_cells":      ROUND_4_MAJOR / "plasma_cell_integrated.h5ad",
    "TCD4_cells":        CLEAN / "TCD4.h5ad",
    "TCD8_cells":        CLEAN / "TCD8.h5ad",
}

PHASES = {
    "pre":  dict(column="stomach_pre_grouping",  ref="Responsed", test="No-response"),
    "post": dict(column="stomach_post_grouping", ref="Responsed", test="No-response"),
}

SEED = 42
MIN_PCT = 0.10
MIN_CELLS_PER_GROUP = 20
LADDER_CELLS = 40
LADDER_TOL = 0.01


class BadMatrix(RuntimeError):
    pass


def assert_log1p_cp10k(X, source):
    """Integer-ladder test, unchanged from recompute_deg.py."""
    Xc = X.tocsr() if sparse.issparse(X) else np.asarray(X)
    rng = np.random.default_rng(SEED)
    rows = rng.choice(Xc.shape[0], min(LADDER_CELLS, Xc.shape[0]), replace=False)
    worst, tested, integral = 0.0, 0, 0
    for i in rows:
        v = (Xc.data[Xc.indptr[i]:Xc.indptr[i + 1]] if sparse.issparse(Xc)
             else Xc[i][Xc[i] > 0]).astype(np.float64)
        if v.size < 20:
            continue
        tested += 1
        if np.abs(v - np.round(v)).max() < 1e-6:
            integral += 1
            continue
        c = np.expm1(v) / np.expm1(v.min())
        worst = max(worst, float(np.abs(c - np.round(c)).max()))
    if tested == 0:
        raise BadMatrix(f"{source}: too few detected genes to test the matrix")
    if integral > tested / 2:
        raise BadMatrix(f"{source}: .raw holds integer counts, not log1p CP10K.")
    if worst >= LADDER_TOL:
        raise BadMatrix(
            f"{source}: .raw is NOT log1p CP10K - implied counts miss the "
            f"integers by up to {worst:.3f}.")


def load_subset(cell, phase):
    """The cells of one population at one timepoint, from the named slot."""
    path, slot = CELL_SOURCES[cell], SLOT[cell]
    cfg = PHASES[phase]
    adata = sc.read_h5ad(path)
    if slot == "raw":
        if adata.raw is None:
            raise BadMatrix(f"{path.name}: declared .raw, but there is none")
    elif adata.raw is not None:
        raise BadMatrix(
            f"{path.name}: declared .X, but the file carries a .raw as well - "
            "which slot is the log1p matrix is no longer unambiguous. Stopping.")
    m = adata.obs[cfg["column"]].astype(str).isin([cfg["ref"], cfg["test"]])
    if m.sum() == 0:
        return None
    sub = adata[m]
    if slot == "raw":
        out = sub.raw.to_adata()[sub.obs_names]
        out.obs = sub.obs.copy()
    else:
        out = sub.copy()
        out.raw = None
    assert_log1p_cp10k(out.X, f"{path.name} .{slot}")
    return out


def filter_genes(X):
    return (np.asarray((X > 0).mean(axis=0)).ravel() >= MIN_PCT)
