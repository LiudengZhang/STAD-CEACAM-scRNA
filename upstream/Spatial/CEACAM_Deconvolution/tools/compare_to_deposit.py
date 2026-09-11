#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""Compare a reproduced GraphST h5ad against the DEPOSITED output.

The deposit may be given as the deposited .h5ad or as a proportions .csv.
MEASURED 2026-09-02: spot_data.csv -- the file the Figure 3 panels read --
agrees with the deposited .h5ad to 1.0e-16, and disagrees with the
*_proportions_14types.csv files by up to 0.497.  Those CSVs and the
spatial_domains/*.csv files are stale intermediates from an earlier run
(written 2025-12-19 13:29, about an hour BEFORE the batch that produced the
h5ads started at 14:27:57).  So the .h5ad is the correct comparison target;
pass a CSV only to re-measure that discrepancy.

Read-only on the deposit.  Reports the magnitude of any movement in the 14
cell-type proportions and in the spatial-domain assignment; adopts nothing.

usage: compare_to_deposit.py <repro.h5ad> <deposited.h5ad|proportions.csv> <sample_key>
"""
import sys
import numpy as np, pandas as pd, anndata as ad

REPRO, DEPOSIT, SAMPLE = sys.argv[1], sys.argv[2], sys.argv[3]

CELL_TYPES_14 = ["B cells","CD4+ T cells","CD8+ T cells","NK cells","DC cells",
                 "Mast cells","Neutrophils","Plasma cells","Monocytes/Macrophages",
                 "Endothelial cells","Fibroblast","Pericyte",
                 "Epi CEACAM-high","Epi CEACAM-low"]

a = ad.read_h5ad(REPRO, backed="r")
new = pd.DataFrame({c: np.asarray(a.obs[c]) for c in CELL_TYPES_14}, index=a.obs_names)
if DEPOSIT.endswith(".h5ad"):
    b = ad.read_h5ad(DEPOSIT, backed="r")
    old = pd.DataFrame({c: np.asarray(b.obs[c]) for c in CELL_TYPES_14 if c in b.obs.columns},
                       index=b.obs_names)
    old_dom = b.obs["spatial_domain"].astype(str) if "spatial_domain" in b.obs.columns else None
else:
    old = pd.read_csv(DEPOSIT, index_col=0)
    old_dom = None

print(f"sample={SAMPLE}")
print(f"repro spots={new.shape[0]}  deposit spots={old.shape[0]}")
common = new.index.intersection(old.index)
print(f"common barcodes={len(common)}  repro-only={len(new.index.difference(old.index))}  deposit-only={len(old.index.difference(new.index))}")
if len(common) == 0:
    raise SystemExit("REFUSING: no common barcodes")

cols = [c for c in CELL_TYPES_14 if c in old.columns]
missing = [c for c in CELL_TYPES_14 if c not in old.columns]
if missing:
    print(f"columns absent from deposit: {missing}")

n = new.loc[common, cols].to_numpy(float)
o = old.loc[common, cols].to_numpy(float)
d = np.abs(n - o)
print(f"\nPROPORTION_MAX_ABS_DIFF {np.nanmax(d):.6e}")
print(f"PROPORTION_MEAN_ABS_DIFF {np.nanmean(d):.6e}")
print(f"PROPORTION_RMS_DIFF {np.sqrt(np.nanmean(d**2)):.6e}")
print("\nper cell type (max abs, mean abs, deposit mean):")
for j, c in enumerate(cols):
    print(f"  {c:<24} {np.nanmax(d[:,j]):.6e}  {np.nanmean(d[:,j]):.6e}  {np.nanmean(o[:,j]):.6e}")

# CEACAM ratio is what Figure 3F/3G/3K/3L/3M read
if "Epi CEACAM-high" in cols and "Epi CEACAM-low" in cols:
    def ratio(m):
        hi = m[:, cols.index("Epi CEACAM-high")]; lo = m[:, cols.index("Epi CEACAM-low")]
        tot = hi + lo
        return np.where(tot > 0, hi / np.where(tot > 0, tot, 1), np.nan)
    rn, ro = ratio(n), ratio(o)
    ok = ~(np.isnan(rn) | np.isnan(ro))
    print(f"\nCEACAM_RATIO_MAX_ABS_DIFF {np.nanmax(np.abs(rn[ok]-ro[ok])):.6e}")
    print(f"CEACAM_RATIO_MEAN_ABS_DIFF {np.nanmean(np.abs(rn[ok]-ro[ok])):.6e}")
    med_n, med_o = np.nanmedian(rn), np.nanmedian(ro)
    print(f"CEACAM_RATIO_MEDIAN repro={med_n:.6f} deposit={med_o:.6f} delta={med_n-med_o:+.6e}")
    grp_n = rn > med_n; grp_o = ro > med_o
    print(f"CEACAM_HIGH_LOW_GROUP_FLIPS {int((grp_n[ok]!=grp_o[ok]).sum())} of {int(ok.sum())}")

# spatial domain agreement (a relabelling is not necessarily a disagreement, so
# the confusion structure is reported alongside the raw agreement rate)
if old_dom is not None and "spatial_domain" in a.obs.columns:
    nd = a.obs["spatial_domain"].astype(str).loc[common]
    od = old_dom.loc[common]
    print(f"\nSPATIAL_DOMAIN_N_REPRO {nd.nunique()}  N_DEPOSIT {od.nunique()}")
    print(f"SPATIAL_DOMAIN_IDENTICAL_LABEL_FRACTION {(nd.values == od.values).mean():.6f}")
    ct = pd.crosstab(od, nd)
    best = ct.max(axis=1).sum() / len(common)
    print(f"SPATIAL_DOMAIN_BEST_RELABEL_AGREEMENT {best:.6f}")
