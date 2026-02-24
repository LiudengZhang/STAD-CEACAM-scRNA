#!/usr/bin/env python3
"""
CEACAM5/6 IHC quantification using color deconvolution (Ruifrok & Johnston, 2001).

Upgrades over previous HSV-based pipeline:
  - Proper Hematoxylin-DAB color deconvolution via skimage
  - DAB optical density thresholding (not ad-hoc HSV ranges)
  - Tissue masking to exclude background

Analyzes 8 pre-treatment samples (4 R, 4 NR) used in manuscript Figure 2.
Outputs: CSV with per-sample staining area %, R vs NR comparison.
"""

import numpy as np
import pandas as pd
from PIL import Image
from skimage.color import rgb2hed, rgb2gray
from scipy import stats
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from paths import IHC_THUMBNAILS

# ── Paths ──
THUMB_DIR = IHC_THUMBNAILS
OUTPUT_DIR = Path(__file__).parent
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── 8 pre-treatment samples ──
# Patient → (CEACAM5 filename, CEACAM6 filename, Response group)
SAMPLES = {
    'P01': ('P01_CEACAM5.png',  'P01_CEACAM6.png',  'NR'),
    'P02': ('P02_CEACAM5.png',  'P02_CEACAM6.png',  'NR'),
    'P03': ('P03_CEACAM5.png',  'P03_CEACAM6.png',  'R'),
    'P04': ('P04_CEACAM5.png',  'P04_CEACAM6.png',  'R'),
    'P21': ('P21_CEACAM5.png',  'P21_CEACAM6.png',  'R'),
    'P22': ('P22_CEACAM5.png',  'P22_CEACAM6.png',  'R'),
    'P25': ('P25_CEACAM5.png',  'P25_CEACAM6.png',  'NR'),
    'P26': ('P26_CEACAM5.png',  'P26_CEACAM6.png',  'NR'),
}

# DAB optical density threshold for positive staining.
# Ruifrok deconvolution returns DAB channel where higher OD = more DAB.
# 0.02 selected after visual QC across 8 samples at 3 thresholds.
DAB_OD_THRESHOLD = 0.02


def tissue_mask(img_array, thresh=230):
    """Binary mask: tissue (True) vs white background (False)."""
    gray = rgb2gray(img_array)  # 0–1 range
    gray_uint8 = (gray * 255).astype(np.uint8)
    mask = gray_uint8 < thresh
    # Morphological cleanup
    from scipy.ndimage import binary_closing, binary_opening
    mask = binary_closing(mask, iterations=3)
    mask = binary_opening(mask, iterations=2)
    return mask


def quantify_dab(img_path, dab_threshold=DAB_OD_THRESHOLD):
    """
    Quantify DAB staining area using Hematoxylin-DAB color deconvolution.

    Returns:
        dict with tissue_area_px, dab_area_px, staining_pct
    """
    img = np.array(Image.open(img_path).convert('RGB'))

    # Color deconvolution: separate Hematoxylin, Eosin (unused), DAB
    # rgb2hed returns optical density in H, E, D channels
    hed = rgb2hed(img)
    dab_channel = hed[:, :, 2]  # DAB channel (index 2)

    # Tissue mask
    tmask = tissue_mask(img)
    tissue_px = int(np.sum(tmask))

    if tissue_px == 0:
        return {'tissue_area_px': 0, 'dab_area_px': 0, 'staining_pct': 0.0,
                'mean_dab_od': 0.0, 'mean_dab_od_positive': 0.0}

    # DAB positive = OD above threshold, within tissue
    dab_positive = (dab_channel > dab_threshold) & tmask
    dab_px = int(np.sum(dab_positive))

    # Mean DAB OD across tissue and across positive pixels
    dab_in_tissue = dab_channel[tmask]
    mean_od = float(np.mean(dab_in_tissue))
    mean_od_pos = float(np.mean(dab_channel[dab_positive])) if dab_px > 0 else 0.0

    return {
        'tissue_area_px': tissue_px,
        'dab_area_px': dab_px,
        'staining_pct': dab_px / tissue_px * 100,
        'mean_dab_od': mean_od,
        'mean_dab_od_positive': mean_od_pos,
    }


def main():
    print("=" * 60)
    print("CEACAM IHC Quantification — Color Deconvolution")
    print(f"DAB OD threshold: {DAB_OD_THRESHOLD}")
    print("=" * 60)

    rows = []
    for patient, (c5_file, c6_file, group) in SAMPLES.items():
        for marker, fname in [('CEACAM5', c5_file), ('CEACAM6', c6_file)]:
            img_path = THUMB_DIR / marker / fname
            if not img_path.exists():
                print(f"  WARNING: {img_path} not found, skipping")
                continue

            result = quantify_dab(img_path)
            row = {
                'patient': patient,
                'group': group,
                'marker': marker,
                'file': fname,
                **result,
            }
            rows.append(row)
            print(f"  {patient} {marker:8s} ({group}): "
                  f"{result['staining_pct']:.2f}% "
                  f"(DAB OD={result['mean_dab_od']:.3f})")

    df = pd.DataFrame(rows)
    csv_path = OUTPUT_DIR / 'ceacam_ihc_color_deconv_results.csv'
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved: {csv_path}")

    # ── Combined CEACAM5+6 per patient ──
    print("\n" + "=" * 60)
    print("COMBINED CEACAM5 + CEACAM6 STAINING (per patient)")
    print("=" * 60)

    pivot = df.pivot_table(index=['patient', 'group'],
                           columns='marker',
                           values='staining_pct').reset_index()
    pivot['combined'] = pivot['CEACAM5'] + pivot['CEACAM6']

    r_vals = pivot[pivot['group'] == 'R']['combined'].values
    nr_vals = pivot[pivot['group'] == 'NR']['combined'].values

    print(f"\nR  (n={len(r_vals)}): {np.round(r_vals, 2)}")
    print(f"   mean = {r_vals.mean():.2f}%")
    print(f"NR (n={len(nr_vals)}): {np.round(nr_vals, 2)}")
    print(f"   mean = {nr_vals.mean():.2f}%")
    print(f"Fold (NR/R): {nr_vals.mean() / r_vals.mean():.2f}x")

    stat, pval = stats.mannwhitneyu(nr_vals, r_vals, alternative='greater')
    print(f"\nMann-Whitney U (one-tailed, NR > R): U={stat:.1f}, P={pval:.4f}")

    # ── Per-marker comparison ──
    print("\n" + "=" * 60)
    print("PER-MARKER COMPARISON")
    print("=" * 60)
    for marker in ['CEACAM5', 'CEACAM6']:
        sub = df[df['marker'] == marker]
        r = sub[sub['group'] == 'R']['staining_pct'].values
        nr = sub[sub['group'] == 'NR']['staining_pct'].values
        u, p = stats.mannwhitneyu(nr, r, alternative='greater')
        print(f"\n{marker}:")
        print(f"  R:  {np.round(r, 2)} -> mean={r.mean():.2f}%")
        print(f"  NR: {np.round(nr, 2)} -> mean={nr.mean():.2f}%")
        print(f"  Mann-Whitney U={u:.1f}, P={p:.4f}")

    # ── Comparison with old HSV method ──
    print("\n" + "=" * 60)
    print("COMPARISON WITH OLD HSV-BASED VALUES")
    print("=" * 60)
    old_c5_r = np.array([3.1592, 3.6557, 0.0903, 2.8956])
    old_c5_nr = np.array([6.9229, 0.7542, 7.3245, 5.9232])
    old_c6_r = np.array([5.3573, 2.8627, 0.5448, 2.9016])
    old_c6_nr = np.array([4.2248, 13.2432, 7.1831, 3.6423])

    new_c5_r = pivot[pivot['group'] == 'R']['CEACAM5'].values
    new_c5_nr = pivot[pivot['group'] == 'NR']['CEACAM5'].values
    new_c6_r = pivot[pivot['group'] == 'R']['CEACAM6'].values
    new_c6_nr = pivot[pivot['group'] == 'NR']['CEACAM6'].values

    # Sort by patient for alignment
    r_patients = pivot[pivot['group'] == 'R']['patient'].values
    nr_patients = pivot[pivot['group'] == 'NR']['patient'].values

    print(f"\n{'Patient':>8s} {'Marker':>8s} {'Old HSV':>10s} {'New Deconv':>12s} {'Diff':>8s}")
    print("-" * 50)
    for i, p in enumerate(sorted(r_patients)):
        row_p = pivot[pivot['patient'] == p].iloc[0]
        # Find old value by matching patient order
        idx = list(sorted(r_patients)).index(p)
        print(f"{p:>8s} {'C5':>8s} {old_c5_r[idx]:>10.2f} {row_p['CEACAM5']:>12.2f} "
              f"{row_p['CEACAM5'] - old_c5_r[idx]:>+8.2f}")
        print(f"{'':>8s} {'C6':>8s} {old_c6_r[idx]:>10.2f} {row_p['CEACAM6']:>12.2f} "
              f"{row_p['CEACAM6'] - old_c6_r[idx]:>+8.2f}")

    # Save combined summary
    summary_path = OUTPUT_DIR / 'ceacam_ihc_combined_summary.csv'
    pivot.to_csv(summary_path, index=False)
    print(f"\nCombined summary saved: {summary_path}")


if __name__ == '__main__':
    main()
