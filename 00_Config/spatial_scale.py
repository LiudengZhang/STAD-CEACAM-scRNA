"""Convert Visium distances from image pixels to micrometres.

spot_data.csv gives spot x and y in full-resolution image pixels, so every
distance derived from them (distance_to_stroma, distance_to_immune) is in
pixels. The factor is derived from the array geometry rather than from a
per-slide scalefactors_json.json: Visium spots sit on a hexagonal lattice with
a 100 um centre-to-centre pitch, so

    um_per_unit = SPOT_PITCH_UM / median(nearest-neighbour distance)

cohort_um_per_unit() returns one factor for all sections (pooled median
spacing) and raises if the per-section factors differ by more than MAX_SPREAD.

Run `python spatial_scale.py` to print the per-section table.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import SPATIAL_SPOT_DATA  # noqa: E402

__all__ = ["SPOT_PITCH_UM", "MAX_SPREAD", "microns_per_pixel",
           "per_sample_um_per_unit", "cohort_um_per_unit"]

# Visium centre-to-centre spot pitch (a property of the capture area).
SPOT_PITCH_UM = 100.0

# Largest relative spread of per-section factors allowed for one cohort factor.
MAX_SPREAD = 0.05

_cache: dict[str, object] = {}


def microns_per_pixel(xy):
    """Derived from the array geometry rather than a per-slide scale factor."""
    xy = np.asarray(xy, float)
    if xy.ndim != 2 or xy.shape[1] != 2 or len(xy) < 2:
        raise ValueError(f"need an (n, 2) array of at least two spots, got "
                         f"{getattr(xy, 'shape', None)}")
    d, _ = cKDTree(xy).query(xy, k=2)
    nn = float(np.median(d[:, 1]))
    if not np.isfinite(nn) or nn <= 0:
        raise ValueError(
            "median nearest-neighbour spacing is not positive; these "
            "coordinates are not a Visium array")
    return SPOT_PITCH_UM / nn


def _spot_table(spot=None):
    if spot is not None:
        return spot
    if "spot" not in _cache:
        _cache["spot"] = pd.read_csv(SPATIAL_SPOT_DATA)
    return _cache["spot"]


def per_sample_um_per_unit(spot=None):
    """One factor per section, so the spread can be read rather than assumed."""
    df = _spot_table(spot)
    return {s: microns_per_pixel(g[["x", "y"]].values)
            for s, g in df.groupby("sample")}


def cohort_um_per_unit(spot=None):
    """
    The cohort factor, um per pixel.

    Uses the median nearest-neighbour spacing pooled over all spots of all
    sections. Raises if the per-section factors disagree by more than
    `MAX_SPREAD`.
    """
    df = _spot_table(spot)
    per = per_sample_um_per_unit(df)
    lo, hi = min(per.values()), max(per.values())
    spread = (hi - lo) / np.mean(list(per.values()))
    if spread > MAX_SPREAD:
        raise RuntimeError(
            f"the per-section micron factors span {spread:.1%} "
            f"({lo:.6f} to {hi:.6f} um/unit), above MAX_SPREAD "
            f"{MAX_SPREAD:.0%}. One cohort factor is a rescale only while the "
            f"sections agree; at this spread the caller has to be told, not "
            f"averaged for. Sections: "
            + ", ".join(f"{k} {v:.6f}" for k, v in sorted(per.items())))
    nn = []
    for _, g in df.groupby("sample"):
        xy = g[["x", "y"]].values
        d, _i = cKDTree(xy).query(xy, k=2)
        nn.append(d[:, 1])
    return SPOT_PITCH_UM / float(np.median(np.concatenate(nn)))


def main():
    per = per_sample_um_per_unit()
    factor = cohort_um_per_unit()
    df = _spot_table()
    print(f"spot pitch                 {SPOT_PITCH_UM:g} um "
          f"(Visium centre-to-centre)")
    print(f"sections                   {len(per)}   spots {len(df):,}")
    print("\n  section        spots   nn median (units)   um/unit")
    for s, g in df.groupby("sample"):
        xy = g[["x", "y"]].values
        d, _i = cKDTree(xy).query(xy, k=2)
        print(f"  {s:<12} {len(g):>7}   {np.median(d[:, 1]):>16.4f}   "
              f"{per[s]:.6f}")
    lo, hi = min(per.values()), max(per.values())
    print(f"\n  per-section spread       {(hi - lo) / np.mean(list(per.values())):.2%} "
          f"({lo:.6f} to {hi:.6f})")
    print(f"  COHORT FACTOR            {factor:.9f} um per array unit")
    print(f"  one array unit           {factor:.6f} um; "
          f"1000 units = {1000 * factor:.2f} um")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
