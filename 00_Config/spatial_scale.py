"""The one place a Visium distance in array units becomes micrometres.

WHY THIS MODULE EXISTS
----------------------
`spot_data.csv` carries `x` and `y` in the full-resolution image's own pixels,
and every distance derived from them - `distance_to_stroma`,
`distance_to_immune`, and the panels and tests built on those - is in those
pixels. Nothing in that path converts. Until 2026-09-10 the Results quoted the
Fig. 3I and 3J mean differences as `371 um` and `179 um` when they were 371 and
179 array units, while both printed panels' y axes read `Distance (a.u.)` and
the published page's text layer contained no `um` anywhere. The figure agreed
with the data; the sentence did not agree with either.

The conversion itself was already in the tree, in
`04_Revision_Analyses/05_R1.6_Spatial_Confounders/scripts/spatial_distance_map.py`,
where it put Figure R2's axes into microns for the response letter and went no
further. Under RULES.md rule 5 a quantity that appears in a figure, the body
and a supplementary table has one owner, so it is here now and that script
imports it. Do not write a second derivation.

HOW THE FACTOR IS DERIVED
-------------------------
Not from a per-slide `scalefactors_json.json`, which is not in this repository,
but from the array geometry itself. A Visium capture area is a fixed hexagonal
lattice with a 100 um centre-to-centre pitch, so the median nearest-neighbour
spacing of a section's spots *is* 100 um, whatever the scan magnification:

    um_per_unit = SPOT_PITCH_UM / median(nearest-neighbour distance)

ONE FACTOR FOR THE COHORT, NOT TEN
----------------------------------
Measured on the ten GSE251950 sections the nearest-neighbour spacing is 347.0
to 354.0 units and the derived factor 0.282476 to 0.288180 um/unit - a spread
of 1.99 per cent, carried almost entirely by one section (sample_04 at 354.0
where the other nine sit at 347.0-348.0).

`cohort_um_per_unit()` returns one factor for the whole cohort, and that is a
decision, not an average taken quietly:

  * The published quantities are cohort-level - a paired difference of
    sample-level means over ten sections, its bootstrap confidence interval,
    and two boxplots of those same ten pairs. Multiplying every section by one
    constant is a pure linear rescale, so every rank test returns exactly the P
    value it returned before and every percentile bootstrap interval is the old
    interval times the factor. Nothing but the unit moves, which is the whole
    of what was authorised.
  * Converting each section by its own factor is defensible too and gives
    106.93 and 51.79 um against 106.54 and 51.47 um here - a 0.4 and 0.6 per
    cent difference. It is not used, because it would re-weight the ten pairs
    against each other, which is an edit to the statistic rather than to its
    unit, and because the interval would then have to be re-bootstrapped from a
    shared RNG stream that cannot be re-entered at that row.
  * `per_sample_um_per_unit()` is here so the spread can be looked at rather
    than assumed, and `cohort_um_per_unit()` REFUSES to return a single factor
    if that spread ever exceeds `MAX_SPREAD` - at which point one factor stops
    being a rescale and the caller has to be told, not quietly averaged for.

Run `python spatial_scale.py` for the derivation, the per-section table, and a
self-test that is shown failing before it is believed.
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
           "per_sample_um_per_unit", "cohort_um_per_unit", "self_test"]

#: Visium centre-to-centre spot pitch. A property of the capture area, not of
#: any one slide or scan.
SPOT_PITCH_UM = 100.0

#: The most the per-section factors may disagree before one cohort factor stops
#: being a rescale. 5 per cent, against a measured 1.99.
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
    The one factor the figure, the sweep and the Results all use.

    Pooled: the median nearest-neighbour spacing over every spot of every
    section, so a large section does not get one section's worth of vote and a
    small one the same. Raises rather than averaging if the per-section factors
    disagree by more than `MAX_SPREAD`.
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


# --------------------------------------------------------------- self-test
def _lattice(pitch, n=12, origin=(0.0, 0.0)):
    """A hexagonal patch with a known centre-to-centre spacing."""
    pts = []
    for r in range(n):
        for c in range(n):
            pts.append((origin[0] + pitch * (c + 0.5 * (r % 2)),
                        origin[1] + pitch * r * np.sqrt(3) / 2))
    return np.asarray(pts, float)


def self_test():
    """
    Show the derivation right, then show it wrong. A factor that cannot come
    out wrong is not a measurement of anything.
    """
    fails = []

    # 1. A lattice at a known pitch returns the factor that pitch implies.
    for pitch, want in ((347.0, SPOT_PITCH_UM / 347.0),
                        (1.0, SPOT_PITCH_UM),
                        (354.0, SPOT_PITCH_UM / 354.0)):
        got = microns_per_pixel(_lattice(pitch))
        ok = abs(got - want) < 1e-9
        print(f"  lattice pitch {pitch:>6.1f} -> {got:.9f} um/unit "
              f"(want {want:.9f})  {'ok' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"lattice {pitch}")

    # 2. MUTATION. Halve the lattice pitch and the factor must double. If the
    #    check accepts the old factor for the new geometry it is reading
    #    nothing off the array and every number downstream is a typed constant.
    mutant = microns_per_pixel(_lattice(347.0 / 2))
    caught = abs(mutant - SPOT_PITCH_UM / 347.0) > 1e-9
    print(f"  MUTANT half-pitch lattice -> {mutant:.9f} um/unit; "
          f"differs from the unmutated factor: {'yes' if caught else 'NO'}")
    if not caught:
        fails.append("mutation not caught: the factor ignores the geometry")

    # 3. MUTATION. Coordinates that are not an array at all must raise, not
    #    return a number that would silently rescale a figure.
    for bad, why in ((np.zeros((5, 2)), "every spot at the same place"),
                     (np.zeros((1, 2)), "one spot")):
        try:
            microns_per_pixel(bad)
            print(f"  MUTANT {why}: returned a factor  FAIL")
            fails.append(f"mutation not caught: {why}")
        except ValueError as exc:
            print(f"  MUTANT {why}: refused ({str(exc)[:48]}...)  ok")

    # 4. MUTATION. A cohort whose sections disagree beyond MAX_SPREAD must be
    #    refused rather than averaged.
    wide = pd.concat([
        pd.DataFrame(_lattice(347.0), columns=["x", "y"]).assign(sample="a"),
        pd.DataFrame(_lattice(500.0), columns=["x", "y"]).assign(sample="b"),
    ])
    try:
        cohort_um_per_unit(wide)
        print("  MUTANT sections 44% apart: returned one factor  FAIL")
        fails.append("mutation not caught: wide spread averaged")
    except RuntimeError as exc:
        print(f"  MUTANT sections 44% apart: refused "
              f"({str(exc)[:48]}...)  ok")

    return 1 if fails else 0


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
    print("\nself-test")
    return self_test()


if __name__ == "__main__":
    raise SystemExit(main())
