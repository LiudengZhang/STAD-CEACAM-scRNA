"""
Internal specimen identifier -> the de-identified study sample ID.

Every analysis output that names a specimen must name it the way Supplementary
Table 1 does (P01-M1, P32-P1, ...). The crosswalk is deliberately not stored in
this repository and is not hard-coded anywhere: the dataset objects carry both
labels in .obs - `sample` holds the identifier the sequencing core assigned and
`Sample ID` holds the study ID - so the mapping is resolved at run time from the
object the analysis is already reading.

Nothing here falls back. A specimen that does not map, a crosswalk that is not
one-to-one, or a study ID that does not look like a study ID raises. The
messages report counts and formats only, never an identifier, because keeping
the internal labels out of files and logs is the whole point.
"""

import re

import numpy as np
import pandas as pd

INTERNAL_COL = "sample"
STUDY_COL = "Sample ID"

# P01-M1, P32-P1, ... as printed in Supplementary Table 1.
STUDY_ID_RE = re.compile(r"^P\d+-\w+$")


def sample_id_map(obs):
    """
    {internal identifier: study sample ID}, read from an AnnData .obs.

    Pass the obs of the object the analysis actually reads, and pass it before
    any subsetting, so the crosswalk is checked against every specimen in it
    rather than against the handful a particular contrast happens to keep.
    """
    missing = [c for c in (INTERNAL_COL, STUDY_COL) if c not in obs.columns]
    if missing:
        raise ValueError(
            f"the specimen crosswalk needs .obs columns {missing!r}; without "
            "them the study sample IDs cannot be resolved")

    pairs = obs[[INTERNAL_COL, STUDY_COL]].astype(str).drop_duplicates()

    n_fanout = int((pairs.groupby(INTERNAL_COL)[STUDY_COL].nunique() > 1).sum())
    if n_fanout:
        raise ValueError(
            f"{n_fanout} internal identifier(s) carry more than one "
            f"'{STUDY_COL}'; the crosswalk is not one-to-one")

    n_collision = int(
        (pairs.groupby(STUDY_COL)[INTERNAL_COL].nunique() > 1).sum())
    if n_collision:
        raise ValueError(
            f"{n_collision} study sample ID(s) are shared by more than one "
            "internal identifier; the crosswalk is not one-to-one")

    n_bad = sum(1 for v in pairs[STUDY_COL] if not STUDY_ID_RE.match(v))
    if n_bad:
        raise ValueError(
            f"{n_bad} of {len(pairs)} '{STUDY_COL}' values do not match "
            f"{STUDY_ID_RE.pattern}, so they are not publishable study sample "
            "IDs")

    return dict(zip(pairs[INTERNAL_COL], pairs[STUDY_COL]))


def to_study_ids(values, mapping):
    """
    Map internal identifiers to study sample IDs, as a plain ndarray.

    An ndarray rather than a Series so that assigning the result back onto a
    column cannot silently re-align on an index and blank a row. Anything
    unmapped raises: an output that quietly kept a raw identifier is the
    failure this module exists to prevent.
    """
    s = pd.Series(list(values), dtype=object).astype(str)
    out = s.map(mapping)
    n_missing = int(out.isna().sum())
    if n_missing:
        raise ValueError(
            f"{n_missing} of {len(s)} specimen labels have no entry in the "
            "crosswalk, so they cannot be de-identified")
    return np.asarray(out, dtype=object)
