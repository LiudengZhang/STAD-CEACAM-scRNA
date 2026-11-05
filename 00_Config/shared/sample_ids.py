"""
Map internal specimen identifiers to de-identified study sample IDs
(P01-M1, P32-P1, ... as in Supplementary Table 1).

The crosswalk is not stored in the repository: the dataset objects carry both
labels in .obs (`sample`, the internal identifier, and `Sample ID`, the study
ID), so the mapping is read from the object being analysed. An unmapped
specimen, a crosswalk that is not one-to-one, or a malformed study ID raises.
Error messages report counts only, never identifiers.
"""

import re

import numpy as np
import pandas as pd

INTERNAL_COL = "sample"
STUDY_COL = "Sample ID"

STUDY_ID_RE = re.compile(r"^P\d+-\w+$")


def sample_id_map(obs):
    """
    {internal identifier: study sample ID}, read from an AnnData .obs.

    Pass the obs of the object being analysed, before any subsetting, so every
    specimen is checked.
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

    An ndarray, so assigning it to a column cannot re-align on an index.
    Anything unmapped raises.
    """
    s = pd.Series(list(values), dtype=object).astype(str)
    out = s.map(mapping)
    n_missing = int(out.isna().sum())
    if n_missing:
        raise ValueError(
            f"{n_missing} of {len(s)} specimen labels have no entry in the "
            "crosswalk, so they cannot be de-identified")
    return np.asarray(out, dtype=object)
