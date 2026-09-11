"""
Reviewer 1 point R1.3b - how far do the pre-treatment results depend on the
response label of a single patient?

cohort_audit_report.txt flags two pre-treatment specimens whose recorded RECIST
trajectory does not follow the Methods definition of non-response (SD or PD):

    P26-P1   labelled NR, RECIST trajectory PR, best response PR
    P02-P1   labelled NR, RECIST "follow-up loss", no evaluable response

Both labels are the ones assigned by the treating team and recorded in Table S1,
and both were re-checked against the clinical records and retained. The pre-
treatment comparison rests on four samples per group, though, so a reviewer is
entitled to ask what happens if either call is wrong. This script answers that
directly by repeating every pre-treatment comparison with each specimen
reclassified or excluded.

Outputs (04_Revision_Analyses/01_R1.3_Cohort_Pairing/outputs/)
  response_label_sensitivity.csv       one row per comparison per scenario
  response_label_sensitivity_report.txt
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import NEW_ANALYSES  # noqa: E402

OUT = NEW_ANALYSES / "01_R1.3_Cohort_Pairing" / "outputs"
CEACAM = NEW_ANALYSES / "04_R1.5_CEACAM5_vs_CEACAM6" / "outputs"

# The two specimens the cohort audit could not reconcile with the SD/PD rule.
RECLASSIFY = "P26"
NOT_EVALUABLE = "P2"


def mann_whitney(nr, r):
    """Two-sided Mann-Whitney with rank-biserial r; NaN if either group is empty."""
    if len(nr) == 0 or len(r) == 0:
        return float("nan"), float("nan")
    u, p = mannwhitneyu(nr, r, alternative="two-sided")
    return p, 2 * u / (len(nr) * len(r)) - 1


def scenarios(df):
    """The four ways the two flagged specimens can be handled."""
    keep = df[df.patient != NOT_EVALUABLE]
    return [
        ("as published", df.assign(g=df.group)),
        (f"{RECLASSIFY} reclassified as responder",
         df.assign(g=np.where(df.patient == RECLASSIFY, "R", df.group))),
        (f"P02 excluded (no evaluable response)", keep.assign(g=keep.group)),
        (f"{RECLASSIFY} reclassified and P02 excluded",
         keep.assign(g=np.where(keep.patient == RECLASSIFY, "R", keep.group))),
    ]


def main():
    audit = pd.read_csv(OUT / "cohort_audit.csv")
    frac = pd.read_csv(CEACAM / "ceacam_state_fractions.csv")
    ihc = pd.read_csv(CEACAM / "ihc_per_marker_values.csv")

    # The scRNA table is keyed by the study sample ID; the audit carries that
    # column beside the study patient ID printed in the paper.
    crosswalk = dict(zip(audit["Sample"].astype(str), audit["Patient ID"]))
    frac["patient"] = frac["sample"].map(crosswalk)
    if frac["patient"].isna().any():
        raise SystemExit("unmapped scRNA samples: "
                         f"{frac.loc[frac.patient.isna(), 'sample'].tolist()}")
    # The IHC sheet zero-pads the patient IDs; the audit does not.
    ihc["patient"] = ihc["patient"].str.replace(r"^P0", "P", regex=True)

    comparisons = [
        (frac, "Double positive",
         "scRNA CEACAM5/6 double-positive epithelial fraction (Fig. 2E, Fig. S7B)"),
        (ihc, "Summed (published)", "IHC CEACAM5 + CEACAM6, summed (Fig. 2N)"),
        (ihc, "CEACAM5", "IHC CEACAM5 alone (Fig. S7C)"),
        (ihc, "CEACAM6", "IHC CEACAM6 alone (Fig. S7C)"),
    ]

    rows = []
    for table, column, label in comparisons:
        for name, d in scenarios(table):
            nr, r = d.loc[d.g == "NR", column], d.loc[d.g == "R", column]
            p, effect = mann_whitney(nr, r)
            rows.append(dict(
                Comparison=label, Scenario=name,
                **{"n (NR)": len(nr), "n (R)": len(r),
                   "Mean (NR)": round(nr.mean(), 4) if len(nr) else None,
                   "Mean (R)": round(r.mean(), 4) if len(r) else None,
                   "P, two-sided": round(p, 4),
                   "Rank-biserial r": round(effect, 3)}))

    result = pd.DataFrame(rows)
    result.to_csv(OUT / "response_label_sensitivity.csv", index=False)

    lines = [
        "SENSITIVITY OF THE PRE-TREATMENT RESULTS TO RESPONSE CLASSIFICATION",
        "Reviewer 1, point R1.3b", "=" * 96, "",
        "Two pre-treatment specimens do not follow the Methods definition of",
        "non-response (SD or PD):", "",
        "   P26-P1   labelled NR, RECIST trajectory PR, best response PR",
        "   P02-P1   labelled NR, RECIST follow-up loss, no evaluable response", "",
        "Both labels are as recorded by the treating team in Table S1 and are",
        "retained. This table shows what each pre-treatment comparison returns if",
        "either call is handled differently.", "",
    ]
    for label in result.Comparison.unique():
        lines += [label, "-" * 96]
        for _, row in result[result.Comparison == label].iterrows():
            lines.append(
                f"   {row.Scenario:<42} n = {row['n (NR)']} vs {row['n (R)']}   "
                f"NR {row['Mean (NR)']:<9} R {row['Mean (R)']:<9} "
                f"P = {row['P, two-sided']:<7} r = {row['Rank-biserial r']:+.3f}")
        lines.append("")
    lines += [
        "READING", "-" * 96,
        "At four samples per group a single reassignment moves the pre-treatment",
        "CEACAM result substantially. This is stated in the Results and the table",
        "is deposited as Table S8, so the dependence is visible rather than",
        "left for a reader to discover.", ""]
    (OUT / "response_label_sensitivity_report.txt").write_text("\n".join(lines),
                                                              encoding="utf-8")
    print("\n".join(lines))
    print(f"wrote {OUT}/response_label_sensitivity.csv and _report.txt")


if __name__ == "__main__":
    main()
