"""
WP1 / Reviewer 1 points R1.0, R1.3a, R1.3b.

Answers three questions the reviewer asked directly:
  1. How many samples are longitudinally paired (same patient, pre AND post)?
  2. Do the post-treatment non-responders represent primary non-responders,
     patients with an initial PR who then progressed, or both?
  3. What is the exact composition of the two response-labelled cohorts?

Sources
  ST1  : submission-tree/04_Manuscript/04_Tables/ST1_patient_sample_characteristics.csv
  CLIN : the hospital's paraffin-block worksheet, reached through the
         STAD_CLINICAL_XLSX environment variable. It carries patient
         identifiers and is not deposited.
         (holds the RECIST 1.1 trajectory column, e.g. "PR-SD", "SD-PD", and the
         treatment regimen with cycle counts; keyed by an internal scRNA ID that
         does not appear in ST1, so it is joined on the clinical fingerprint
         Age + Sex + Differentiation + sampling procedure + treatment phase)

Outputs (04_Revision_Analyses/01_R1.3_Cohort_Pairing/outputs/)
  cohort_audit.csv          one row per response-labelled stomach sample
  pairing_summary.csv       per-patient timepoint coverage
  cohort_audit_report.txt   plain-text answers for the response letter
"""

from pathlib import Path
import os
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import B_CELLS_H5AD, MANUSCRIPT  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

ST1 = MANUSCRIPT / "04_Tables" / "ST1_patient_sample_characteristics.csv"
H5AD_FOR_MAPPING = B_CELLS_H5AD

# The RECIST trajectory lives in the hospital's own spreadsheet, which carries
# patient identifiers and is therefore not deposited with the rest of the data.
# It has no default location: point STAD_CLINICAL_XLSX at it to re-run this
# audit. Everything downstream of it is summarised in the outputs below, which
# are released.
#
# 2026-09-11: it is NOT summarised in Supplementary Table 1 any more, and this
# comment said it was. ST1 used to carry three columns built from the audit's
# recist_raw / recist_best / recist_change; on the author's ruling it now
# reports the collaborators' DEFINITIVE adjudicated call in a single column,
# "RECIST 1.1 response", and does not report the original assessment at all.
# The two are not the same reading of the same patients - for the six
# post-treatment non-responders the trajectory's best responses are SD in four
# and PD in two, while the adjudicated calls are SD in two and PD in four - so
# the outputs here and ST1 must not be quoted against each other. This script
# is otherwise unaffected: it reads no RECIST column out of ST1, only Patient
# ID, Sample, Age, Sex, cTNM stage, Differentiation, Stomach site, Anatomical
# site, Biopsy method and R/NR Grouping.
CLIN = Path(os.environ["STAD_CLINICAL_XLSX"]) if "STAD_CLINICAL_XLSX" in os.environ \
    else None

PHASE_FROM_TIMING = {
    "Before Neoadjuvant treatment": "Pre",
    "After Neoadjuvant treatment": "Post",
}


def norm_procedure(s):
    """The two tables spell the sampling procedure slightly differently."""
    s = str(s).strip().lower()
    if "surgery" in s or "surgical" in s:
        return "Surgery"
    if "gastroscop" in s:
        return "Gastroscopic biopsy"
    return s


def classify_trajectory(recist):
    """
    Split the RECIST 1.1 field into best response and whether the disease
    subsequently worsened. "PR-SD" means an initial partial response that later
    became stable disease; "SD-PD" means stable disease that later progressed.
    """
    r = str(recist).strip()
    if r.lower().startswith("follow-up"):
        return pd.Series(["Unknown", "Unknown", "Follow-up lost"])
    parts = [p.strip() for p in r.split("-") if p.strip()]
    best = parts[0]
    last = parts[-1]
    order = {"CR": 0, "PR": 1, "SD": 2, "PD": 3}
    if len(parts) == 1:
        change = "No documented change"
    elif order.get(last, 9) > order.get(best, 9):
        change = "Worsened after best response"
    else:
        change = "Improved after first assessment"
    return pd.Series([best, last, change])


# ---------------------------------------------------------------- load ST1
st1 = pd.read_csv(ST1)
st1.columns = [c.strip() for c in st1.columns]

# ------------------------------------------------- Q1: longitudinal pairing
phases = (
    st1.groupby("Patient ID")["Treatment phase"]
    .apply(lambda s: sorted(set(s.dropna().astype(str))))
    .rename("phases_any_site")
)
stomach = st1[st1["Anatomical site"] == "Stomach"]
phases_stomach = (
    stomach.groupby("Patient ID")["Treatment phase"]
    .apply(lambda s: sorted(set(s.dropna().astype(str))))
    .rename("phases_stomach")
)
pairing = pd.concat([phases, phases_stomach], axis=1)
pairing["paired_any_site"] = pairing["phases_any_site"].apply(
    lambda v: isinstance(v, list) and {"Pre", "Post"} <= set(v)
)
pairing["paired_stomach"] = pairing["phases_stomach"].apply(
    lambda v: isinstance(v, list) and {"Pre", "Post"} <= set(v)
)
pairing = pairing.reset_index()
pairing.to_csv(OUT / "pairing_summary.csv", index=False)

n_paired_any = int(pairing["paired_any_site"].sum())
n_paired_stomach = int(pairing["paired_stomach"].sum())
paired_ids = pairing.loc[pairing["paired_any_site"], "Patient ID"].tolist()

# ------------------------------- the 19 response-labelled stomach specimens
labelled = stomach[stomach["R/NR Grouping"].notna()].copy()
labelled["procedure_key"] = labelled["Biopsy method"].map(norm_procedure)

# ------------------------------------------------- load the clinical table
if CLIN is None or not CLIN.exists():
    # A designed skip, not a failure. The table is deliberately not deposited,
    # so every reviewer running the capsule reaches this line; exiting 1 made
    # the driver record the deposit as having a broken script and made the
    # whole run non-zero for something that is working as intended. The message
    # is unchanged and still goes to stderr; only the status changes.
    print(
        "SKIPPED: the RECIST trajectory table is not available.\n"
        "It carries patient identifiers and is not part of the deposited data.\n"
        "Set STAD_CLINICAL_XLSX to a copy to re-run this audit; its results are\n"
        "reported in this directory's outputs/. Supplementary Table 1 reports the\n"
        "separately adjudicated definitive response, not this trajectory.",
        file=sys.stderr)
    sys.exit(0)
clin = pd.read_excel(CLIN)
clin.columns = [str(c).strip() for c in clin.columns]
clin = clin.rename(
    columns={
        "Therapeutic effect (Recist1.1)": "recist_raw",
        "Specimen collection timing": "timing",
        "Treatment regimes": "regimen",
        "Sampling procedure": "procedure",
        "cTNM stage": "cTNM",
        "yp/pTNM stage": "ypTNM",
        "scRNA ID": "internal_id",
    }
)
clin["Treatment phase"] = clin["timing"].map(PHASE_FROM_TIMING)
clin["procedure_key"] = clin["procedure"].map(norm_procedure)
clin["diff_key"] = clin["Differentiation"].astype(str).str.strip().str.lower()
clin["site_key"] = clin["Primary site"].astype(str).str.strip().str.lower()
clin["cTNM_key"] = clin["cTNM"].astype(str).str.strip().str.lower()

labelled["diff_key"] = labelled["Differentiation"].astype(str).str.strip().str.lower()
labelled["site_key"] = labelled["Stomach site"].astype(str).str.strip().str.lower()
labelled["cTNM_key"] = labelled["cTNM stage"].astype(str).str.strip().str.lower()

# The two tables share no key directly: ST1 uses P01-P1 style sample names while
# the clinical sheet uses the internal sequencing IDs. The h5ad obs carries both,
# so it is the authoritative crosswalk. B_cells.h5ad is used because it is the
# smallest object with the full obs schema, and it is opened backed so only the
# metadata is read.
def sample_id_crosswalk():
    import anndata

    a = anndata.read_h5ad(H5AD_FOR_MAPPING, backed="r")
    m = (a.obs[["Sample ID", "sample"]].drop_duplicates()
         .rename(columns={"Sample ID": "Sample", "sample": "internal_id"}))
    m["Sample"] = m["Sample"].astype(str)
    m["internal_id"] = m["internal_id"].astype(str)
    return m.reset_index(drop=True)


def norm_id(s):
    """CA-0923 in the h5ad, CA_0923 in the clinical sheet."""
    return str(s).strip().upper().replace("-", "_")


crosswalk = sample_id_crosswalk()
crosswalk["id_key"] = crosswalk["internal_id"].map(norm_id)
clin["id_key"] = clin["internal_id"].map(norm_id)

CLIN_COLS = ["internal_id", "recist_raw", "regimen", "ypTNM"]

merged = (labelled.reset_index(drop=True)
          .merge(crosswalk[["Sample", "id_key"]], on="Sample", how="left",
                 validate="one_to_one")
          .merge(clin[["id_key"] + CLIN_COLS], on="id_key", how="left",
                 validate="one_to_one"))

unmapped = merged[merged["internal_id"].isna()]["Sample"].tolist()
if unmapped:
    raise SystemExit(
        "These response-labelled samples could not be joined to the clinical "
        f"record via the h5ad crosswalk: {unmapped}")

merged[["recist_best", "recist_last", "recist_change"]] = merged["recist_raw"].apply(
    classify_trajectory
)
merged["mapping_tier"] = "h5ad sample-ID crosswalk"
merged["mapping_ambiguous"] = False
merged["paired_with_other_timepoint"] = merged["Patient ID"].isin(paired_ids)

cols = [
    "Patient ID", "Sample", "Treatment phase", "R/NR Grouping",
    "paired_with_other_timepoint", "Age", "Sex", "Differentiation",
    "Biopsy method", "regimen", "recist_raw", "recist_best", "recist_last",
    "recist_change", "ypTNM", "internal_id", "mapping_tier", "mapping_ambiguous",
]
audit = merged[cols].sort_values(["Treatment phase", "R/NR Grouping", "Patient ID"])
audit.to_csv(OUT / "cohort_audit.csv", index=False)

# ------------------------------------------------------------- the report
unmatched = audit[audit["recist_raw"].isna()]
post = audit[audit["Treatment phase"] == "Post"]
post_nr = post[post["R/NR Grouping"] == "NR"]
post_r = post[post["R/NR Grouping"] == "R"]

lines = []
w = lines.append
w("COHORT AUDIT - Reviewer 1 points R1.0, R1.3a, R1.3b")
w("=" * 78)
w("")
w("Q1. HOW MANY SAMPLES WERE LONGITUDINALLY PAIRED?")
w("-" * 78)
w(f"Total patients in the atlas                     : {st1['Patient ID'].nunique()}")
w(f"Total specimens                                 : {len(st1)}")
w(f"Patients with pre AND post at ANY site          : {n_paired_any}  {paired_ids}")
w(f"Patients with pre AND post in STOMACH           : {n_paired_stomach}")
if paired_ids:
    sub = st1[st1["Patient ID"].isin(paired_ids)]
    w("")
    w("  Paired specimens:")
    for _, r in sub.iterrows():
        w(f"    {r['Patient ID']:>4}  {r['Sample']:<8} {r['Anatomical site']:<10} "
          f"{r['Treatment phase']:<5} {str(r['R/NR Grouping'])}")
w("")
w("  => The response comparisons are between INDEPENDENT patient groups.")
w("     The design is timepoint-stratified and cross-sectional, not")
w("     longitudinally paired.")
w("")
w("")
w("Q2. WERE POST-TREATMENT NON-RESPONDERS PRIMARY NON-RESPONDERS, PATIENTS")
w("    WITH AN INITIAL PR WHO PROGRESSED, OR BOTH?")
w("-" * 78)
if len(unmatched):
    w(f"  WARNING: {len(unmatched)} sample(s) could not be matched to the clinical")
    w("  record and have no RECIST trajectory:")
    for _, r in unmatched.iterrows():
        w(f"    {r['Sample']} ({r['Treatment phase']}, {r['R/NR Grouping']})")
    w("")
w(f"  Post-treatment NON-RESPONDERS (n = {len(post_nr)}):")
for _, r in post_nr.iterrows():
    w(f"    {r['Sample']:<8} RECIST {str(r['recist_raw']):<16} best={r['recist_best']:<8} {r['recist_change']}")
w("")
w(f"  Post-treatment RESPONDERS (n = {len(post_r)}):")
for _, r in post_r.iterrows():
    w(f"    {r['Sample']:<8} RECIST {str(r['recist_raw']):<16} best={r['recist_best']:<8} {r['recist_change']}")
w("")
n_nr_ever_pr = int((post_nr["recist_best"].isin(["CR", "PR"])).sum())
n_r_worsened = int((post_r["recist_change"] == "Worsened after best response").sum())
w(f"  Post-treatment NR whose BEST response was CR/PR : {n_nr_ever_pr}")
w(f"  Post-treatment R who worsened after best response: {n_r_worsened}")
w("")
w("")
w("Q3. IS THE R/NR LABEL CONSISTENT WITH THE RECORDED RECIST RESPONSE?")
w("-" * 78)
w("RECIST 1.1 defines a responder as CR or PR. Any sample labelled R whose best")
w("response is SD/PD, or labelled NR whose best response is CR/PR, is flagged.")
w("")
audit["expected_group"] = audit["recist_best"].map(
    {"CR": "R", "PR": "R", "SD": "NR", "PD": "NR"}
)
audit["label_consistent"] = (
    audit["expected_group"].isna() | (audit["expected_group"] == audit["R/NR Grouping"])
)
bad = audit[~audit["label_consistent"]]
noresp = audit[audit["expected_group"].isna()]
w(f"  Samples whose label contradicts the recorded RECIST response : {len(bad)}")
for _, r in bad.iterrows():
    w(f"    {r['Sample']:<8} {r['Treatment phase']:<5} labelled {r['R/NR Grouping']:<3} "
      f"but RECIST = {r['recist_raw']} (best {r['recist_best']})")
w(f"  Samples with no evaluable RECIST response                    : {len(noresp)}")
for _, r in noresp.iterrows():
    w(f"    {r['Sample']:<8} {r['Treatment phase']:<5} labelled {r['R/NR Grouping']:<3} "
      f"but RECIST = {r['recist_raw']}")
if len(bad) or len(noresp):
    w("")
    w("  *** UNRESOLVED. These rows must be checked against the primary clinical")
    w("      records before the revision is submitted. Do not silently relabel. ***")
w("")
w("")
w("Q4. COHORT COMPOSITION")
w("-" * 78)
w(audit.groupby(["Treatment phase", "R/NR Grouping"]).size().to_string())
w("")
w("All treatment was given in the NEOADJUVANT setting; post-treatment specimens")
w("were taken at surgery or repeat gastroscopy after 1-6 cycles.")

report = "\n".join(lines)
(OUT / "cohort_audit_report.txt").write_text(report, encoding="utf-8")
print(report)
