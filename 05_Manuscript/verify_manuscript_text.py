"""
Check that the sentences in the revised manuscript still agree with the tables
they cite.

verify_numbers.py checks the analysis outputs, and it cannot see the
manuscript. That is how one-tailed P values survived in the Results after every
figure and Table S6 had already been converted: nothing compared the prose with
the table. This script closes that gap.

Each claim below pins the phrase the manuscript uses, with the number left as a
slot that is filled from the table. The phrase then has to occur verbatim in the
clean docx. If the table moves, the filled phrase stops matching; if someone
edits the sentence, the phrase stops matching. Either way this fails and names
the sentence.

Not shipped in the code release: it needs the manuscript, which the release does
not carry.

Run: python verify_manuscript_text.py
Exit status is non-zero if any check fails.
"""

from pathlib import Path
import argparse
import re
import sys

import docx
import pandas as pd

HERE = Path(__file__).resolve().parent

# Defaults are the working-tree layout. The submission archive puts the same two
# inputs elsewhere, so both are settable; run it there with
#   python verify_manuscript_text.py \
#       --docx ../01_Manuscript/Manuscript_CIR-26-0753-ET_revised.docx \
#       --tables ../04_Supplementary_Tables
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--docx", type=Path,
                default=HERE / "01_Main_Text" / "Manuscript_R1_clean.docx",
                help="the clean revised manuscript")
ap.add_argument("--tables", type=Path, default=HERE / "04_Tables",
                help="directory holding ST6, ST8 and ST9")
args = ap.parse_args()

CLEAN = args.docx
ST6 = args.tables / "ST6_two_sided_sensitivity.csv"
ST8 = args.tables / "ST8_response_label_sensitivity.csv"
ST9 = args.tables / "ST9_crosscohort_convergence.csv"
for f in (CLEAN, ST6, ST8, ST9):
    if not f.exists():
        sys.exit(f"not found: {f}\nPass --docx and --tables for this layout.")

failures, checks = [], 0
TEXT = "\n".join(p.text for p in docx.Document(CLEAN).paragraphs)

sweep = pd.read_csv(ST6)
sens = pd.read_csv(ST8)
conv = pd.read_csv(ST9)


def two_tailed(fragment):
    rows = sweep[sweep["Comparison"].str.contains(fragment, regex=False)]
    if len(rows) != 1:
        return None
    return float(rows["P, two-tailed (revised)"].iloc[0])


def scenario_p(comparison_fragment, scenario_fragment):
    rows = sens[sens["Comparison"].str.contains(comparison_fragment, regex=False)
                & sens["Scenario"].str.contains(scenario_fragment, regex=False)]
    if len(rows) != 1:
        return None
    return float(rows["P, two-sided"].iloc[0])


def combined_p(gene, method="Stouffer, unweighted"):
    rows = conv[(conv["Analysis"] == "Combined across cohorts sharing no patients")
                & (conv["Measurement"] == gene) & (conv["Statistic"] == method)]
    return None if len(rows) != 1 else float(rows["P, two-sided"].iloc[0])


def concordance_rho(marker):
    rows = conv[(conv["Analysis"] == "Transcript fraction against protein staining")
                & (conv["Measurement"] == marker)]
    return None if len(rows) != 1 else float(rows["Value"].iloc[0])


def loo_range(gene, column):
    """Range of `column` across the eight leave-one-out refits, published excluded."""
    rows = conv[(conv["Analysis"] == "Leave-one-patient-out")
                & (conv["Measurement"] == gene)
                & (conv["Comparison"] != "none (as published)")]
    return None if len(rows) != 8 else (float(rows[column].min()),
                                        float(rows[column].max()))


def claim(label, template, *values):
    """The template, filled from the table, must appear verbatim in the text."""
    global checks
    checks += 1
    if any(v is None for v in values):
        failures.append(f"{label}: the table row it cites was not found")
        return
    phrase = template.format(*values)
    if phrase not in TEXT:
        failures.append(f"{label}: the manuscript does not contain "
                        f"“{phrase}”")


# --------------------------------------- Results against Table S6 (R1.3c)
claim("CEACAM6 pre-treatment epithelium",
      "two-sided P = {:.3f}, rank-biserial r = 0.88",
      two_tailed("CEACAM6 expression, pre-treatment"))
claim("CEACAM5 pre-treatment epithelium",
      "two-sided P = {:.3f}, r = 0.75",
      two_tailed("CEACAM5 expression, pre-treatment"))
claim("CEACAM5 in PRJEB25780",
      "CEACAM5, two-sided P = {:.3f}, r = 0.37",
      two_tailed("CEACAM5 expression, PRJEB25780"))
claim("CEACAM6 in PRJEB25780",
      "CEACAM6, two-sided P = {:.3f}, r = 0.38",
      two_tailed("CEACAM6 expression, PRJEB25780"))
claim("IHC summed score",
      "two-sided P = {:.3f}, r = 0.88, Mann–Whitney U test, n = 4 per group",
      two_tailed("CEACAM5 + CEACAM6 summed"))
claim("IHC markers separately",
      "(CEACAM5, P = {:.2f}; CEACAM6, P = {:.2f}; Table S7)",
      two_tailed("IHC staining, CEACAM5 only"),
      two_tailed("IHC staining, CEACAM6 only"))
claim("S-MP4 metaprogram",
      "two-sided exact permutation P = {:.3f}",
      two_tailed("S-MP4 score"))
claim("S-MP5 metaprogram",
      "(P = {:.3f}; Fig. 2H, I)",
      two_tailed("S-MP5 score"))
claim("neighborhood epithelial density",
      "two-sided P = {:.3f}, matched-pairs rank-biserial r = 0.85",
      two_tailed("Neighbourhood epithelial density"))
claim("distance to stroma",
      "two-sided P = {:.3f}, Wilcoxon signed-rank test; Fig. 3I",
      two_tailed("Distance to stroma"))
claim("distance to immune-rich regions",
      "P = {:.3f}; Fig. 3J",
      two_tailed("Distance to immune-rich"))
claim("BACH1 regulon",
      "two-sided exact permutation test, P = {:.3f}; Table S6)",
      two_tailed("BACH1 regulon activity"))
claim("NFKB1 regulon",
      "two-sided exact permutation test, P = {:.3f}; Fig. 5E, Table S6)",
      two_tailed("NFKB1 regulon activity"))
claim("PD-L1 in dendritic cells",
      "in dendritic cells (P = {:.3f})",
      two_tailed("post-treatment dendritic cells"))
claim("PD-L1 in the other three cell types",
      "(P = {:.3f}, {:.3f} and {:.3f}, respectively)",
      two_tailed("post-treatment monocytes/macrophages"),
      two_tailed("post-treatment epithelial cells"),
      two_tailed("post-treatment fibroblasts"))

# ------------------------------- Results and Limitations against Table S8
claim("response-label sensitivity, Results",
      "(two-sided P = {:.2f} when P26 is treated as a responder, versus "
      "P = {:.3f} as classified)",
      scenario_p("double-positive", "P26 reclassified as responder"),
      scenario_p("double-positive", "as published"))
claim("response-label sensitivity, Limitations",
      "from P = {:.3f} to P = {:.2f} (Table S8)",
      scenario_p("double-positive", "as published"),
      scenario_p("double-positive", "P26 reclassified as responder"))

# ------------------------------------- Results and Limitations against Table S9
# The sentence now names the weighting and reports the two other combinations,
# so that quoting the unweighted value cannot read as selection. All four values
# come from Table S9, so all four are checked here rather than only the pair.
claim("cross-cohort combination, Results",
      "combined two-sided P = {:.3f} for CEACAM6 and P = {:.3f} for CEACAM5 "
      "(unweighted Stouffer's method; the sqrt(n)-weighted and Fisher "
      "combinations agree, at P = {:.3f} and {:.3f} for CEACAM6 and "
      "P = {:.3f} and {:.3f} for CEACAM5; Table S9)",
      combined_p("CEACAM6"), combined_p("CEACAM5"),
      combined_p("CEACAM6", "Stouffer, weighted by sqrt(n)"),
      combined_p("CEACAM6", "Fisher"),
      combined_p("CEACAM5", "Stouffer, weighted by sqrt(n)"),
      combined_p("CEACAM5", "Fisher"))
claim("transcript-protein concordance, Results",
      "(summed, Spearman \u03c1 = {:.2f}, P = 0.004; CEACAM5, \u03c1 = {:.2f}, "
      "P = 0.007; CEACAM6, \u03c1 = {:.2f}, P = 0.10; Table S9)",
      concordance_rho("Summed"), concordance_rho("CEACAM5"),
      concordance_rho("CEACAM6"))
_lo_p, _hi_p = loo_range("CEACAM6", "P, two-sided") or (None, None)
_lo_g, _hi_g = loo_range("CEACAM6", "Value") or (None, None)
claim("leave-one-out stability, Limitations",
      "(CEACAM6, P between {:.3f} and {:.3f}, Hedges' g between {:.2f} and "
      "{:.2f}; Table S9)", _lo_p, _hi_p, _lo_g, _hi_g)

# The immunohistochemistry is not an independent cohort, and the Results must
# not say it is; ST1 and ST5 list the same eight patients.
checks += 1
if "immunohistochemical analysis on independent samples" in TEXT:
    failures.append("the Results still describe the immunohistochemistry as "
                    "performed on independent samples; Tables S1 and S5 list the "
                    "same eight patients")

# ------------------------------------------------ no one-tailed P survives
# The Methods statement and the Table S6 legend are the only two places the
# words belong; anywhere else means a sentence was missed in the conversion.
ALLOWED_ONE_TAILED = (
    "Statistical analyses were performed using Python",
    "Supplementary Table 6. Two-sided sensitivity analysis.",
)
checks += 1
for par in docx.Document(CLEAN).paragraphs:
    if not re.search(r"[Oo]ne-tail|[Oo]ne-sided", par.text):
        continue
    if not any(par.text.startswith(a) for a in ALLOWED_ONE_TAILED):
        failures.append("a one-tailed test is still declared outside Methods and "
                        f"the Table S6 legend: “{par.text[:90]}…”")

# -------------------------------- every supplementary item cited and legended
legends = {int(m.group(1)) for m in
           (re.match(r"Figure S(\d+)\.", p.text) for p in
            docx.Document(CLEAN).paragraphs) if m}
tab_legends = {int(m.group(1)) for m in
               (re.match(r"Supplementary Table (\d+)\.", p.text) for p in
                docx.Document(CLEAN).paragraphs) if m}
cited_figs = {int(m) for m in re.findall(r"Fig(?:s?\.|ure)? ?S(\d+)", TEXT)}
cited_tabs = {int(m) for m in re.findall(r"Table S(\d+)", TEXT)}

for label, have, want in (
        ("supplementary figure legends", legends, set(range(1, 10))),
        ("supplementary table legends", tab_legends, set(range(1, 11))),
        ("supplementary figures cited", cited_figs, set(range(1, 10))),
        ("supplementary tables cited", cited_tabs, set(range(1, 11)))):
    checks += 1
    missing = sorted(want - have)
    if missing:
        failures.append(f"{label}: missing {missing}")
    extra = sorted(have - want)
    if extra:
        failures.append(f"{label}: refers to items that do not exist {extra}")

# ------------------------- every panel cited by its own label, not "S9C, D"
# A compound citation reads correctly but is invisible to any tool, editorial or
# otherwise, that searches for the panel label. Each panel gets its full label.
PANELS = ("S2D S2E S2H S3E S7A S7B S7C S8A S8B "
          "S9A S9B S9C S9D S9E").split()
first_legend = next(i for i, par in enumerate(docx.Document(CLEAN).paragraphs)
                    if re.match(r"Figure S1\.", par.text))
BODY = "\n".join(par.text for par in
                 docx.Document(CLEAN).paragraphs[:first_legend])
checks += 1
uncited = [pn for pn in PANELS if not re.search(rf"\b{pn}\b", BODY)]
if uncited:
    failures.append("panels never cited by their own label in the body text: "
                    f"{uncited}")

print(f"Checked {checks} claims in the manuscript against the tables it cites.")
if failures:
    print(f"\n{len(failures)} FAILED:\n")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print("Every sentence agrees with its table.")
