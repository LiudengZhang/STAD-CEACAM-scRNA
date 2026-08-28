"""
Acceptance test for patch_figure_annotations.py.

Reports, for each main figure, whether significance stars or exact P values
appear in the vector text. It is the check the pre-submission review shipped,
rewritten against PyMuPDF because that is what this environment has, and
extended to look at the drawing layer as well: Figure 2's annotations are vector
outlines with no text layer at all, so a text-only scan calls it clean whether
or not it has been patched.

Figure 4 carries label corrections rather than P values, so it is checked
differently: the strings that were replaced must be gone and the replacements
present.

Exit code 0 when every affected panel carries an exact P value and Figure 4
carries the corrected labels.

    python check_figure_annotations.py [--dir DIR]
"""

import argparse
import re
import sys
from pathlib import Path

import fitz

# The panels whose statistics were converted to two-sided, by printed letter.
AFFECTED = {
    "Figure_2": ["E cluster proportion", "H S-MP4", "I S-MP5", "K CEACAM5/6",
                 "L PRJEB25780", "N IHC"],
    "Figure_3": ["H epithelial density", "I distance to stroma",
                 "J distance to immune"],
    "Figure_5": ["D BACH1 regulon", "E NFKB1 regulon", "J PD-L1 four panels",
                 "L IL-6/JAK/STAT3"],
}
EXPECTED_P = {"Figure_2": 8, "Figure_3": 3, "Figure_5": 9}

# Stars that are not defects. Neither of these came from a one-tailed test, and
# neither is a two-group comparison where an exact value would fit: the first is
# a correlation between two markers at n = 31, the second is a nine-cell-type
# bar chart whose star key is given in the legend. They are counted so that a
# star appearing anywhere else still fails the check.
ALLOWED_STARS = {
    "Figure_2": (1, "panel D, Spearman rho between CEACAM5 and CEACAM6"),
    "Figure_3": (5, "panel N, immune recruitment across nine cell types"),
    "Figure_5": (0, ""),
}

# Figure 4: (string that must be absent, string that must be present, how many).
# "CD16" is deliberately not in the absent list - Mono_CD16 is a state name and
# keeps it; only the italic gene label was wrong.
LABEL_CHECKS = {
    "Figure_4": [("Mac_IL1B", "MoMac_IL1B", 6), (None, "FCGR3A", 1)],
}

STAR = re.compile(r"(?<![A-Za-z0-9])\*+(?![A-Za-z0-9])")
PVAL = re.compile(r"P\s*[=<]\s*0?\.\d+")


def scan(path):
    doc = fitz.open(path)
    text = "".join(page.get_text() for page in doc)
    doc.close()
    return len(STAR.findall(text)), len(PVAL.findall(text)), len(text.strip())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="Main_Figures/_patched")
    a = ap.parse_args()
    root = Path(a.dir)

    print(f"{'figure':<12}{'stars':>7}{'P values':>10}{'expected':>10}   verdict")
    print("-" * 66)
    bad = []
    for name in sorted(AFFECTED):
        path = root / f"{name}.pdf"
        if not path.exists():
            print(f"{name:<12}{'-':>7}{'-':>10}{'-':>10}   file not found")
            bad.append(name)
            continue
        stars, pvals, _ = scan(path)
        want = EXPECTED_P[name]
        allowed, reason = ALLOWED_STARS[name]
        if stars > allowed:
            verdict, ok = f"{stars - allowed} unexpected star(s)", False
        elif pvals < want:
            verdict, ok = f"only {pvals} of {want} exact P values", False
        else:
            verdict, ok = "exact P values", True
        print(f"{name:<12}{stars:>7}{pvals:>10}{want:>10}   {verdict}")
        if allowed and ok:
            print(f"{'':<12}{allowed:>7} allowed: {reason}")
        if not ok:
            bad.append(name)

    for name in sorted(LABEL_CHECKS):
        path = root / f"{name}.pdf"
        if not path.exists():
            print(f"{name:<12}{'-':>7}{'-':>10}{'-':>10}   file not found")
            bad.append(name)
            continue
        doc = fitz.open(path)
        text = "".join(page.get_text() for page in doc)
        doc.close()
        problems = []
        for old, new, count in LABEL_CHECKS[name]:
            got = len(re.findall(rf"(?<![A-Za-z0-9_]){re.escape(new)}", text))
            if got != count:
                problems.append(f"{new} appears {got} times, expected {count}")
            if old is not None and re.search(rf"(?<![A-Za-z0-9_]){re.escape(old)}",
                                             text):
                problems.append(f"{old} is still present")
        verdict = "corrected labels" if not problems else "; ".join(problems)
        print(f"{name:<12}{'n/a':>7}{'n/a':>10}{'n/a':>10}   {verdict}")
        if problems:
            bad.append(name)

    print()
    for name in bad:
        if name in AFFECTED:
            print(f"{name}: panels to check -> " + "; ".join(AFFECTED[name]))
    if bad:
        print("\nRegenerate with: python patch_figure_annotations.py")
        return 1
    print("Every converted panel prints its exact two-sided P value, and "
          "Figure 4 carries the corrected labels.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
