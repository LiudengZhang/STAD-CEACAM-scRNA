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

The framing wording is checked across the whole shipped set and in both
directions: the qualifier withdrawn after a reviewer objected that the claim
was too broad must appear on no figure, and the narrower claim approved in its
place must appear on Figure 6 exactly once. Forbidding the old wording alone
would pass an edit that dropped the approved label as well.

One carried-over supplementary figure is checked the same way. Supplementary
Figure S2 panel D printed the one-sided exact Mann-Whitney value (P = 0.03);
patch_S2_two_sided.py replaces it with the two-sided value Table S6 reports
(P = 0.057). The shipped S2 must print the latter and not the former. The
supplementary figures sit beside the main ones in every layout this runs in -
03_Supplementary_Figures/ next to the archive's 02_Figures/, or
Supplementary_Fixes/_patched/ next to Main_Figures/_patched/ in the working
tree - so their directory is found from --dir, or given with --supp-dir. A
missing S2 is a failure, not a skip.

Exit code 0 when every affected panel carries an exact P value, Figure 4
carries the corrected labels, the framing wording is the approved one and S2
prints the two-sided value.

    python check_figure_annotations.py [--dir DIR] [--supp-dir DIR]
"""

import argparse
import hashlib
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

# The framing wording, asserted across the whole shipped set rather than inside
# one figure, and asserted in both directions.
#
# FRAMING_FORBIDDEN is the qualifier withdrawn after a reviewer objected that
# the claim was too broad. It must appear on no figure, so it is searched for
# on every one of them; a check that looked only where it used to be would pass
# a page that had reacquired it somewhere else. Matched case-insensitively,
# because the retired schematic label set it capitalised.
#
# FRAMING_REQUIRED is the narrower claim that replaced it, with the figure that
# carries it and how many times. Forbidding the old wording alone leaves an
# edit that quietly drops the approved label indistinguishable from a pass, so
# the replacement is required by name and by count.
FRAMING_FORBIDDEN = ("pan-TME",)
FRAMING_REQUIRED = {"Figure_6": [("Multi-Lineage", 1)]}

# Supplementary figures carried over from the submission and patched in place,
# same tuple form. Matched with a digit boundary on both sides, so that the
# replacement "P = 0.057" can never be read as the old "P = 0.03" still present.
SUPP_LABEL_CHECKS = {
    "S2_CEACAM_Metaprogram_Validation": [("P = 0.03", "P = 0.057", 1)],
}

STAR = re.compile(r"(?<![A-Za-z0-9])\*+(?![A-Za-z0-9])")
PVAL = re.compile(r"P\s*[=<]\s*0?\.\d+")

# Content loss is checked because a panel swap can silently delete a neighbour.
# The redaction removes any image whose bounding box touches its slot, and the
# Figure 2D slot overlaps panel F's Milo neighbourhood graph by 2.12 pt in the
# white gutter between them. A figure can therefore ship with a blank panel
# while every P value and label string this file checks is still intact.
#
# Images are matched by the md5 of their own bytes, not by where they sit on the
# page. That matters because Figure 1's panel A is replaced and the page grows,
# so every placement shifts; matching on position would report all eight of its
# images as lost. Matching on content reports zero, correctly, and needs no
# allow-list - which is the point, since an allow-list is another thing to get
# wrong.
#
# Only substantive images are compared. The threshold is measured, not chosen:
# across the six submitted figures the largest text-label raster is 542.5 pt^2
# and the smallest panel raster is 1026.5 pt^2, so anything in between separates
# them. Small rasters are legitimately removed when an annotation is redrawn as
# text, and counting those would make the check cry wolf on every run.
SUBSTANTIVE_PT2 = 800.0
# The submitted figures. They are the reference this check needs and they are
# NOT deposited: 00_GROUND_TRUTH/ exists only in the working tree, so in the
# release and in the submission archive this path resolves to nothing. If that
# were fatal, check_no_content_lost() would print "file not found" six times and
# return 1 for anyone running the check from the deposit rather than the working
# tree. --submitted overrides the location; when the directory is absent the
# content-loss check is skipped and says so, while the P-value and label checks
# - which need only --dir - stay fatal.
SUBMITTED = Path(__file__).parent.parent / "00_GROUND_TRUTH" / "figures"

# How each figure was made is read from SHIPPED_MANIFEST.csv. It sits beside the
# figures in the working tree, where --dir points into Main_Figures/, and beside
# this script in the deposit, where --dir points at a figures directory that has
# no manifest next to it. Both places are looked in, because the same file has
# to run in both layouts.
MANIFEST_NAME = "SHIPPED_MANIFEST.csv"

# Which mechanisms are compared against the submitted page, and what a missing
# image means for each. A figure that derives from that page carries the page's
# own rasters, so one that has gone missing is published content dropped:
# "patched" is the submitted page edited in place, and "supplied" is a corrected
# version of it handed over whole, which shares every raster with it and can
# drop one just as quietly. A figure redrawn at the size it prints at draws its
# own rasters and shares none, so the comparison has no meaning there.
COMPARED = {"patched": "deleted by patching",
            "supplied": "absent from the supplied page"}
NOT_COMPARED = {"slotted": "redrawn; its rasters are its own"}


def find_manifest(root):
    """The build manifest, from the figures directory or from beside this file."""
    for cand in (root.parent / MANIFEST_NAME,
                 Path(__file__).resolve().parent / MANIFEST_NAME):
        if cand.is_file():
            return cand
    return None


def substantive_images(path):
    """{md5 of image bytes} for every image drawn larger than the threshold."""
    doc = fitz.open(path)
    page = doc[0]
    out = set()
    for img in page.get_images(full=True):
        rects = page.get_image_rects(img[0])
        if not rects:
            continue
        if max(r.width * r.height for r in rects) < SUBSTANTIVE_PT2:
            continue
        out.add(hashlib.md5(doc.extract_image(img[0])["image"]).hexdigest())
    doc.close()
    return out


def check_no_content_lost(root, submitted):
    """Every substantive image in the submitted figure must survive.

    Asked of the figures that derive from the submitted page - the ones patched
    in place and the ones supplied whole as a corrected version of it. A figure
    redrawn at the size it prints at draws its own rasters and shares none with
    the submitted page, so the comparison has no meaning there; what replaces it
    is the content gate against a harvest of the script that made the shipped
    page, and the string comparison against that page. Which figure is made
    which way is read from SHIPPED_MANIFEST.csv rather than assumed, and a
    mechanism this does not recognise is reported rather than guessed at in
    either direction.

    Returns (bad, ran). ran is False when the reference tree is absent, so the
    caller can report a skip rather than six failures or, worse, a silent pass.

    The manifest is required whenever the check runs. Without it every redrawn
    figure would be compared as though it had been patched, and every one of its
    rasters counted as lost - a false failure that reads exactly like a real
    one. Treating an absent manifest as "assume patched" is the same silent
    wrong answer in the other direction, so it is refused rather than guessed.
    """
    if not submitted.is_dir():
        print(f"\ncontent-loss check SKIPPED: no submitted figures at "
              f"{submitted}.\n  This tree is not deposited. Point --submitted at "
              f"a copy to run it;\n  the P value and label checks above did run.")
        return [], False
    manifest = find_manifest(root)
    if manifest is None:
        raise SystemExit(
            f"content-loss check cannot run: no {MANIFEST_NAME} beside "
            f"{root.parent} or beside {Path(__file__).resolve().parent}. It "
            f"records which figures are patched and which are redrawn, and "
            f"without it a redrawn figure is misread as one that lost its "
            f"content.")
    import csv as _csv
    with manifest.open(newline="") as fh:
        built_by = {r["figure"]: r["built_by"] for r in _csv.DictReader(fh)}
    print(f"\n{'figure':<12}{'submitted':>11}{'shipped':>9}{'lost':>6}   verdict")
    print("-" * 66)
    bad = []
    for i in range(1, 7):
        src = submitted / f"Figure {i}.pdf"
        dst = root / f"Figure_{i}.pdf"
        if not src.exists() or not dst.exists():
            print(f"{'Figure_' + str(i):<12}{'-':>11}{'-':>9}{'-':>6}   file not found")
            bad.append(f"Figure_{i}")
            continue
        want, got = substantive_images(src), substantive_images(dst)
        how = built_by.get(str(i))
        if how in NOT_COMPARED:
            print(f"{'Figure_' + str(i):<12}{len(want):>11}{len(got):>9}"
                  f"{'n/a':>6}   {NOT_COMPARED[how]}")
            continue
        if how not in COMPARED:
            print(f"{'Figure_' + str(i):<12}{len(want):>11}{len(got):>9}"
                  f"{'-':>6}   built by {how!r}, which this check does not "
                  f"know how to read")
            bad.append(f"Figure_{i}")
            continue
        lost = want - got
        verdict = ("intact" if not lost
                   else f"{len(lost)} image(s) {COMPARED[how]}")
        print(f"{'Figure_' + str(i):<12}{len(want):>11}{len(got):>9}{len(lost):>6}   {verdict}")
        if lost:
            bad.append(f"Figure_{i}")
    return bad, True


def check_framing_labels(root):
    """The retired qualifier appears nowhere; the approved label appears where
    it was approved, exactly as often as it was approved.

    Both halves read the vector text layer, so a figure that extracts no text
    at all is reported rather than passed: absence cannot be asserted of a page
    whose strings were never read.
    """
    print(f"\n{'figure':<12}{'pan-TME':>9}{'approved':>10}   verdict")
    print("-" * 66)
    bad = []
    for i in range(1, 7):
        name = f"Figure_{i}"
        path = root / f"{name}.pdf"
        if not path.exists():
            print(f"{name:<12}{'-':>9}{'-':>10}   file not found")
            bad.append(name)
            continue
        doc = fitz.open(path)
        text = "".join(page.get_text() for page in doc)
        doc.close()
        problems = []
        if not text.strip():
            problems.append("no text layer, so no string can be asserted "
                            "absent from it")
        forbidden = 0
        for f in FRAMING_FORBIDDEN:
            n = text.lower().count(f.lower())
            forbidden += n
            if n:
                problems.append(f"{f} appears {n} times, expected 0")
        approved = 0
        for wanted, count in FRAMING_REQUIRED.get(name, []):
            got = len(re.findall(rf"(?<![A-Za-z0-9_-]){re.escape(wanted)}"
                                 rf"(?![A-Za-z0-9_])", text))
            approved += got
            if got != count:
                problems.append(f"{wanted} appears {got} times, expected {count}")
        wanted_here = " ".join(w for w, _ in FRAMING_REQUIRED.get(name, []))
        verdict = ("; ".join(problems) if problems
                   else (f"approved wording: {wanted_here}" if wanted_here
                         else "retired qualifier absent"))
        print(f"{name:<12}{forbidden:>9}{approved:>10}   {verdict}")
        if problems:
            bad.append(name)
    return bad


def scan(path):
    doc = fitz.open(path)
    text = "".join(page.get_text() for page in doc)
    doc.close()
    return len(STAR.findall(text)), len(PVAL.findall(text)), len(text.strip())


def supplementary_dir(root, given):
    """Where the shipped supplementary figures are, for this layout of --dir."""
    if given:
        return Path(given)
    for cand in (root.parent / "03_Supplementary_Figures",
                 root.parent.parent / "Supplementary_Fixes" / "_patched"):
        if cand.is_dir():
            return cand
    return None


def check_supplementary_labels(root, given):
    """The patched supplementary panels print the new string, never the old."""
    supp = supplementary_dir(root, given)
    bad = []
    for name in sorted(SUPP_LABEL_CHECKS):
        path = supp / f"{name}.pdf" if supp else None
        if path is None or not path.exists():
            print(f"{name[:12]:<12}{'-':>7}{'-':>10}{'-':>10}   file not found"
                  + (f" under {supp}" if supp else "; no supplementary directory "
                     "beside --dir, pass --supp-dir"))
            bad.append(name)
            continue
        doc = fitz.open(path)
        text = "".join(page.get_text() for page in doc)
        doc.close()
        problems = []
        for old, new, count in SUPP_LABEL_CHECKS[name]:
            got = len(re.findall(rf"(?<![A-Za-z0-9_.]){re.escape(new)}(?![0-9])", text))
            if got != count:
                problems.append(f"{new} appears {got} times, expected {count}")
            if old is not None and re.search(
                    rf"(?<![A-Za-z0-9_.]){re.escape(old)}(?![0-9])", text):
                problems.append(f"{old} is still present")
        verdict = "two-sided value" if not problems else "; ".join(problems)
        print(f"{name[:12]:<12}{'n/a':>7}{'n/a':>10}{'n/a':>10}   {verdict}")
        if problems:
            bad.append(name)
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="Main_Figures/_shipped")
    ap.add_argument("--supp-dir", default=None,
                    help="directory of the shipped supplementary figures; by "
                         "default found beside --dir (03_Supplementary_Figures/ "
                         "or Supplementary_Fixes/_patched/)")
    ap.add_argument("--submitted", default=None,
                    help="directory of the submitted figures, named "
                         "'Figure N.pdf'. Defaults to the working tree's "
                         "00_GROUND_TRUTH/figures, which is not deposited.")
    a = ap.parse_args()
    root = Path(a.dir)
    submitted = Path(a.submitted) if a.submitted else SUBMITTED

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

    bad += [n for n in check_framing_labels(root) if n not in bad]

    bad += check_supplementary_labels(root, a.supp_dir)

    lost, content_check_ran = check_no_content_lost(root, submitted)
    bad += [n for n in lost if n not in bad]

    print()
    for name in bad:
        if name in AFFECTED:
            print(f"{name}: panels to check -> " + "; ".join(AFFECTED[name]))
    if bad:
        print("\nRegenerate with: python patch_figure_annotations.py "
              "(main figures) or Supplementary_Fixes/patch_S2_two_sided.py (S2)")
        return 1
    if content_check_ran:
        print("Every converted panel prints its exact two-sided P value, "
              "Figure 4 carries the corrected labels, the framing wording is the "
              "approved one, S2D prints the two-sided value, and no substantive "
              "image was lost to a redaction.")
    else:
        print("Every converted panel prints its exact two-sided P value, "
              "Figure 4 carries the corrected labels, the framing wording is the "
              "approved one and S2D prints the two-sided value. The "
              "content-loss check did not run; see the note above.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
