"""
WP3 / Reviewer 1 point R1.3c - convert the affected panel scripts.

Copies every panel script that used a one-tailed test into
03_Final_Panels/, switches the test to two-sided, and replaces
star-based significance annotation with the exact P value, as the reviewer
asked ("report exact P values rather than only significance thresholds").

Only the statistics and the annotation change; layout, colours and data paths
are untouched, so the revised panels remain drop-in replacements.

Run:  python patch_panels_two_sided.py          (writes the patched copies)
      python patch_panels_two_sided.py --check  (reports without writing)
"""

from pathlib import Path
import re
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import FINAL_PANELS, PREPARATION, REVISED_PANELS  # noqa: E402

DEST = REVISED_PANELS / "Main_Figures"

# Panel scripts that carried a one-tailed test, from twosided_sweep.csv.
PANEL_SCRIPTS = [
    "02_Figure_2/02_J/create_ceacam6_boxplot.py",
    "02_Figure_2/02_K/create_ceacam5_boxplot.py",
    "02_Figure_2/02_N/create_ceacam6_prjeb25780_boxplot.py",
    "02_Figure_2/02_O/create_ceacam5_prjeb25780_boxplot.py",
    "02_Figure_2/02_M/create_ihc_combined_boxplot.py",
    "02_Figure_2/02_G/create_mp45_horizontal_boxplot.py",
    "03_Figure_3/03_H/create_spatial_boxplots.py",
    "03_Figure_3/03_I/create_spatial_boxplot_stroma.py",
    "03_Figure_3/03_J/create_spatial_boxplot_immune.py",
    "05_Figure_5/05_I/create_cd274_momac_boxplot.py",
    "05_Figure_5/05_J/create_cd274_epithelial_boxplot.py",
    "05_Figure_5/05_K/create_cd274_fibroblast_boxplot.py",
    "05_Figure_5/05_DC_CD274/create_cd274_dc_boxplot.py",
    "05_Figure_5/05_TF/create_tf_4group_panels.py",
    # Figure 5L. Missed on the first pass and caught by the leftover scan below;
    # the text quotes no P value for this panel, so only the panel changes.
    "05_Figure_5/05_IL6_CD4/create_il6_stat3_cd4t_boxplot.py",
]

# Preparation scripts outside 03_Final_Panels that also used one-tailed tests.
PREP_SCRIPTS = [
    "IHC/quantify_ceacam_ihc.py",
    "BayesPrism/step3_plot_results.py",
    "Metaprogram_Permutation/mp4_pre_r_permutation_analysis.py",
]

# 1. the tests themselves
TEST_SUBS = [
    (re.compile(r"""alternative\s*=\s*['"]greater['"]"""), "alternative='two-sided'"),
    (re.compile(r"""alternative\s*=\s*['"]less['"]"""), "alternative='two-sided'"),
]

# 2. star annotation -> exact P. Each pattern rewrites the whole assignment.
STAR_SUBS = [
    (re.compile(
        r"""(\w+)\s*=\s*'\*\*\*'\s*if\s*(\w+)\s*<\s*0\.001\s*else\s*'\*\*'"""
        r"""\s*if\s*\2\s*<\s*0\.01\s*else\s*'\*'\s*if\s*\2\s*<\s*0\.05\s*else\s*'ns'"""),
     r"\1 = f'P = {\2:.3f}' if \2 >= 0.001 else 'P < 0.001'"),
]

# 2b. Text that appears ON the figure or in the console describing the tail.
# A panel that plots a two-sided P but is captioned "one-tailed" is worse than
# the original, so these are rewritten wherever they occur.
CAPTION_SUBS = [
    (re.compile(r"Wilcoxon signed-rank \(one-tailed\)"),
     "Wilcoxon signed-rank (two-sided)"),
    (re.compile(r"Mann-Whitney U \(one-tailed, NR > R\)"),
     "Mann-Whitney U (two-sided)"),
    (re.compile(r"One-tailed Mann-Whitney U test \(NR > R\)"),
     "Two-sided Mann-Whitney U test"),
    (re.compile(r"Mann-Whitney U \(NR > R\)"), "Mann-Whitney U (two-sided)"),
    (re.compile(r"# One-tailed Mann-Whitney U \(NR > R\)"),
     "# Two-sided Mann-Whitney U"),
    (re.compile(r"One-tailed Mann-Whitney U test \(NR > R\), n=4 R, n=4 NR"),
     "Two-sided Mann-Whitney U test, n=4 R, n=4 NR"),
    (re.compile(r"Mann-Whitney \(one-sided: NR > R\)"),
     "Mann-Whitney U (two-sided)"),
    (re.compile(r"P-value \(one-sided NR>R\)"), "P-value (two-sided)"),
    # Left behind by the substitution two rules above, which fires first.
    (re.compile(r"# One-tailed Mann-Whitney U \(two-sided\)"),
     "# Two-sided Mann-Whitney U"),
]

# 2c. A machine-specific path that would ship in the code release. The IHC
# thumbnails already have a name in paths.py; the script predates it.
PATH_SUBS = [
    (re.compile(
        r'THUMB_DIR = Path\("/[^"]*Experimental_Data[^"]*"\s*\n\s*"[^"]*"\)'),
     'import sys\n'
     '# Depth is that of the code release, 02_Preparation_for_Panels/IHC/.\n'
     'sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))\n'
     'from paths import IHC_THUMBNAILS  # noqa: E402\n\n'
     'THUMB_DIR = IHC_THUMBNAILS'),
    # Panel 2M builds its own path to the IHC table instead of taking the name
    # paths.py already defines, so it never sees an overridden input root and
    # fails wherever the prepared results are not inside the code tree.
    (re.compile(
        r'IHC_CSV = Path\(__file__\)\.resolve\(\)\.parents\[3\] / '
        r'"02_Preparation_for_Panels" / "IHC" / '
        r'"ceacam_ihc_color_deconv_results\.csv"'),
     'from paths import IHC_COLOR_DECONV_CSV  # noqa: E402\n'
     'IHC_CSV = IHC_COLOR_DECONV_CSV'),
]

# 3. multi-line if/elif star ladders (spatial panels)
LADDER = re.compile(
    r"""(?P<ind>[ \t]*)if\s+(?P<p>\w+)\s*<=?\s*0\.0001:\s*\n"""
    r"""(?:[ \t]*(?P<v>\w+)\s*=\s*'\*\*\*\*'\s*\n)"""
    r"""(?:[ \t]*elif\s+(?P=p)\s*<=?\s*0\.001:\s*\n[ \t]*(?P=v)\s*=\s*'\*\*\*'\s*\n)"""
    r"""(?:[ \t]*elif\s+(?P=p)\s*<=?\s*0\.01:\s*\n[ \t]*(?P=v)\s*=\s*'\*\*'\s*\n)"""
    r"""(?:[ \t]*elif\s+(?P=p)\s*<=?\s*0\.05:\s*\n[ \t]*(?P=v)\s*=\s*'\*'\s*\n)"""
    r"""(?:[ \t]*else:\s*\n[ \t]*(?P=v)\s*=\s*'ns'\s*\n)""")

HEADER = (
    "# REVISED FOR CIR-26-0753-ET, reviewer 1 point R1.3c:\n"
    "# the test is two-sided and the annotation reports the exact P value.\n"
    "# Original one-tailed version: submission-tree/03_Final_Panels/{rel}\n")


def patch_text(text, rel):
    n_test = n_star = 0
    for pat, rep in TEST_SUBS:
        text, k = pat.subn(rep, text)
        n_test += k
    for pat, rep in STAR_SUBS:
        text, k = pat.subn(rep, text)
        n_star += k

    def ladder_rep(m):
        return (f"{m.group('ind')}{m.group('v')} = "
                f"f'P = {{{m.group('p')}:.3f}}' if {m.group('p')} >= 0.001 "
                f"else 'P < 0.001'\n")

    text, k = LADDER.subn(ladder_rep, text)
    n_star += k

    # Captions run before the provenance header is inserted below, so the
    # header's own "one-tailed version" wording is never rewritten.
    for pat, rep in CAPTION_SUBS:
        text, k = pat.subn(rep, text)
        n_star += k

    for pat, rep in PATH_SUBS:
        text, k = pat.subn(rep, text)
        n_star += k

    if n_test:
        lines = text.split("\n")
        insert_at = 0
        if lines and lines[0].startswith("#!"):
            insert_at = 1
        lines.insert(insert_at, HEADER.format(rel=rel).rstrip("\n"))
        text = "\n".join(lines)
    return text, n_test, n_star


def main():
    check = "--check" in sys.argv
    total_t = total_s = 0
    touched = []

    for rel in PANEL_SCRIPTS:
        src = FINAL_PANELS / rel
        if not src.exists():
            print(f"  MISSING {src}")
            continue
        text = src.read_text(encoding="utf-8")
        new, nt, ns = patch_text(text, rel)
        total_t += nt
        total_s += ns
        if nt or ns:
            touched.append((rel, nt, ns))
        if not check and (nt or ns):
            dst = DEST / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_text(new, encoding="utf-8")

    for rel in PREP_SCRIPTS:
        src = PREPARATION / rel
        if not src.exists():
            print(f"  MISSING {src}")
            continue
        text = src.read_text(encoding="utf-8")
        new, nt, ns = patch_text(text, f"../02_Preparation_for_Panels/{rel}")
        total_t += nt
        total_s += ns
        if nt or ns:
            touched.append((f"prep/{rel}", nt, ns))
        if not check and (nt or ns):
            dst = DEST / "_preparation" / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_text(new, encoding="utf-8")

    print(f"{'Would patch' if check else 'Patched'} {len(touched)} scripts:")
    for rel, nt, ns in touched:
        print(f"   {rel:<58} tests={nt}  annotations={ns}")
    print(f"\nTotal: {total_t} one-tailed tests converted, "
          f"{total_s} star annotations replaced with exact P.")
    if not check:
        print(f"Written to {DEST}")

    # Nothing should be left one-tailed in the revised tree.
    if not check:
        # Archived scripts are not part of any current figure and are left as
        # they are, so the check that matters is over the live ones.
        # Two things have to be true afterwards, and the first version of this
        # audit only checked the first: no test is still one-tailed, AND no
        # panel still labels its result with a star. Figure 5D/E passed the
        # one-tailed check while continuing to print '*', because the star
        # ladder there is written across two lines and the substitution regex
        # is single-line. R1.3c asks for exact P values, so the label matters
        # as much as the test.
        patterns = (
            ("one-tailed test", re.compile(r"""alternative\s*=\s*['"](greater|less)['"]""")),
            ("star label", re.compile(r"""['"]\*\*\*['"]\s*if""")),
        )
        for what, rx in patterns:
            leftover, archived = [], []
            for f in DEST.rglob("*.py"):
                if rx.search(f.read_text(encoding="utf-8")):
                    (archived if "_archived" in f.parts or
                     f.name.startswith("_") else leftover).append(f)
            print(f"Remaining {what}s in live panel scripts: {len(leftover)}")
            for f in leftover:
                print(f"   {f.relative_to(DEST)}")
            if archived:
                print(f"({len(archived)} archived/exploratory scripts still carry "
                      f"a {what}; they generate no current figure)")


if __name__ == "__main__":
    main()
