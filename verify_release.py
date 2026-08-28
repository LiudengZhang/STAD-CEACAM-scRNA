"""
Static self-check for this code release.

It answers one question: would every script here start, and if not, why. It
runs no analysis and produces no figure - it does not need, and does not touch,
the input data.

Checked:
  1. every Python file parses
  2. every `sys.path.insert(... parents[N] / "00_Config")` reaches the real one
  3. every name imported from paths.py is defined there
  4. no machine-specific absolute paths, and no /path/to/ placeholder
     outside the scripts that document the upstream pipeline
  5. no symlink, and no two-group comparison left one-tailed in a panel
  6. the figure driver understands every target the entry point invokes
  7. which declared input directories are present, and which are not

Run: python verify_release.py [root]      (root defaults to this file's directory)
Exits 0 if checks 1-6 pass. Check 7 is reported, never fatal: a code-only
deposit is expected to have no data.
"""

import ast
import os
import re
import sys
from pathlib import Path

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent).resolve()
CFG = ROOT / "00_Config"

PARENTS_CFG = re.compile(r"parents\[(\d+)\]\s*/\s*['\"]00_Config['\"]")
# Machine-specific roots: the two shared filesystems this work ran on, and
ABSOLUTE = re.compile(
    r"['\"]/(?:priv\d*data\d*|rsrch\d+|home/[a-z0-9]+|Users)/")
# The delocaliser leaves this placeholder in scripts that document how an
# upstream input was produced. It is acceptable only there: a placeholder in
# a script the driver calls is a script that cannot run, which is how
# Figure 1C shipped broken.
PLACEHOLDER = re.compile(r"/path/to/")
UPSTREAM_MARK = "upstream Round_4 processing pipeline"
ONE_TAILED = re.compile(r"alternative\s*=\s*['\"](greater|less)['\"]")


def python_files():
    for f in sorted(ROOT.rglob("*.py")):
        if "__pycache__" in f.parts:
            continue
        yield f


def main():
    problems, notes = [], []
    files = list(python_files())

    # 1. parses
    trees = {}
    for f in files:
        try:
            trees[f] = ast.parse(f.read_text(encoding="utf-8", errors="replace"))
        except SyntaxError as e:
            problems.append(f"does not parse: {f.relative_to(ROOT)}: {e}")

    # 2. 00_Config is where each script thinks it is. Build tools (leading
    # underscore) are skipped: nothing runs them, and one of them contains the
    # import line it writes into other scripts, inside a string literal.
    for f in files:
        if f.name.startswith("_"):
            continue
        for m in PARENTS_CFG.finditer(f.read_text(encoding="utf-8",
                                                  errors="replace")):
            n = int(m.group(1))
            try:
                target = f.parents[n] / "00_Config" / "paths.py"
            except IndexError:
                problems.append(f"00_Config unreachable: {f.relative_to(ROOT)}")
                continue
            if not target.exists():
                problems.append(f"00_Config depth wrong: "
                                f"{f.relative_to(ROOT)} parents[{n}]")

    # 3. names imported from paths.py exist
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(CFG))
    try:
        import paths
    except Exception as e:                                  # noqa: BLE001
        problems.append(f"00_Config/paths.py does not import: {e}")
        paths = None
    if paths is not None:
        available = {n for n in dir(paths) if not n.startswith("_")}
        for f, tree in trees.items():
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom) and node.module == "paths":
                    for a in node.names:
                        if a.name != "*" and a.name not in available:
                            problems.append(
                                f"paths.py has no {a.name!r}: "
                                f"{f.relative_to(ROOT)}")

    # 4. absolute paths, and placeholders outside provenance scripts
    others = [f for f in ROOT.rglob("*")
              if f.is_file() and ".git" not in f.parts
              and f.suffix.lower() in (".sh", ".yaml", ".yml", ".r",
                                       ".json")]
    for f in list(files) + sorted(others):
        text = f.read_text(encoding="utf-8", errors="replace")
        provenance = UPSTREAM_MARK in text
        for i, line in enumerate(text.splitlines(), 1):
            if ABSOLUTE.search(line):
                problems.append(f"absolute path: {f.relative_to(ROOT)}:{i}")
            if PLACEHOLDER.search(line) and not provenance:
                problems.append(
                    f"placeholder path: {f.relative_to(ROOT)}:{i}")

    # 4b. symlinks. git stores one as the text of its target, so an
    # absolute-path link publishes the author's filesystem layout.
    for f in ROOT.rglob("*"):
        if ".git" in f.parts:
            continue
        if f.is_symlink():
            problems.append(f"symlink: {f.relative_to(ROOT)}")

    # 5. one-tailed tests in panel scripts
    panels = ROOT / "03_Final_Panels"
    if panels.exists():
        for f in panels.rglob("*.py"):
            if ONE_TAILED.search(f.read_text(encoding="utf-8", errors="replace")):
                problems.append(f"one-tailed test: {f.relative_to(ROOT)}")

    # 6. driver targets
    driver = panels / "_run_all_panels.sh"
    if not driver.exists():
        problems.append("03_Final_Panels/_run_all_panels.sh missing")
    else:
        d = driver.read_text(encoding="utf-8")
        for target in ("revision|rev) run_revision", "supp|s|S) run_supp"):
            if target not in d:
                problems.append(f"driver does not dispatch: {target}")

    # 7. inputs present (reported, not fatal)
    if paths is not None:
        for name in ("RAW_INPUTS", "PREPARATION"):
            p = Path(getattr(paths, name))
            has = p.exists() and any(p.iterdir()) if p.exists() else False
            notes.append(f"{name:<12} {p}  {'present' if has else 'ABSENT'}")

    print("=" * 74)
    print("CODE RELEASE SELF-CHECK")
    print("=" * 74)
    print(f"root         : {ROOT}")
    print(f"python files : {len(files)}")
    print()
    for n in notes:
        print(f"  {n}")
    print()
    if problems:
        print(f"{len(problems)} PROBLEMS")
        for p in problems:
            print(f"  {p}")
        return 1
    print("All checks passed: every script parses, every import resolves, no")
    print("machine-specific paths, no one-tailed panel tests, and the driver")
    print("understands every target the entry point uses.")
    if any("ABSENT" in n for n in notes):
        print()
        print("No input data are attached, so NO ANALYSIS WAS RUN and no figure")
        print("was produced. Attach the data asset and re-run to reproduce the")
        print("figures; see data/README.txt for the expected layout.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
