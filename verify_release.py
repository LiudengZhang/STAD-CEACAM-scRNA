# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
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
     outside the scripts that say why they carry one - the ones that document
     the upstream pipeline, and the ones whose interpreter and library paths
     under a home directory were delocalised
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
UPSTREAM_MARK = "upstream processing pipeline"
# The second reason a file may legitimately carry a placeholder: the
# home-directory rewrite. update_release.py prepends
# LOCAL_HOME_NOTE, which contains this phrase, wherever that rewrite fires and
# a '#' comment line is legal. It is a separate mark from UPSTREAM_MARK on
# purpose: the four shell drivers it applies to are the revision analyses' own
# drivers, and claiming they document the upstream pipeline to buy a pass here
# would be a false statement in a deposited file. Kept in step with
# update_release.py:LOCAL_HOME_MARK.
LOCAL_HOME_MARK = "paths under the author's home directory"
# The third reason: the upstream pipelines. Those ran on
# two named filers and on a cluster login node, and update_release.py's
# MACHINE_ROOT_NOTE says so. It is separate from the other two marks for the
# same reason they are separate from each other: each has to be true of the file
# it is prepended to, and neither "documents the upstream pipeline" nor "paths
# under the author's home directory" describes /path/to/machine or /path/to/machine. Kept in
# step with update_release.py:MACHINE_ROOT_MARK.
MACHINE_ROOT_MARK = "specific to the machines this pipeline"
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

    # 2b. every MODULE-LEVEL import resolves to something the release carries.
    #
    # This exists because of what it would have caught and did not. Seven
    # shipped drivers began with
    #     sys.path.insert(0, str(ROOT / "10_Reproduction"))
    #     import compare_panel_content
    # and 10_Reproduction/ is a verification tree the release deliberately does
    # not ship. Every one of them died at that import on every capsule run,
    # while this file reported "every import resolves" and exited 0 - because
    # nothing here looked at imports other than `from paths import ...`, and
    # check 2 skips any file whose name begins with "_" on the reasoning that
    # nothing runs them. _driver_base.py is imported by eight scripts on every
    # run, so that reasoning was false and the "_" exemption does not apply
    # here.
    #
    # A name is accepted if it is the standard library, or a module the release
    # itself carries, or a dependency environment.yml DECLARES. Declared, not
    # installed: the check must give the same answer on a machine that happens
    # to have a package as on one that does not, or it measures the host
    # instead of the record.
    #
    # Scope is the code a reviewer runs. upstream/ ships as a record of what
    # produced the deposited intermediates and the entry point never calls it.
    import sysconfig                                            # noqa: PLC0415

    def _declared_dependencies():
        names = set()
        for y in ROOT.rglob("environment.yml"):
            for line in y.read_text(encoding="utf-8",
                                    errors="replace").splitlines():
                line = line.strip().lstrip("-").strip()
                if not line or line.startswith("#") or line.endswith(":"):
                    continue
                pkg = re.split(r"[=<>!\[ ]", line, 1)[0].strip()
                if pkg:
                    names.add(pkg.lower().replace("-", "_"))
        # Import name differs from distribution name for these.
        # Import name differs from distribution name for these; kept
        # lowercase because every comparison below is lowercased.
        names.update({"sklearn", "skimage", "yaml", "pil", "cv2", "fitz",
                      "dateutil", "pkg_resources", "setuptools", "attr",
                      "mpl_toolkits", "importlib_metadata",
                      "typing_extensions", "pertpy"})
        return names

    RUNNABLE = ("03_Final_Panels", "04_Revision_Analyses", "00_Config")
    declared = _declared_dependencies()
    shipped_mods = {f.stem for f in files} | {
        d.name for d in ROOT.rglob("*") if d.is_dir()
        and (d / "__init__.py").exists()}
    stdlib = set(getattr(sys, "stdlib_module_names", ())) | set(
        sysconfig.get_config_vars().get("TZPATH", "").split(":"))
    for f, tree in trees.items():
        rel = f.relative_to(ROOT)
        if not rel.parts or rel.parts[0] not in RUNNABLE:
            continue
        # `work/` is a record of how an output was produced, like upstream/.
        # _run_all_panels.sh globs 04_Revision_Analyses/*/scripts/*.py and
        # never reaches it, so an import there cannot break a reviewer's run.
        # Measured while writing this check: work/build_pseudobulk.py imports
        # h5py, which environment.yml does not declare - it arrives
        # transitively with anndata. Reported, not silently accepted: if that
        # script ever moves under scripts/, this check will say so.
        if "work" in rel.parts:
            continue
        for node in tree.body:                        # module level only
            if isinstance(node, ast.Import):
                tops = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                tops = ([node.module.split(".")[0]]
                        if node.level == 0 and node.module else [])
            else:
                continue
            for t in tops:
                low = t.lower().replace("-", "_")
                if (t in stdlib or low in declared or t in shipped_mods
                        or low in {m.lower() for m in shipped_mods}):
                    continue
                problems.append(
                    f"module-level `import {t}` resolves to nothing the "
                    f"release ships and environment.yml does not declare it: "
                    f"{rel}")

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
        explained = (UPSTREAM_MARK in text or LOCAL_HOME_MARK in text
                     or MACHINE_ROOT_MARK in text)
        for i, line in enumerate(text.splitlines(), 1):
            if ABSOLUTE.search(line):
                problems.append(f"absolute path: {f.relative_to(ROOT)}:{i}")
            if PLACEHOLDER.search(line) and not explained:
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
