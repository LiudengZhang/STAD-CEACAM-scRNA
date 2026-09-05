"""
Gate PROVENANCE.csv against the files it describes.

The failure this exists to prevent: a panel script whose output disagrees with
the printed figure, sitting in the tree with nothing to say so, until someone
concludes the paper is wrong. That has happened twice, both times on Figure 5A,
both times retracted. See 00_GROUND_TRUTH/README.md.

Eight checks.

  1. Structure      every row names a directory and a script that exist.
  2. Coverage       every printed panel of Figures 1-6 and S7-S11 has a row, and
                    every panel directory is either mapped to a printed panel or
                    recorded as an orphan.
  3. Known broken   a row marked reproduces_published = no must carry a
                    KNOWN_BROKEN.md beside its script and a warning in the
                    script's own docstring. A panel may be broken; it may not be
                    silently broken.
  4. Placement      a row marked yes that carries printed_rect_mm must actually
                    appear at that rectangle of the shipped figure: the panel's
                    own text tokens are compared against the text inside the
                    rectangle of Main_Figures/_patched/Figure_N.pdf. This is the
                    check that a whole-page comparison cannot make - a text
                    sweep against the whole figure scores Figure 5A at 0 percent
                    missing, because "Hypoxia" appears in panel N.
  5. Guard rails    CLAUDE.md, 00_GROUND_TRUTH/README.md and the two tree
                    READMEs exist, and the ground-truth figures are read-only.
  6. Freshness      a verdict must be newer than the script it judges. A panel
                    edited after it was verified carries a verdict about a file
                    that no longer exists, which reads exactly like a verdict
                    about the file that does.
  7. Lookup         every panel directory that holds a panel script or a drawn
                    artefact is named in some row's source_dir. CLAUDE.md's
                    second rule says to look every panel up here rather than
                    read its letter off the directory name; a directory that
                    cannot be looked up leaves no way to obey that but guessing.
                    Check 2 already does this for the main figures by directory
                    name; this extends it to the supplementary trees and matches
                    on the path the rows actually carry.
  8. No unknowns    reproduces_published must be yes, no or na. "unknown" is not
                    an answer to whether a panel reproduces - it is the absence
                    of one, and it passes every other check in this file while
                    saying nothing at all.
  9. Inputs         every data file a panel script names must exist, and must not
                    live in an archive, a backup or a temp workspace. Figure 5A
                    was recorded not reproducible for a month because its input
                    had been quietly repointed from GSEA/MoMac_mast_prerank_gsea.csv
                    to a different run under GSEA/post/, and nothing in the tree
                    looked at where a script reads from. Scripts this cannot
                    resolve anything from are listed rather than passed over: a
                    blind spot that reports itself is a blind spot, one that stays
                    quiet reads as a clean bill of health.

                    Two blind spots were closed on 2026-09-01, both of them
                    shapes that let the Figure 5A repoint through:

                    a. Directory reads. A script that names its input with
                       .glob("*.csv"), .iterdir(), os.listdir(), a formatted
                       path or a directory variable joined at run time named no
                       file literal, so nothing was looked at. A whole archived
                       GSEA run could be read a file at a time in silence. Such
                       reads are now resolved against the filesystem where the
                       directory and the pattern are both static, and where
                       either is not, the read is reported as UNRESOLVED. An
                       input this check cannot evaluate has to be visible; it
                       must not be absent.

                    b. Same-directory reads. Any data path resolving inside a
                       script's own directory was assumed to be something the
                       script writes, and skipped. That is wrong for a panel
                       script that reads a prepared table sitting beside it: a
                       stale or missing one was invisible. Read and write are
                       now told apart by how the path is used - read_csv,
                       read_table, open(..., "r") against to_csv, savefig,
                       open(..., "w") - and not by where it sits. A path only
                       written, or written and then read back, is still an
                       output.

                    Positive controls for both live in the session scratch, not
                    in the tree: a script that globs
                    07_Archive/.../12_R1.8_DEG_Recompute/outputs/gsea, and one
                    that reads a nonexistent prepared table beside itself. Both
                    passed this check in silence before the change.
 11. Upstream       the same archive guard, applied to the pipelines that produce
     inputs         what the panels read. Check 9 looks only at panel scripts, so
                    an input repointed one step upstream - in a NicheNet or NMF
                    config, in a BayesPrism step script - passed it untouched.
                    That is the Figure 5A fault with one more link in the chain.
                    Two roots are read: 02_Preparation_for_Panels/ (the five imported
                    pipelines) and the Round_5 preparation tree whose configs
                    they reference.

                    Only the archive guard is applied here, never the existence
                    check. 02_Preparation_for_Panels/ holds scripts and no outputs on purpose
                    - import_pipelines.py says so in its own docstring, because
                    a second copy of an intermediate can drift from the one the
                    panels read - so every relative output path in it is
                    legitimately absent and 81 of them would fail an existence
                    check that means nothing.

                    Python is read with _inputs_of/_Resolver, exactly as check 9
                    reads a panel script. YAML, R and shell carry their paths as
                    plain text and are scanned for absolute paths. A file that
                    itself lives in an archive is skipped and counted, not
                    failed: an archived config naming archived data is history,
                    not a live repoint.

 10. Broken expires a row marked reproduces_published = no must name a mechanism
                    in its note - which input, which run, which step - and its
                    KNOWN_BROKEN.md must carry a "Retest-by: YYYY-MM-DD" line that
                    has not passed. The marker exists so a machine cannot overwrite
                    a human finding; it must not also let a wrong finding harden
                    into a fact. The Figure 5A marker did exactly that.

Run: python verify_panel_provenance.py
Exit status is non-zero if anything fails.
"""

from pathlib import Path
import ast
import csv
import datetime
import re
import sys
import unicodedata

try:
    import fitz
except ImportError:
    sys.exit("PyMuPDF (fitz) is required: conda run -n Liudeng_Python_310 ...")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
GROUND = ROOT / "00_GROUND_TRUTH"
PATCHED = HERE / "Main_Figures" / "_patched"
MANIFEST = HERE / "PROVENANCE.csv"
RELEASE = ROOT / "05_Code_Release" / "github_repo" / "03_Final_Panels"

MM = 72 / 25.4
# All five figures resolve into the working tree. Figures 1 and 4 pointed at
# RELEASE until 2026-09-01, because they were the two with no working-tree panel
# directories - update_release.py says so in its own docstring, and filled them
# from Round_5. So this file was checking a *generated* copy of eleven scripts:
# a stale release would have been validated as though it were the source, and
# the check would have reported "all inputs live" about code nobody was editing.
# Both directories were copied in from Round_5/03_Final_Panels that day, verified
# byte-identical file by file, so nothing this file judges changed value.
FIG_DIR = {"1": HERE / "Main_Figures" / "01_Figure_1",
           "2": HERE / "Main_Figures" / "02_Figure_2",
           "3": HERE / "Main_Figures" / "03_Figure_3",
           "4": HERE / "Main_Figures" / "04_Figure_4",
           "5": HERE / "Main_Figures" / "05_Figure_5"}
PRINTED_PANELS = {"1": "ABC", "2": "ABCDEFGHIJKLMN", "3": "ABCDEFGHIJKLMN",
                  "4": "ABCDEFGH", "5": "ABCDEFGHIJKLMN"}
BROKEN_MARKER = "KNOWN_BROKEN.md"
RETEST = re.compile(r"^Retest-by:\s*(\d{4})-(\d{2})-(\d{2})\s*$", re.M)
# Suffixes that make a path a data input rather than a drawn artefact.
DATA_SUFFIX = {".csv", ".tsv", ".txt", ".h5ad", ".pkl", ".xlsx", ".xls",
               ".loom", ".rds", ".json", ".gmt", ".parquet", ".npy", ".npz"}
# Directory names that mean "this is not a live input".
STALE_MARKERS = ("_archived", "_Archived", "99_Archive", "98_Temp_workspace",
                 "backup", "_BACKUP", "Backup", "_Old", "_old", "07_Archive",
                 "release_candidate")
# Check 11. The upstream pipelines that produce what the panels read.
UPSTREAM = ROOT / "02_Preparation_for_Panels"
# An absolute path written into a config or a non-Python pipeline step. The two
# prefixes are the ones this project's storage actually uses; LOCAL_PATH in
# verify_release_sync.py names the same pair.
ABS_PATH = re.compile(r"/(?:priv18data1|rsrch\d+)/[A-Za-z0-9_./+~-]+")
# YAML has no AST here and R and shell have no Python one; their paths are text.
TEXT_INPUT_PATTERNS = ("*.yaml", "*.yml", "*.R", "*.sh")
# How a path is used tells read from write. Location does not: a table prepared
# beside a panel script and read by it is an input, and check 9 skipped every
# one of those as an output until 2026-09-01.
READ_CALLS = {"read_csv", "read_table", "read_excel", "read_tsv", "read_fwf",
              "read_parquet", "read_json", "read_pickle", "read_hdf",
              "read_feather", "read_stata", "read_sas", "read_h5ad",
              "read_10x_h5", "read_10x_mtx", "read_loom", "read_mtx",
              "read_h5mu", "read_gmt", "read_umi_tools", "loadtxt",
              "genfromtxt", "load", "load_npz", "read_visium"}
WRITE_CALLS = {"to_csv", "to_excel", "to_parquet", "to_pickle", "to_json",
               "to_hdf", "to_feather", "to_stata", "savefig", "savetxt",
               "save", "savez", "savez_compressed", "write_h5ad", "write_csv",
               "write_loom", "write", "dump", "tofile", "makedirs"}
# Methods that read or write the path they are called on, not an argument.
READ_METHODS = {"read_text", "read_bytes"}
WRITE_METHODS = {"write_text", "write_bytes", "mkdir", "touch", "unlink"}
# Calls that read a whole directory rather than a named file.
DIR_METHODS = {"glob", "rglob", "iterdir"}
READ, WRITE = "read", "write"
# Patterns that name a drawn artefact rather than a data input; a directory of
# these is figure assembly, which NO_DATA_INPUTS already accounts for.
ARTEFACT_PATTERNS = tuple(f"*{s}" for s in (".pdf", ".svg", ".png", ".jpg",
                                            ".jpeg", ".tif", ".tiff", ".eps"))
# Scripts that draw nothing from a data file - they compose finished PDFs.
NO_DATA_INPUTS = {"assemble_figure_2.py", "assemble_figure_3.py",
                  "assemble_figure_5.py", "assemble_new_supplementaries.py",
                  "build_figure.py", "compare_orientations.py",
                  "create_ihc_representative_2x2.py", "create_gsea_4panels.py"}
# What makes a directory a panel directory rather than a bag of inputs.
ARTEFACT_SUFFIXES = {".pdf", ".svg", ".png"}
# Words a panel may carry that the printed figure legitimately renders as
# outlines or drops: the panel letter is redrawn by the patcher, and matplotlib
# writes minus signs the printed file does not use.
TOKEN_MIN_LEN = 4

failures, notes = [], []


def fail(msg):
    failures.append(msg)


def note(msg):
    notes.append(msg)


def toks(text):
    text = unicodedata.normalize("NFKC", text)
    text = (text.replace("−", "-").replace("α", "alpha")
                .replace("κ", "k").replace("ρ", "rho"))
    return {w.lower() for w in re.findall(r"[A-Za-z][A-Za-z0-9/+-]{%d,}" %
                                          (TOKEN_MIN_LEN - 1), text)}


def resolve(figure, rel):
    """
    Absolute path for a source_dir / source_script / source_file entry.

    Entries are relative to the figure's own panel root for the main figures,
    and to this directory for everything else. Figure 1's panel A is the one
    exception: it was redrawn this round and lives in the working tree, not in
    the generated release where the rest of Figure 1's panels sit.
    """
    if not rel or rel.startswith("<"):
        return None
    if rel.startswith("_panel_1A"):
        return HERE / "Main_Figures" / rel
    if rel.startswith(("00_GROUND_TRUTH/", "04_Revision_Analyses/")):
        # Twenty-one supplementary panels are written by an analysis module
        # rather than by a script beside them.
        return ROOT / rel
    if rel.startswith(("Supplementary_New/", "Supplementary_Fixes/")):
        return HERE / rel
    return FIG_DIR.get(figure, HERE) / rel


def panel_dir(row):
    """Absolute directory for a row's first source_dir, or None."""
    return resolve(row["figure"], row["source_dir"].split("; ")[0])


def load():
    if not MANIFEST.exists():
        sys.exit(f"{MANIFEST} not found - run build_provenance.py first")
    with open(MANIFEST) as fh:
        return list(csv.DictReader(fh))


def check_structure(rows):
    for r in rows:
        if r["build_path"] in ("schematic", "carried_over"):
            continue
        d = panel_dir(r)
        tag = f"Figure {r['figure']} panel {r['printed_panel'] or '-'}"
        if d is None:
            fail(f"{tag}: no source_dir")
            continue
        if not d.is_dir():
            fail(f"{tag}: source_dir does not exist: {d}")
            continue
        for rel in filter(None, r["source_script"].split("; ")):
            p = resolve(r["figure"], rel)
            if p is None or not p.exists():
                fail(f"{tag}: source_script missing: {p or rel}")


def check_coverage(rows):
    seen = {(r["figure"], r["printed_panel"]) for r in rows}
    for fig, letters in PRINTED_PANELS.items():
        for c in letters:
            if (fig, c) not in seen:
                fail(f"Figure {fig} panel {c} has no row in PROVENANCE.csv")
    for s in sorted((HERE / "Supplementary_New").glob("S*/")):
        figname = s.name.split("_")[0]
        for d in sorted(p for p in s.iterdir()
                        if p.is_dir() and not p.name.startswith("_")):
            letter = d.name.split("_")[-1]
            if (figname, letter) not in seen:
                fail(f"{figname} panel {letter} has no row in PROVENANCE.csv")
    for fig, root in FIG_DIR.items():
        if not root.is_dir():
            continue
        mapped = {p for r in rows if r["figure"] == fig
                  for p in r["source_dir"].split("; ")}
        for d in sorted(p.name for p in root.iterdir()
                        if p.is_dir() and not p.name.startswith("_")):
            if d not in mapped:
                fail(f"Figure {fig}: directory {d} is in neither a panel row "
                     f"nor an orphan row")


def check_known_broken(rows):
    for r in rows:
        if r["reproduces_published"] != "no":
            continue
        tag = f"Figure {r['figure']} panel {r['printed_panel']}"
        d = panel_dir(r)
        if d is None or not (d / BROKEN_MARKER).exists():
            fail(f"{tag} is marked reproduces_published=no but has no "
                 f"{BROKEN_MARKER} beside its script - a broken panel must not "
                 f"be silently broken")
            continue
        warned = False
        for rel in filter(None, r["source_script"].split("; ")):
            p = resolve(r["figure"], rel)
            if p and p.exists() and BROKEN_MARKER in p.read_text(errors="replace"):
                warned = True
        if not warned:
            fail(f"{tag}: {BROKEN_MARKER} exists but no source_script points at "
                 f"it from its own docstring")
        else:
            note(f"{tag}: recorded as not reproducible, marker and warning both "
                 f"present")


def check_placement(rows):
    for r in rows:
        rect = r["printed_rect_mm"]
        if r["reproduces_published"] != "yes" or not rect:
            continue
        fig = r["figure"]
        tag = f"Figure {fig} panel {r['printed_panel']}"
        shipped = PATCHED / f"Figure_{fig}.pdf"
        if not shipped.exists():
            fail(f"{tag}: {shipped} not found - run patch_figure_annotations.py")
            continue
        d = panel_dir(r)
        src = None
        for cand in (d.glob("*.pdf") if d and d.is_dir() else []):
            if r["source_file"] and cand.stem == Path(r["source_file"]).stem:
                src = cand
        if src is None:
            stems = [c for c in (d.glob("*.pdf") if d and d.is_dir() else [])]
            src = stems[0] if len(stems) == 1 else None
        if src is None:
            fail(f"{tag}: cannot identify the panel PDF in {d}")
            continue

        x0, y0, x1, y1 = (float(v) for v in rect.split(","))
        page = fitz.open(shipped)[0]
        in_rect = toks(page.get_text(clip=fitz.Rect(x0 * MM, y0 * MM,
                                                    x1 * MM, y1 * MM)))
        panel = toks("\n".join(p.get_text() for p in fitz.open(src)))
        if not panel:
            note(f"{tag}: panel carries no text layer; placement not checked")
            continue
        missing = sorted(panel - in_rect)
        if missing:
            fail(f"{tag}: {len(missing)} of {len(panel)} words in {src.name} are "
                 f"absent from the shipped figure at {rect} mm - the panel in "
                 f"the paper is not this file. Missing: {', '.join(missing[:8])}")
        else:
            note(f"{tag}: all {len(panel)} words of {src.name} found inside "
                 f"{rect} mm of Figure_{fig}.pdf")


def check_lookup(rows):
    """Every panel directory must be findable in PROVENANCE.csv.

    Rows name the main-figure directories by bare name (02_D) and the
    supplementary ones by path (Supplementary_New/S9_Mechanism_Specificity/S9_C),
    so each tree is matched on the form its own rows use.
    """
    mapped = {p.strip() for r in rows
              for p in r["source_dir"].split(";") if p.strip()}
    checked = 0
    def scan(root, by_path):
        nonlocal checked
        if not root.is_dir():
            return
        for d in sorted(p for p in root.iterdir()
                        if p.is_dir() and not p.name.startswith("_")):
            if not any(f.suffix.lower() in ARTEFACT_SUFFIXES
                       or f.name.startswith(("create_", "generate_"))
                       for f in d.iterdir() if f.is_file()):
                continue
            checked += 1
            key = str(d.relative_to(HERE)) if by_path else d.name
            if key not in mapped:
                fail(f"{d.relative_to(HERE)} produces a panel artefact but no "
                     f"row names it in source_dir - there is no way to look its "
                     f"printed letter up, and the directory name is not it")
    for root in sorted((HERE / "Main_Figures").glob("0*_Figure_*")):
        scan(root, by_path=False)
    for root in sorted((HERE / "Supplementary_New").glob("S*")):
        scan(root, by_path=True)
    scan(HERE / "Supplementary_Fixes", by_path=True)
    note(f"{checked} panel directories are named in PROVENANCE.csv")


def check_no_unknowns(rows):
    """A row that answers "unknown" answers nothing."""
    for r in rows:
        if r["reproduces_published"] != "unknown":
            continue
        fail(f"Figure {r['figure']} panel {r['printed_panel'] or '-'} "
             f"({r['source_dir'] or 'no source_dir'}): "
             f"reproduces_published is 'unknown'"
             + (f" with verdict {r['verdict']!r}" if r["verdict"] else "")
             + " - decide whether it reproduces, does not, or has nothing to "
               "compare against, and record that")


def _path_constants():
    """The Path constants panel scripts import from 00_Config/paths.py."""
    sys.path.insert(0, str(ROOT / "00_Config"))
    import paths
    return {k: v for k, v in vars(paths).items()
            if k.isupper() and isinstance(v, Path)}


class _Resolver:
    """Enough of pathlib to follow how a panel script names a file.

    Handles CONST / "name.csv", Path(__file__).parent, .resolve(),
    .parents[n], .joinpath(), os.path.join(), an f-string whose every part
    resolves, and module-level names bound to any of those.
    Anything else resolves to None and is reported as unresolved rather than
    guessed at.

    The f-string and os.path.join forms were added on 2026-09-03 for the glob
    patterns of SUPPLEMENTARY_AUDIT.md fault 3: glob.glob(f"{DIR}/*.csv") and
    glob.glob(os.path.join(str(DIR), "*.csv")) named a directory that check 9
    could not see at all - not as an input, not as unresolved, not anywhere.
    An f-string resolves only if EVERY part of it resolves; one run-time piece
    makes the whole thing None, so nothing is guessed at.
    """

    def __init__(self, script, consts):
        self.env = dict(consts)
        self.env["__file__"] = script

    def resolve(self, node):
        if isinstance(node, ast.Name):
            return self.env.get(node.id)
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return Path(node.value)
        if isinstance(node, ast.Attribute):
            base = self.resolve(node.value)
            return base.parent if base is not None and node.attr == "parent" else None
        if isinstance(node, ast.Subscript):
            f = node.value
            if isinstance(f, ast.Attribute) and f.attr == "parents":
                base, idx = self.resolve(f.value), node.slice
                if base is not None and isinstance(idx, ast.Constant) \
                        and isinstance(idx.value, int):
                    try:
                        return base.parents[idx.value]
                    except IndexError:
                        return None
            return None
        if isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Name) and f.id in ("Path", "str") and node.args:
                return self.resolve(node.args[0])
            if isinstance(f, ast.Attribute) and f.attr in ("resolve", "absolute"):
                return self.resolve(f.value)
            if isinstance(f, ast.Attribute) and f.attr == "joinpath":
                base = self.resolve(f.value)
                for a in node.args:
                    part = self.resolve(a)
                    if base is None or part is None:
                        return None
                    base = base / part
                return base
            # os.path.join(...). Guarded on the receiver so that "".join() and
            # every other str.join in the tree stays unresolved.
            if isinstance(f, ast.Attribute) and f.attr == "join" \
                    and isinstance(f.value, ast.Attribute) \
                    and f.value.attr == "path" \
                    and isinstance(f.value.value, ast.Name) \
                    and f.value.value.id == "os":
                base = None
                for a in node.args:
                    part = self.resolve(a)
                    if part is None:
                        return None
                    base = part if base is None else base / part
                return base
            return None
        if isinstance(node, ast.JoinedStr):
            parts = []
            for v in node.values:
                if isinstance(v, ast.Constant) and isinstance(v.value, str):
                    parts.append(v.value)
                    continue
                inner = v.value if isinstance(v, ast.FormattedValue) else v
                # A conversion or a format spec would change the text; refuse
                # rather than resolve to something the script does not use.
                if isinstance(v, ast.FormattedValue) and (
                        v.conversion not in (-1, None) or v.format_spec is not None):
                    return None
                got = self.resolve(inner)
                if got is None:
                    return None
                parts.append(str(got))
            return Path("".join(parts)) if parts else None
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
            left, right = self.resolve(node.left), self.resolve(node.right)
            return None if left is None or right is None else left / right
        return None


def _io_contexts(tree):
    """How each path expression in a script is used: read, written, or neither.

    Returns ({id(node): {"read"|"write"}}, {variable name: {...}}). The second
    map is what lets `TABLE = HERE / "prepared.csv"` at the top of a file be
    judged by the `pd.read_csv(TABLE)` two hundred lines below it. Location
    cannot make that call; use can.
    """
    by_node, by_name = {}, {}

    def mark(node, kind):
        if node is None:
            return
        by_node.setdefault(id(node), set()).add(kind)
        if isinstance(node, ast.Name):
            by_name.setdefault(node.id, set()).add(kind)

    for n in ast.walk(tree):
        if not isinstance(n, ast.Call):
            continue
        f = n.func
        attr = f.attr if isinstance(f, ast.Attribute) else None
        name = attr or getattr(f, "id", "")
        arg = n.args[0] if n.args else None
        if arg is None:
            for kw in n.keywords:
                if kw.arg in ("filepath_or_buffer", "path", "path_or_buf",
                              "filename", "fname", "file", "io", "filepath"):
                    arg = kw.value
                    break
        if name == "open":
            mode = ""
            if len(n.args) > 1 and isinstance(n.args[1], ast.Constant):
                mode = str(n.args[1].value)
            for kw in n.keywords:
                if kw.arg == "mode" and isinstance(kw.value, ast.Constant):
                    mode = str(kw.value.value)
            written = any(c in mode for c in "wax")
            # Path.open() names its receiver; open(path), Image.open(path) and
            # fitz.open(path) name their argument. Tell them apart by the first
            # positional: a mode string, or nothing, means the receiver is the
            # path.
            receiver = attr == "open" and (
                arg is None or (isinstance(arg, ast.Constant)
                                and isinstance(arg.value, str)
                                and len(arg.value) <= 4
                                and set(arg.value) <= set("rwaxbt+U")))
            mark(f.value if receiver else arg, WRITE if written else READ)
        elif attr in READ_METHODS:
            mark(f.value, READ)
        elif attr in WRITE_METHODS:
            mark(f.value, WRITE)
        elif name in READ_CALLS:
            mark(arg, READ)
        elif name in WRITE_CALLS:
            mark(arg, WRITE)
    return by_node, by_name


def _dir_read(node, r):
    """(directory, pattern) for a call that reads a whole directory, else None.

    `directory` is None when the base cannot be resolved statically and
    `pattern` is None when it is built at run time; either way the call is
    still returned, so that check_inputs can report it as unresolved. A
    directory read that resolves to nothing must be visible - the shape this
    check could not see is a whole archived run read one file at a time.
    """
    f = node.func
    attr = f.attr if isinstance(f, ast.Attribute) else None
    name = attr or getattr(f, "id", "")

    def const(a):
        return (a.value if isinstance(a, ast.Constant)
                and isinstance(a.value, str) else None)

    # glob.glob("/abs/dir/*.csv") and glob.iglob - directory and pattern in one.
    # The argument is not always a literal: glob.glob(f"{DIR}/*.csv") and
    # glob.glob(os.path.join(str(DIR), "*.csv")) name a directory just as
    # plainly, and until 2026-09-03 both resolved to nothing at all. _Resolver
    # is asked second, so a literal still behaves exactly as before.
    if attr in ("glob", "iglob") and isinstance(f.value, ast.Name) \
            and f.value.id == "glob":
        s = const(node.args[0]) if node.args else None
        if s is None and node.args:
            got = r.resolve(node.args[0])
            s = str(got) if got is not None else None
        if s is None:
            return None, None
        p = Path(s)
        return (p.parent if p.is_absolute() else None), p.name
    if name in ("listdir", "scandir"):
        return (r.resolve(node.args[0]) if node.args else None), "*"
    if attr in DIR_METHODS:
        base = r.resolve(f.value)
        if attr == "iterdir":
            return base, "*"
        pat = const(node.args[0]) if node.args else None
        if pat is None and node.args:
            # Path.glob(f"...{const}...") - a pattern built out of pieces that
            # all resolve is still a pattern this check can test.
            got = r.resolve(node.args[0])
            pat = str(got) if got is not None else None
        if pat is None:
            return base, None
        return base, ("**/" + pat if attr == "rglob" else pat)
    return None


def _longest_path(node, r):
    """The deepest absolute path anywhere inside an expression, or None.

    A read path built at run time - DIR / name, an f-string, a .format() - can
    still name the directory it comes out of, and that is the part worth
    checking against the archive markers.
    """
    best = None
    for n in ast.walk(node):
        p = r.resolve(n)
        if p is not None and p.is_absolute() \
                and (best is None or len(p.parts) > len(best.parts)):
            best = p
    return best


def _describe(node, limit=90):
    try:
        text = " ".join(ast.unparse(node).split())
    except Exception:
        return type(node).__name__
    return text if len(text) <= limit else text[:limit - 3] + "..."


def _inside(p, own):
    try:
        p.resolve().relative_to(own)
        return True
    except ValueError:
        return False


def _inputs_of(script, consts):
    """What a panel script reads.

    Returns (files, dirs, unresolved, error).

      files       {absolute path: (line, why)} - named data files it reads
      dirs        [(line, directory|None, pattern|None)] - directory reads
      unresolved  [(line, description, nearest resolvable path|None)]

    A data-suffixed path resolving outside the script's own directory is an
    input, as it always was. One resolving inside it is no longer assumed to be
    an output: it is judged by use. Read makes it an input; write, or write
    then read back, makes it an output; neither leaves it alone.
    """
    r = _Resolver(script, consts)
    try:
        tree = ast.parse(script.read_text(errors="replace"))
    except SyntaxError as e:
        return None, None, None, f"{e}"
    by_node, by_name = _io_contexts(tree)
    own = script.parent.resolve()
    files, dirs, unresolved = {}, [], []
    seen_dirs, seen_unres = set(), set()

    # Names bound by a loop over a directory read are the files that read
    # already accounts for. Reporting each of them again would bury it.
    loop_bound = set()
    for n in ast.walk(tree):
        it = getattr(n, "iter", None)
        if it is None or not isinstance(n, (ast.For, ast.AsyncFor,
                                            ast.comprehension)):
            continue
        if any(isinstance(c, ast.Call) and _dir_read(c, r) is not None
               for c in ast.walk(it)):
            for t in ast.walk(n.target):
                if isinstance(t, ast.Name):
                    loop_bound.add(t.id)

    def context(node):
        c = set(by_node.get(id(node), ()))
        if isinstance(node, ast.Name):
            c |= by_name.get(node.id, set())
        return c

    # {name: the expression it was last bound to} for bindings that do not
    # resolve, so that a read of a bare `f` can be reported as the path
    # expression behind it rather than as the letter f.
    pending = {}

    # Everything inside a call that has already been accounted for as a
    # directory read. glob.glob(str(DIR / "*.tsv")) would otherwise be reported
    # twice: once as the directory, and once as a "file" called *.tsv that does
    # not exist. The directory report is the true one.
    globbed = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and _dir_read(node, r) is not None:
            for sub in ast.walk(node):
                if sub is not node:
                    globbed.add(id(sub))

    for node in ast.walk(tree):
        if id(node) in globbed:
            continue
        if isinstance(node, ast.Assign):
            val = r.resolve(node.value)
            for t in node.targets:
                if not isinstance(t, ast.Name):
                    continue
                if val is not None:
                    r.env[t.id] = val
                elif isinstance(node.value, (ast.BinOp, ast.Call, ast.Subscript,
                                             ast.JoinedStr, ast.Attribute)):
                    pending[t.id] = node.value
        if not isinstance(node, (ast.Name, ast.BinOp, ast.Call, ast.Attribute,
                                 ast.Subscript, ast.JoinedStr)):
            continue
        line = getattr(node, "lineno", 0)

        if isinstance(node, ast.Call):
            d = _dir_read(node, r)
            if d is not None:
                directory, pattern = d
                if pattern and pattern.endswith(ARTEFACT_PATTERNS):
                    continue          # figure assembly, not a data input
                key = (str(directory), pattern)
                if key not in seen_dirs:
                    seen_dirs.add(key)
                    dirs.append((line, directory, pattern))
                continue

        p = r.resolve(node)
        if p is not None:
            if p.is_absolute() and p.suffix.lower() in DATA_SUFFIX:
                if not _inside(p, own):
                    files[str(p)] = (line, "named")
                else:
                    c = context(node)
                    if READ in c and WRITE not in c:
                        files[str(p)] = (line, "read from beside the script")
            continue

        # Unresolvable, and used as a read: name it rather than drop it.
        if READ not in context(node):
            continue
        open_names = {n.id for n in ast.walk(node)
                      if isinstance(n, ast.Name) and r.env.get(n.id) is None}
        if open_names and open_names <= loop_bound:
            continue
        shown_node = node
        if isinstance(node, ast.Name) and node.id in pending:
            shown_node = pending[node.id]
        near = _longest_path(shown_node, r)
        desc = _describe(shown_node)
        key = (desc, str(near))
        if key not in seen_unres:
            seen_unres.add(key)
            unresolved.append((line, desc, near))
    return files, dirs, unresolved, None


def check_inputs(rows):
    """Where every panel script actually reads from."""
    consts = _path_constants()
    scripts, seen = [], set()
    for r in rows:
        for rel in filter(None, r["source_script"].split("; ")):
            p = resolve(r["figure"], rel)
            if p and p.suffix == ".py" and p.exists() and p not in seen:
                seen.add(p)
                scripts.append(p)
    checked = stale = open_ended = dirs_checked = 0
    for s in scripts:
        files, dirs, unres, err = _inputs_of(s, consts)
        if err:
            fail(f"{s.name} does not parse: {err}")
            continue
        if not files and not dirs and not unres:
            if s.name not in NO_DATA_INPUTS:
                open_ended += 1
                note(f"{s.name}: no input path could be resolved - not checked")
            continue
        for path, (line, why) in sorted(files.items()):
            checked += 1
            if not Path(path).exists():
                fail(f"{s.name}:{line} reads {path}, which does not exist"
                     + (" - it is read from beside the script, not written "
                        "by it" if why != "named" else ""))
            elif any(m in path for m in STALE_MARKERS):
                stale += 1
                fail(f"{s.name}:{line} reads a live input out of an archive, "
                     f"backup or temp workspace: {path} - move it to a live "
                     f"location or say in PROVENANCE.csv why it belongs there")
        for line, directory, pattern in dirs:
            if directory is None:
                open_ended += 1
                note(f"{s.name}:{line} reads a directory it names at run time "
                     f"- UNRESOLVED, nothing under it is checked")
                continue
            shown = pattern or "<pattern built at run time>"
            if any(m in str(directory) for m in STALE_MARKERS):
                stale += 1
                fail(f"{s.name}:{line} reads a live input out of an archive, "
                     f"backup or temp workspace: {directory}/{shown} - move it "
                     f"to a live location or say in PROVENANCE.csv why it "
                     f"belongs there")
                continue
            if not directory.is_dir():
                fail(f"{s.name}:{line} reads {directory}/{shown}, whose "
                     f"directory does not exist")
                continue
            dirs_checked += 1
            if pattern is None:
                open_ended += 1
                note(f"{s.name}:{line} globs {directory} with a pattern built "
                     f"at run time - the directory is live, the files it "
                     f"matches are UNRESOLVED")
                continue
            if not sorted(directory.glob(pattern)):
                fail(f"{s.name}:{line} reads {directory}/{pattern}, which "
                     f"matches no file")
        for line, desc, near in unres:
            if near is not None and any(m in str(near) for m in STALE_MARKERS):
                stale += 1
                fail(f"{s.name}:{line} builds a read path under an archive, "
                     f"backup or temp workspace: {near} - {desc}")
                continue
            # SUPPLEMENTARY_AUDIT.md fault 3, row 5: PREPARATION / "GSEA" /
            # phase / f"{name}_hallmark.csv" names a real directory even though
            # the leaf is built at run time. The leaf cannot be checked; the
            # directory can, and until 2026-09-03 its absence was written into
            # a note that nothing reads. A read whose own resolvable root is
            # gone is the Figure 5A fault, not an open question.
            if near is not None and not near.exists():
                fail(f"{s.name}:{line} builds a read path under {near}, which "
                     f"does not exist - {desc}")
                continue
            open_ended += 1
            where = f" (under {near})" if near is not None else ""
            note(f"{s.name}:{line} reads a path it builds at run time - "
                 f"UNRESOLVED: {desc}{where}")
    if not stale:
        note(f"{checked} panel inputs and {dirs_checked} input directories "
             f"across {len(scripts)} scripts exist and none is read out of an "
             f"archive, backup or temp workspace"
             + (f" ({open_ended} unresolved)" if open_ended else ""))


def _stale_marker(text):
    """The archive marker a path carries, or None."""
    for m in STALE_MARKERS:
        if m in str(text):
            return m
    return None


def check_upstream_inputs():
    """Check 11: the archive guard, one step upstream of the panels.

    Check 9 asks where a *panel script* reads from. Everything a panel script
    reads was written by something else, and that something else has inputs of
    its own. A prior model, a reference matrix or a deconvolution table
    repointed into an archive up there reaches the printed page just as surely,
    and nothing in this file looked at it until 2026-09-03.

    Reuses check 9's machinery rather than growing a second one: _inputs_of and
    _Resolver read the Python, STALE_MARKERS decides what an archive is.
    """
    consts = _path_constants()
    prep = consts.get("PREPARATION")
    py_roots = [("02_Preparation_for_Panels", UPSTREAM)]
    text_roots = [("02_Preparation_for_Panels", UPSTREAM),
                  ("Round_5/02_Preparation_for_Panels", prep)]

    stale = archived = 0
    n_py = n_text = 0

    for label, root in py_roots:
        if root is None or not root.is_dir():
            fail(f"check 11 cannot reach {label} - no upstream script was read, "
                 f"and a check that reads nothing reports nothing")
            continue
        for p in sorted(root.rglob("*.py")):
            if "__pycache__" in p.parts:
                continue
            if _stale_marker(p):
                archived += 1
                continue
            files, dirs, unres, err = _inputs_of(p, consts)
            if err:
                fail(f"upstream {p.relative_to(root)} does not parse: {err}")
                continue
            n_py += 1
            hits = [(line, path) for path, (line, _) in files.items()]
            hits += [(line, f"{d}/{pat or '<pattern built at run time>'}")
                     for line, d, pat in dirs if d is not None]
            hits += [(line, near) for line, _, near in unres if near is not None]
            for line, where in hits:
                marker = _stale_marker(where)
                if marker:
                    stale += 1
                    fail(f"upstream {p.relative_to(root)}:{line} reads a live "
                         f"input out of an archive, backup or temp workspace "
                         f"({marker}): {where} - move it to a live location or "
                         f"record in PROVENANCE.csv why it belongs there")

    for label, root in text_roots:
        if root is None or not root.is_dir():
            fail(f"check 11 cannot reach {label} - no upstream config was read, "
                 f"and a check that reads nothing reports nothing")
            continue
        for pat in TEXT_INPUT_PATTERNS:
            for p in sorted(root.rglob(pat)):
                if "__pycache__" in p.parts:
                    continue
                if _stale_marker(p):
                    archived += 1
                    continue
                n_text += 1
                text = p.read_text(errors="replace")
                for s in sorted(set(ABS_PATH.findall(text))):
                    marker = _stale_marker(s)
                    if not marker:
                        continue
                    line = text[:text.index(s)].count("\n") + 1
                    stale += 1
                    fail(f"upstream {label}/{p.relative_to(root)}:{line} names "
                         f"an input inside an archive, backup or temp workspace "
                         f"({marker}): {s} - move it to a live location or "
                         f"record in PROVENANCE.csv why it belongs there")

    if not n_py or not n_text:
        fail(f"check 11 read {n_py} upstream scripts and {n_text} upstream "
             f"configs - it cannot have found anything, so its silence means "
             f"nothing")
        return
    note(f"{n_py} upstream scripts and {n_text} upstream configs scanned for "
         f"archived inputs ({archived} skipped as themselves archived, "
         f"{stale} live reads out of an archive)")


def check_broken_expires(rows):
    """A "no" must name a mechanism and come up for re-test."""
    today = datetime.date.today()
    for r in rows:
        if r["reproduces_published"] != "no":
            continue
        tag = f"Figure {r['figure']} panel {r['printed_panel']}"
        if len(r["note"]) < 80:
            fail(f"{tag} is marked reproduces_published=no with a note too short "
                 f"to name a mechanism - say which input, which run or which step "
                 f"differs, not merely that the output does")
        d = panel_dir(r)
        marker = d / BROKEN_MARKER if d else None
        if marker is None or not marker.exists():
            continue                       # check 3 already failed this row
        m = RETEST.search(marker.read_text(errors="replace"))
        if not m:
            fail(f"{tag}: {BROKEN_MARKER} carries no 'Retest-by: YYYY-MM-DD' "
                 f"line - a marker that never expires lets a wrong finding "
                 f"harden into a fact, which is what happened here")
            continue
        due = datetime.date(*(int(x) for x in m.groups()))
        if due < today:
            fail(f"{tag}: {BROKEN_MARKER} came up for re-test on {due} and has "
                 f"not been revisited - re-run the panel against the printed "
                 f"figure and either renew the date or lift the marker")
        else:
            note(f"{tag}: recorded not reproducible, re-test due {due}")


def check_guard_rails():
    for p in (ROOT / "CLAUDE.md", GROUND / "README.md",
              HERE / "Main_Figures" / "README.md",
              HERE / "Supplementary_New" / "README.md"):
        if not p.exists():
            fail(f"missing guard-rail document: {p}")
    figs = sorted((GROUND / "figures").glob("*.pdf"))
    if not figs:
        fail(f"{GROUND / 'figures'} holds no figures")
    writable = [f.name for f in figs if f.stat().st_mode & 0o222]
    if writable:
        fail(f"ground-truth figures are writable: {', '.join(writable)}")


def check_freshness(rows):
    """A verdict older than the script it judges is not evidence about it.

    Nothing else in this file would notice. The verdict, its date and its
    evidence all stay in place when the script underneath them changes, and the
    row goes on asserting something that was measured against a file that has
    since been edited.
    """
    import datetime
    checked = 0
    for r in rows:
        if not r["verdict"] or not r["verdict_date"]:
            continue
        d = panel_dir(r)
        script = r["source_script"].split(";")[0].strip()
        if not d or not script:
            continue
        path = d.parent / script if (d.parent / script).exists() else d / Path(script).name
        if not path.exists():
            continue
        # Figures 1 and 4 resolve into the generated release, where every file
        # carries the mtime of the last update_release.py run rather than of an
        # edit. Judging those would fail the whole check every time the release
        # is rebuilt. The Round_5 original they are generated from is the file
        # the verdict was actually made against.
        if RELEASE in path.parents:
            source = (ROOT.parent / "Round_5" / "03_Final_Panels"
                      / path.relative_to(RELEASE))
            if not source.exists():
                continue
            path = source
        try:
            judged = datetime.date.fromisoformat(r["verdict_date"])
        except ValueError:
            fail(f"Figure {r['figure']} panel {r['printed_panel']}: "
                 f"verdict_date {r['verdict_date']!r} is not a date")
            continue
        edited = datetime.date.fromtimestamp(path.stat().st_mtime)
        checked += 1
        if edited > judged:
            fail(f"Figure {r['figure']} panel {r['printed_panel']}: "
                 f"{path.name} was edited on {edited} but its verdict is dated "
                 f"{judged} - re-run the panel and re-judge it, or the verdict "
                 f"describes a file that no longer exists")
    note(f"{checked} verdicts are newer than the scripts they judge")


def main():
    rows = load()
    check_structure(rows)
    check_coverage(rows)
    check_known_broken(rows)
    check_placement(rows)
    check_guard_rails()
    check_freshness(rows)
    check_lookup(rows)
    check_no_unknowns(rows)
    check_inputs(rows)
    check_upstream_inputs()
    check_broken_expires(rows)

    print(f"PROVENANCE.csv: {len(rows)} rows")
    for tag in ("yes", "no", "unknown", "na"):
        n = sum(1 for r in rows if r["reproduces_published"] == tag)
        print(f"  reproduces_published={tag:<8} {n:>3}")
    for n in notes:
        print(f"  ok    {n}")
    if failures:
        print(f"\n{len(failures)} FAILURE(S):")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print("\nAll provenance checks pass.")


if __name__ == "__main__":
    main()
