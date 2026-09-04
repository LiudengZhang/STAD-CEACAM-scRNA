#!/usr/bin/env python3
"""
Shared support for the restyle drivers.

Twenty-one of the twenty-six S7-S11 panels are drawn by fourteen analysis
scripts under `02_New_Analyses/*/scripts/`, which in the same run compute the
statistics and write the CSVs `verify_numbers.py` checks. A driver restyles
those panels without the statistics being touched, duplicated or re-run:

  1. it **imports** the analysis module, so every constant, colour map, label
     helper and lineage table comes from the one file that defines it, and
     `main()` never runs (the analysis scripts all guard it with
     `if __name__ == "__main__"`);
  2. it **reads the CSVs that `main()` already wrote**, so nothing recomputes,
     no h5ad is opened, and no GSEA table is read;
  3. it holds a **restyled copy of the drawing code only** - the part that
     turns an already-computed table into marks on paper;
  4. it **proves** that copy draws the same numbers as the original, by calling
     the analysis module's own `_panel_*` function on the same frames and
     diffing every plotted value with `compare_panel_content`.

Point 2 is what keeps the S10E input question open rather than silently
answering it. `adaptive_immune_resource.py:146` reads the damaged-matrix GSEA
tables and line 147 skips a missing file in silence; a driver that read the
CSV never reaches either line, so re-running the restyle neither acts on the
author's pending ruling nor swallows a failure. The same holds for
`design_defence.py`'s bare `except Exception` around the SCENIC read.

Point 4 is the gate. `check` fails loudly if one value moves.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[4]          # Round_7_major_revision
RESTYLED = Path(__file__).resolve().parents[1]      # Supplementary_New/_restyled
NEW_ANALYSES = ROOT / "02_New_Analyses"

sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(ROOT / "10_Reproduction"))

import panel_style_cns as style                      # noqa: E402
import compare_panel_content as cpc                  # noqa: E402


def analysis(rel):
    """Import an analysis module by path without running its main().

    `rel` is relative to 02_New_Analyses, e.g.
    "09_R2.2_Adaptive_Immune/scripts/adaptive_immune_resource.py".
    """
    path = NEW_ANALYSES / rel
    if not path.exists():
        raise FileNotFoundError(path)
    name = "analysis_" + path.stem
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)                     # module level only
    return mod


def outputs_of(mod):
    """The `outputs/` directory the analysis module writes its tables to."""
    return Path(mod.OUT)


def panel_dir(figure, panel):
    """Where a restyled panel is written. Mirrors Version A's layout."""
    d = RESTYLED / figure / panel
    d.mkdir(parents=True, exist_ok=True)
    return d


def require(path, what):
    """A table the driver reads must exist. No silent skip, ever.

    The audit found two silent-degradation paths in these scripts - a bare
    `except Exception` and an `if not f.exists(): continue` - either of which
    turns a missing input into a quietly incomplete panel that still exits 0.
    A driver must not add a third, so every read goes through here.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"{what} is missing: {p}\n"
            f"    The driver reads the tables the analysis already wrote and "
            f"never recomputes them. Run the analysis script itself if this "
            f"file needs to be produced - but read RESTYLE_REPORT.md first: "
            f"two of these scripts carry faults awaiting the author's ruling.")
    if p.stat().st_size == 0:
        raise ValueError(f"{what} is empty: {p}")
    return p


def check(label, original, restyled, nd=9):
    """Prove the restyled drawing plots what the original plots.

    `original` and `restyled` are zero-argument callables that each build a
    figure. Neither is allowed to write to disk; `cpc.capture` enforces that.
    Returns 0 when the content matches.
    """
    plt.close("all")
    with cpc.capture() as cap_old:
        original()
    old = list(cap_old.figures) or [plt.figure(n) for n in plt.get_fignums()]
    old_content = cpc.harvest_figures(old, nd)

    plt.close("all")
    with cpc.capture() as cap_new:
        restyled()
    new = list(cap_new.figures) or [plt.figure(n) for n in plt.get_fignums()]
    new_content = cpc.harvest_figures(new, nd)

    plt.close("all")
    return cpc.report(cpc.diff(old_content, new_content, "figures"), label)


def save(fig, figure, panel, stem):
    """Write a restyled panel into its directory, at exactly its own size."""
    out = panel_dir(figure, panel) / stem
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING {stem}: ink outside the canvas (l,r,b,t mm) "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, out)
    w, h = fig.get_size_inches()
    print(f"  {stem:<44} {w * 25.4:6.1f} x {h * 25.4:5.1f} mm")
    return out


def run(panels, argv=None):
    """Driver entry point. `panels` maps label -> (build, original) callables.

    build()    draws and saves the restyled panel, returning nothing
    original() draws the analysis module's own version, for the check
    """
    argv = sys.argv[1:] if argv is None else argv
    only = [a for a in argv if not a.startswith("-")]
    do_check = "--check" in argv
    rc = 0
    for label, (build, original) in panels.items():
        if only and label not in only:
            continue
        if do_check:
            rc |= check(label, original, build)
        else:
            build()
    return rc
