#!/usr/bin/env python3
"""
Shared support for the supplementary panel drivers.

Most of the supplementary panels are drawn by analysis scripts under
`04_Revision_Analyses/*/scripts/`, which in the same run compute the statistics and
write the CSVs `verify_numbers.py` checks. A driver draws one of those panels
at print size without the statistics being touched, duplicated or re-run:

  1. it **imports** the analysis module, so every constant, colour map, label
     helper and lineage table comes from the one file that defines it, and
     `main()` never runs (the analysis scripts all guard it with
     `if __name__ == "__main__"`);
  2. it **reads the CSVs that `main()` already wrote**, so nothing recomputes,
     no h5ad is opened, and no GSEA table is read;
  3. it holds a **copy of the drawing code only** - the part that
     turns an already-computed table into marks on paper;
  4. it **proves** that copy draws the same numbers as the original, by calling
     the analysis module's own `_panel_*` function on the same frames and
     diffing every plotted value with `compare_panel_content`.

Point 2 keeps an open input question open rather than silently answering it: a
driver that reads a written table never reaches the analysis's own file
discovery, so re-running a driver neither acts on a pending ruling nor swallows
a missing input.

Point 4 is the gate. `check` fails loudly if one value moves, and it runs both
sides with every write diverted, so a comparison cannot rewrite an output.
"""

from __future__ import annotations

import contextlib
import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SUPP = HERE.parent                                  # Supplementary_New
ROOT = HERE.parents[2]                              # project root
NEW_ANALYSES = ROOT / "04_Revision_Analyses"

sys.path.insert(0, str(ROOT / "00_Config"))
sys.path.insert(0, str(ROOT / "10_Reproduction"))

import panel_style_cns as style                      # noqa: E402


def _cpc():
    """`compare_panel_content`, imported only when the gate is actually run.

    It used to be imported at module scope, which made a VERIFICATION TOOL a
    DRAWING DEPENDENCY. The module lives in 10_Reproduction/, which the code
    release does not ship - it is not part of the deposited pipeline - so in
    the capsule every one of the eight drivers died at
    `import compare_panel_content` before drawing a single panel, and
    Supplementary Figures S7-S9 could not be rebuilt at all. Measured in the
    reviewer simulation of 2026-09-10: seven drivers with exactly that
    ModuleNotFoundError, and assemble_new_supplementaries.py failing behind
    them for want of a panel SVG.

    Author's ruling 2026-09-10: make the import lazy. Drawing must not need it;
    `check()` may, and says so plainly when it is missing rather than failing
    somewhere further down.
    """
    try:
        import compare_panel_content as cpc
    except ModuleNotFoundError as e:                      # noqa: PERF203
        raise SystemExit(
            f"the content gate needs 10_Reproduction/compare_panel_content.py, "
            f"which is not on the path ({e}).\n"
            f"    It is a verification tool and is deliberately not part of the "
            f"code release, so --check cannot be run from the deposited "
            f"capsule. Run it in the working tree, where 10_Reproduction/ "
            f"sits beside 03_Final_Panels/. Drawing the panels needs none of "
            f"this and works either way.") from e
    return cpc

MM_TO_INCH = style.MM_TO_INCH

#: The type specification every panel in this set is set at: axis labels and
#: titles at 7 pt, tick labels, legends and annotations at 6 pt.
TYPE_SPEC = dict(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)

#: Paper left between the drawn ink and the edge of the panel canvas.
PAD_MM = 0.6

#: Where the authorised label changes are declared. A string a panel prints in
#: a shortened form so that it fits the slot it is drawn in has to be named
#: here; the gate then reports it by name, and an undeclared string change
#: still fails.
LABELS = ROOT / "00_Config" / "shared" / "labels.py"


def aliases():
    """The authorised label changes, loaded from the declaration file.

    Only the gate needs these. `run()` resolves them lazily, so a driver hands
    it THIS FUNCTION rather than its result - calling it while the panels dict
    is built would put the gate back on the drawing path that _cpc() exists to
    keep it off.
    """
    return _cpc().load_aliases(LABELS)


def apply_style():
    """Put the type specification into matplotlib's rcParams."""
    return style.apply(**TYPE_SPEC)


@contextlib.contextmanager
def _no_mkdir():
    """Let nothing under this block create a directory.

    Two things need it. `capture` diverts every write, but a panel function's
    first statement is usually `d.mkdir(parents=True, exist_ok=True)` for the
    directory it would have saved into; left alone, comparing against a panel
    whose tree has been retired quietly recreates that tree as empty
    directories, and the provenance coverage check then asks for rows to
    describe them. And one analysis module creates its own panel directory at
    import time, so merely importing it for its constants put the directory
    back.
    """
    real = Path.mkdir
    Path.mkdir = lambda self, *a, **k: None
    try:
        yield
    finally:
        Path.mkdir = real


def analysis(rel):
    """Import an analysis module by path without running its main().

    `rel` is relative to 04_Revision_Analyses, e.g.
    "09_R2.2_Adaptive_Immune/scripts/adaptive_immune_resource.py".
    """
    path = NEW_ANALYSES / rel
    if not path.exists():
        raise FileNotFoundError(path)
    name = "analysis_" + path.stem
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    with _no_mkdir():
        spec.loader.exec_module(mod)                 # module level only
    return mod


def outputs_of(mod):
    """The `outputs/` directory the analysis module writes its tables to."""
    return Path(mod.OUT)


def panel_dir(figure, panel, root=None):
    """Where a panel is written. `figure` may be None for a loose panel."""
    d = (root or SUPP)
    if figure:
        d = d / figure
    d = d / panel
    d.mkdir(parents=True, exist_ok=True)
    return d


def require(path, what):
    """A table the driver reads must exist. No silent skip, ever.

    Two silent-degradation paths exist in the analysis scripts - a bare
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
            f"file needs to be produced.")
    if p.stat().st_size == 0:
        raise ValueError(f"{what} is empty: {p}")
    return p


def _clear_figure_title(fig, gap_mm=0.8):
    """Lower the axes box until a figure-level title sits clear of the axes.

    `fit_margins` fits the union of the ink into the canvas; it has no notion
    of two pieces of ink overlapping each other. A figure-level title is pinned
    near the top of the canvas, so fitting the axes to that same edge prints
    the axis titles underneath it and nothing reports an overflow, because
    neither has left the page. The overlap is measured here and the axes box is
    lowered by exactly that much, which can only move ink inward.
    """
    tops = [t for t in fig.texts
            if t.get_visible() and t.get_position()[1] > 0.5]
    if not tops or not fig.axes:
        return
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi, (w_in, h_in) = fig.dpi, fig.get_size_inches()
    floor = min(t.get_window_extent(renderer=r).y0 for t in tops)
    ceiling = -float("inf")
    for ax in fig.axes:
        boxes = [ax.get_window_extent(renderer=r)]
        if ax.title.get_text():
            boxes.append(ax.title.get_window_extent(renderer=r))
        ceiling = max([ceiling] + [b.y1 for b in boxes])
    drop = ceiling - (floor - gap_mm * MM_TO_INCH * dpi)
    if drop <= 0:
        return
    sp = fig.subplotpars
    fig.subplots_adjust(top=(sp.top * h_in - drop / dpi) / h_in)


def fit(fig, **subplot_kw):
    """Set the margins from the rendered ink, and prove the panel is clean.

    A margin typed as a number is a guess about how wide a tick label will be,
    and at 1:1 a guess that is 2 mm short does not shrink away - the label runs
    off the canvas and the assembler's nested viewport clips it silently. So
    the margins are measured, and both gates are hard failures: no ink outside
    the canvas, and nothing under the panel letter.
    """
    if subplot_kw:
        fig.subplots_adjust(**subplot_kw)
    style.fit_margins(fig, pad_mm=PAD_MM)
    _clear_figure_title(fig)
    over = style.overflow_mm(fig)
    if max(over) > 0.0:
        raise RuntimeError(
            f"ink outside the canvas (left, right, bottom, top mm): {over}")
    intruders = style.letter_clear(fig)
    if intruders:
        raise RuntimeError(f"the panel-letter cell is not clear: {intruders}")
    return over


def check(label, original, redrawn, nd=9, alias=None, keep_axes=None):
    """Prove the redrawn panel plots what the original plots.

    `original` and `redrawn` are zero-argument callables that each build a
    figure. Every write either side makes is diverted, so a comparison cannot
    reach the CSVs `verify_numbers.py` checks. `alias`, when given, is the
    declared set of authorised label changes; each one that fires is printed by
    name and an undeclared string change still fails. Returns 0 when the
    content matches.
    """
    cpc = _cpc()
    plt.close("all")
    with _no_mkdir(), cpc.capture(writes="divert") as cap_old:
        original()
    old = list(cap_old.figures) or [plt.figure(n) for n in plt.get_fignums()]
    old_content = cpc.harvest_figures(old, nd)
    if keep_axes is not None:
        old_content = [[fig[i] for i in keep_axes] for fig in old_content]

    plt.close("all")
    with _no_mkdir(), cpc.capture(writes="divert") as cap_new:
        redrawn()
    new = list(cap_new.figures) or [plt.figure(n) for n in plt.get_fignums()]
    new_content = cpc.harvest_figures(new, nd)

    plt.close("all")
    diverted = len(cap_old.diverted) + len(cap_new.diverted)
    if diverted:
        print(f"  ({diverted} write(s) diverted, nothing reached disk)")
    rc = cpc.report(cpc.diff(old_content, new_content, "figures", alias), label)
    if alias is not None:
        alias.report()
    return rc


def save(fig, figure, panel, stem, root=None):
    """Write a panel into its directory, at exactly the size it prints at."""
    out = panel_dir(figure, panel, root) / stem
    over = style.overflow_mm(fig)
    if max(over) > 0.0:
        raise RuntimeError(
            f"{stem}: ink outside the canvas (left, right, bottom, top mm) "
            f"{over}. A panel drawn 1:1 has no assembler shrink left to hide "
            f"it, so this is a failure and not a warning.")
    intruders = style.letter_clear(fig)
    if intruders:
        raise RuntimeError(f"{stem}: the panel-letter cell is not clear: "
                           f"{intruders}")
    style.save_panel(fig, out)
    w, h = fig.get_size_inches()
    print(f"  {stem:<44} {w * 25.4:6.1f} x {h * 25.4:5.1f} mm  "
          f"overflow {over}  letter cell clear")
    return out


def run(panels, argv=None):
    """Driver entry point. `panels` maps label -> (build, original) callables.

    build()    draws and saves the panel, returning nothing
    original() draws the analysis module's own version, for the check

    A third element, when present, is the declared label changes to allow; a
    fourth is the axis indices of the original that the panel retains.
    """
    argv = sys.argv[1:] if argv is None else argv
    only = [a for a in argv if not a.startswith("-")]
    do_check = "--check" in argv
    rc = 0
    for label, spec in panels.items():
        build, original = spec[0], spec[1]
        # A driver may hand a CALLABLE here instead of a resolved alias set.
        # Resolving it reaches compare_panel_content, so it is resolved only
        # when the gate is actually being run; the drawing path never touches
        # it. A plain (already-resolved) value still works.
        alias = spec[2] if len(spec) > 2 else None
        keep = spec[3] if len(spec) > 3 else None
        if only and label not in only:
            continue
        if do_check:
            if callable(alias):
                alias = alias()
            rc |= check(label, original, build, alias=alias, keep_axes=keep)
        else:
            build()
    return rc
