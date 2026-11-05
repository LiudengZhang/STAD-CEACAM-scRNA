#!/usr/bin/env python3
"""Shared support for the supplementary panel drivers.

Supplementary panels S4 A-C, S5 E, S7 A-E and S9 A-E come from analysis
scripts under `04_Revision_Analyses/*/scripts/`. Each draw_*.py driver imports the
analysis module for its constants (without running `main()`), reads the CSVs
the analysis wrote, and draws the panel at print size. `--check` compares the
drawing with the module's own `_panel_*` function; it needs
10_Reproduction/compare_panel_content.py, which is not part of the release.

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
    """Import compare_panel_content (needed only for --check)."""
    try:
        import compare_panel_content as cpc
    except ModuleNotFoundError as e:                      # noqa: PERF203
        raise SystemExit(
            f"--check needs 10_Reproduction/compare_panel_content.py, which is "
            f"not on the path ({e}). It is a verification tool and is not part "
            f"of the code release. Drawing needs none of it.") from e
    return cpc

MM_TO_INCH = style.MM_TO_INCH

#: Axis labels and titles 7 pt, tick labels, legends and annotations 6 pt.
TYPE_SPEC = dict(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)

#: Paper left between the drawn ink and the edge of the panel canvas.
PAD_MM = 0.6

#: Declared label changes, used by --check.
LABELS = ROOT / "00_Config" / "shared" / "labels.py"

def aliases():
    """Declared label changes. Pass the function itself to `run()`; it is
    resolved only under --check."""
    return _cpc().load_aliases(LABELS)

def apply_style():
    """Put the type specification into matplotlib's rcParams."""
    return style.apply(**TYPE_SPEC)

def despine(ax):
    """Hide the top and right spines."""
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

def finish(ax):
    """Tick marks RULE_PT wide and 2 pt long, then despine."""
    ax.tick_params(axis="both", width=style.RULE_PT, length=2)
    despine(ax)

@contextlib.contextmanager
def _no_mkdir():
    """Disable Path.mkdir inside the block (one analysis module creates a
    directory at import)."""
    real = Path.mkdir
    Path.mkdir = lambda self, *a, **k: None
    try:
        yield
    finally:
        Path.mkdir = real

def analysis(rel):
    """Import an analysis module by path without running its main(). `rel` is
    relative to 04_Revision_Analyses."""
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
    """Return `path`, raising if it is missing or empty."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"{what} is missing: {p}\n    The driver reads the tables the "
            f"analysis already wrote and never recomputes them; run the "
            f"analysis script itself to produce this file.")
    if p.stat().st_size == 0:
        raise ValueError(f"{what} is empty: {p}")
    return p

def _clear_figure_title(fig, gap_mm=0.8):
    """Lower the axes until a figure-level title clears them by `gap_mm`."""
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
    """Set the margins from the rendered ink. Raises if ink leaves the canvas
    or covers the panel-letter cell."""
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
    """Compare the plotted values and strings of `redrawn` with `original`,
    with all file writes diverted. `alias` lists allowed label changes."""
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
    """Write a panel into its directory at print size. Raises if ink leaves
    the canvas or covers the panel-letter cell."""
    out = panel_dir(figure, panel, root) / stem
    over = style.overflow_mm(fig)
    if max(over) > 0.0:
        raise RuntimeError(
            f"{stem}: ink outside the canvas (left, right, bottom, top mm) "
            f"{over}.")
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
    """Driver entry point. `panels` maps label -> (build, original[, alias[,
    keep_axes]]): build() draws and saves, original() draws the analysis's own
    version for the check, alias is the declared label changes to allow,
    keep_axes the axis indices of the original the panel retains."""
    argv = sys.argv[1:] if argv is None else argv
    only = [a for a in argv if not a.startswith("-")]
    do_check = "--check" in argv
    rc = 0
    for label, spec in panels.items():
        build, original = spec[0], spec[1]
        # alias may be a callable, resolved only under --check.
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
