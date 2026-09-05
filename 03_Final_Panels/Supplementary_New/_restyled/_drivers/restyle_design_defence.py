#!/usr/bin/env python3
"""
Restyle S7C - design defence (Reviewer 1 points R1.0 and R1.3b).

Analysis: 04_Revision_Analyses/01_R1.3_Cohort_Pairing/scripts/design_defence.py

NOTHING IS RECOMPUTED HERE. The driver imports that module for its constants
and reads the two tables its `main()` already wrote. It never opens MoMac.h5ad
and it never reaches `acquired()`, which is where the bare `except Exception`
around the SCENIC read lives (design_defence.py:196). That handler turns a
failed AUCell read into a silently narrower table - `nfkb1_regulon` simply
stops being one of `cols` and the panel is drawn as though the column had never
been asked for. A restyle must not inherit that silence, so the driver asserts
the column is present and names the fallback if it is not. See `_frames()`.

    conda run -n Liudeng_Python_310 python restyle_design_defence.py
    conda run -n Liudeng_Python_310 python restyle_design_defence.py --check
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _driver_base as base
import panel_style_cns as style

A = base.analysis("01_R1.3_Cohort_Pairing/scripts/design_defence.py")
OUT = base.outputs_of(A)
FIG = "S7_Cohort_Statistics"

# Printed box, millimetres. Version A's box was 171 x 42; the two-line x label,
# the two-line title on the right axes and 7 pt covariate names need the height.
W, H = 171.0, 54.0

# Mark sizes, rescaled so they keep their size relative to the type.
# Version A drew this panel at 4x and the assembler fitted it at 0.2625; its
# smallest type was the 5 pt tick/covariate label, i.e. 5.25 pt on paper.
A_FIT = 0.2625
A_TYPE = 5.0 * A.SCALE * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2


def frames():
    """The two tables the analysis already wrote. No computation."""
    bal = pd.read_csv(base.require(
        OUT / "covariate_balance.csv", "S7C covariate balance"))
    d = pd.read_csv(base.require(
        OUT / "acquired_resistance_subgroup.csv",
        "S7C acquired-resistance subgroup"), index_col=0)

    # `cols` exactly as design_defence.acquired() derives it (line 199).
    cols = [c for c in ("il1b_fraction", "nfkb1_regulon") if c in d.columns]

    # The analysis wraps its SCENIC read in `except Exception` and, when it
    # fires, writes this CSV without `nfkb1_regulon` and carries on. That is a
    # missing input, not a narrower panel, and it must not pass quietly here.
    if "nfkb1_regulon" not in d.columns:
        raise RuntimeError(
            "acquired_resistance_subgroup.csv has no `nfkb1_regulon` column: "
            f"columns are {list(d.columns)}.\n"
            "    design_defence.acquired() reads the SCENIC AUCell matrix "
            "inside a bare `except Exception` (design_defence.py:184-197). "
            "When that read fails it prints '(regulon activity unavailable)' "
            "into the report, drops the column and still exits 0 - so this CSV "
            "is the visible symptom of a failed SCENIC read, not a design "
            "choice. Fix the input; do not restyle around it.")

    # `_panel` reads `bal.dropna(subset=["smd"])` and `bal["covariate"]`, and
    # `d["subgroup"]` / `d["il1b_fraction"]`. Fail here rather than in a plot.
    for col in ("smd", "covariate"):
        if col not in bal.columns:
            raise RuntimeError(f"covariate_balance.csv has no `{col}` column: "
                               f"columns are {list(bal.columns)}")
    for col in ("subgroup", "il1b_fraction"):
        if col not in d.columns:
            raise RuntimeError(f"acquired_resistance_subgroup.csv has no "
                               f"`{col}` column: columns are {list(d.columns)}")
    return dict(bal=bal, d=d, cols=cols)


# --------------------------------------------------------------------------
# S7C - covariate balance (left) and the acquired-resistance subgroup (right).
# Layout unchanged: 1 x 2.
# --------------------------------------------------------------------------
def draw_C(fr, save=True):
    style.apply()
    bal, d = fr["bal"], fr["d"]
    fig, axes = style.subplots_mm(W, H, 1, 2)

    ax = axes[0]
    b = bal.dropna(subset=["smd"])
    y = np.arange(len(b))
    ax.barh(y, b["smd"], color="#4d4d4d", edgecolor="#333", linewidth=0.4,
            height=0.6)
    for v in (-0.5, 0.5):
        ax.axvline(v, color="#B2182B", linestyle="--", linewidth=0.5)
    ax.axvline(0, color="#666", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(b["covariate"])
    ax.set_xlabel("Standardised mean difference\n(pre vs post group)")
    ax.set_xlim(-1.2, 1.2)
    ax.set_title("Covariate balance")

    ax = axes[1]
    if "il1b_fraction" in d.columns:
        groups = ["Sustained PR", "PR then SD (emerging resistance)"]
        colors = ["#2166AC", "#B2182B"]
        for i, (g, c) in enumerate(zip(groups, colors)):
            v = d.loc[d["subgroup"] == g, "il1b_fraction"].dropna() * 100
            ax.scatter([i] * len(v), v, s=40 * A.SCALE * AREA, c=c, zorder=3,
                       edgecolors="white", linewidths=0.5 * A.SCALE * MARK)
            if len(v):
                ax.plot([i - 0.22, i + 0.22], [v.mean()] * 2, color=c,
                        linewidth=1.5 * MARK)
        ax.set_xticks(range(2))
        ax.set_xticklabels(["Sustained\nPR", "PR then SD"])
        ax.set_ylabel("IL-1$\\beta$+ state\n(% of mono/macrophages)")
        ax.set_xlim(-0.5, 1.5)
        ax.set_title("Post-treatment responders only\n(descriptive, n = 3 vs 2)")

    for ax in axes:
        ax.tick_params(axis="both", width=0.6, length=2)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)

    style.margins_mm(fig, left=28, right=2, top=9, bottom=14, wspace=0.42)
    if save:
        base.save(fig, FIG, "S7_C", "S7_C_design_defence")
    return fig


def main():
    fr = frames()
    return base.run({
        "S7_C": (lambda: draw_C(fr),
                 lambda: A._panel(fr["bal"], fr["d"], fr["cols"])),
    })


if __name__ == "__main__":
    sys.exit(main())
