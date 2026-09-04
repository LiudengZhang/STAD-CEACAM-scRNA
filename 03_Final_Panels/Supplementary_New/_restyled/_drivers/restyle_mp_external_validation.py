#!/usr/bin/env python3
"""
Restyle S8D - metaprogram external validation in PRJEB25780 (Reviewer 1, R1.4).

Analysis: 02_New_Analyses/03_R1.4_MP_Direction_PrePost/scripts/mp_external_validation.py

This panel was blocked at the end of stage 2 and is now unblocked. `_panel`
draws 45 per-sample points per programme out of `scores`, which `main()` built
in memory and never wrote down, so the panel could not be redrawn without
redoing the signature score - the log2, the within-cohort z, the >=5-gene
filter and the row mean. Recomputing it would have put the scoring code in two
files, which is the one thing these drivers exist to avoid. The analysis now
writes `mp_external_sample_scores.csv` beside its summary, so the driver reads
it like every other driver reads its table.

Nothing is recomputed here. In particular:

  scores  <- outputs/mp_external_sample_scores.csv, one column per programme,
             indexed by sample. `_panel` wants a dict of Series, which is what
             the columns of that frame are.
  v       <- outputs/mp_external_validation.csv, the per-programme summary that
             carries the P values printed in the titles.
  meta    <- the stored TIGER metadata TSV, reindexed to the sample list the
             scores file already fixes. The intersection with the expression
             matrix and the response_NR filter that produced that list are NOT
             re-run - they are read, in the analysis's own row order. The
             BayesPrism epithelial expression table is never opened.

    conda run -n Liudeng_Python_310 python restyle_mp_external_validation.py --check
    conda run -n Liudeng_Python_310 python restyle_mp_external_validation.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import _driver_base as base
import panel_style_cns as style

A = base.analysis("03_R1.4_MP_Direction_PrePost/scripts/mp_external_validation.py")
OUT = base.outputs_of(A)
FIG = "S8_CEACAM_Metaprogram"

W, H = 171.0, 54.0

# Version A drew 11.0 x 4.0 cm at SCALE 4 and the assembler fitted it into a
# 171 x 42 mm box at 0.2625, so its smallest type - tick labels at 5 * SCALE -
# printed at 5.25 pt.
A_FIT = 0.2625
A_TYPE = 5.0 * A.SCALE * A_FIT
MARK = (style.tick_pt() / A_TYPE) * A_FIT
AREA = MARK ** 2


def frames():
    """The tables the analysis already wrote. No score is recomputed."""
    sc = pd.read_csv(
        base.require(OUT / "mp_external_sample_scores.csv",
                     "S8D per-sample metaprogram scores"),
        index_col="sample")
    v = pd.read_csv(base.require(OUT / "mp_external_validation.csv",
                                 "S8D per-programme summary"))

    meta = pd.read_csv(base.require(A.TIGER_META, "PRJEB25780 metadata"),
                       sep="\t").set_index("sample_id")
    missing = sc.index.difference(meta.index)
    if len(missing):
        raise KeyError(
            f"{len(missing)} scored sample(s) are absent from the TIGER "
            f"metadata: {list(missing)[:5]}. The scores file and the metadata "
            f"have drifted apart; do not draw this panel until they agree.")
    # Reindexed to the sample list the analysis already settled on, in its row
    # order. This is a lookup, not the filter being re-run.
    meta = meta.loc[sc.index]

    # `_panel` splits on this column and on nothing else. If it ever stopped
    # holding exactly the two labels the analysis kept, the panel would quietly
    # draw fewer points rather than fail.
    labels = set(meta["response_NR"].dropna().unique())
    if not labels <= {"R", "N"}:
        raise ValueError(
            f"response_NR carries {sorted(labels)} for the scored samples; the "
            f"analysis keeps only R and N. The panel would silently drop rows.")

    scores = {c: sc[c] for c in sc.columns}
    return dict(v=v, scores=scores, meta=meta)


def draw_D(fr, save=True):
    style.apply()
    v, scores, meta = fr["v"], fr["scores"], fr["meta"]
    keys = sorted(scores)
    fig, axes = style.subplots_mm(W, H, 1, len(keys))
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(0)
    for ax, name in zip(axes, keys):
        s = scores[name]
        r = s[meta["response_NR"] == "R"].values
        nr = s[meta["response_NR"] == "N"].values
        bp = ax.boxplot([r, nr], positions=[0, 1], widths=0.55, patch_artist=True,
                        showfliers=False,
                        boxprops=dict(linewidth=0.8 * MARK),
                        whiskerprops=dict(linewidth=0.8 * MARK),
                        capprops=dict(linewidth=0.8 * MARK),
                        medianprops=dict(color="black", linewidth=1.2 * MARK))
        bp["boxes"][0].set_facecolor(A.COLOR_R); bp["boxes"][0].set_alpha(0.55)
        bp["boxes"][1].set_facecolor(A.COLOR_NR); bp["boxes"][1].set_alpha(0.55)
        for i, (vals, c) in enumerate(((r, A.COLOR_R), (nr, A.COLOR_NR))):
            ax.scatter(i + rng.uniform(-0.1, 0.1, len(vals)), vals,
                       s=7 * A.SCALE * AREA, c=c, zorder=3, edgecolors="white",
                       linewidths=0.3 * A.SCALE * MARK)
        p = v.loc[v["program"] == name, "p_two_tailed"]
        short = name.replace(" (stomach)", "").replace("NMF ", "")
        ax.set_title(f"{short}\nP = {p.iloc[0]:.3f}" if len(p) else short)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["R", "NR"])
        ax.tick_params(axis="both", width=0.6, length=2)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    axes[0].set_ylabel("Signature score\n(PRJEB25780 epithelium)")
    style.margins_mm(fig, left=17, right=2, top=10, bottom=8, wspace=0.48)
    if save:
        base.save(fig, FIG, "S8_D", "S8_D_metaprogram_external_validation")
    return fig


def main():
    fr = frames()
    return base.run({
        "S8_D": (lambda: draw_D(fr),
                 lambda: A._panel(fr["v"], fr["scores"], fr["meta"])),
    })


if __name__ == "__main__":
    sys.exit(main())
