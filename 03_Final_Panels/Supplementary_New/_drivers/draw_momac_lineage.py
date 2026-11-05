#!/usr/bin/env python3
"""Supplementary Figure S9 A and B: monocyte versus macrophage lineage.

A: per-state mean monocyte and macrophage signature scores, read from
momac_lineage_scores.csv written by
04_Revision_Analyses/06_R1.7_MoMac_Lineage_Markers/scripts/momac_lineage.py
B: marker dot plot from MoMac.h5ad over 21 monocyte and macrophage genes.
Writes S9_MoMac_Identity_NFkB/S9_A and S9_B.

Run: python draw_momac_lineage.py [--check [S9_A|S9_B]]"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import scanpy as sc

import _driver_base as base
import panel_style_cns as style

base.apply_style()

A = base.analysis("06_R1.7_MoMac_Lineage_Markers/scripts/momac_lineage.py")
OUT = base.outputs_of(A)
FIG = "S9_MoMac_Identity_NFkB"

# Panel sizes in mm.
A_W, A_H = 58.0, 58.0
B_W, B_H = 108.1, 58.0

# Marker sizes for A, scaled to the type size.
A_FIT = 0.3100
A_TYPE = 4.5 * A.SCALE * A_FIT
A_MARK = (style.tick_pt() / A_TYPE) * A_FIT
A_AREA = A_MARK ** 2

#: Printed state names (as in Figure 4), by cluster.
SHOWN = {
    'C0_Mac_Classic_TREM2':           'Mac_TREM2',
    'C1_Mono_Classic_CD14':           'Mono_CD14',
    'C2_MoMac_Intermediate_HLA-DRA':  'MoMac_Inter',
    'C3_Mac_Inflam_IL1B':             'MoMac_IL1B',
    'C4_Mono_Alternative_CD16':       'Mono_CD16',
    'C5_Mac_Prolif_MKI67':            'Mac_Prolif',
    'C6_Mac_Metallothionein_MT1G':    'Mac_MT1G',
}

def frames():
    """The per-state summary table the analysis wrote."""
    summary = pd.read_csv(base.require(
        OUT / "momac_lineage_scores.csv", "S10A lineage scores"), index_col=0)
    for col in ("monocyte_score", "macrophage_score"):
        if col not in summary.columns:
            raise RuntimeError(f"momac_lineage_scores.csv has no `{col}` "
                               f"column: columns are {list(summary.columns)}")
    return dict(summary=summary)

#: The gene columns of S9 B, in order.
PUBLISHED_MONOCYTE = ["FCN1", "VCAN", "S100A8", "S100A9", "S100A12", "SELL",
                      "FCGR3A", "LYZ", "CSF3R"]
PUBLISHED_MACROPHAGE = ["CD68", "CD163", "MRC1", "C1QA", "C1QB", "C1QC",
                        "APOE", "APOC1", "TREM2", "MERTK", "MSR1", "SELENOP"]

_MATRIX = {}

def matrix():
    """MoMac.h5ad and the two gene lists, loaded once. Raises if a gene is
    missing from the object.
    """
    if not _MATRIX:
        ad = sc.read_h5ad(base.require(A.MOMAC_H5AD, "S10B MoMac matrix"))
        present = set(ad.var_names)
        want = PUBLISHED_MONOCYTE + PUBLISHED_MACROPHAGE
        gone = [g for g in want if g not in present]
        if gone:
            raise SystemExit(
                f"S9 B: the panel shows {len(want)} gene columns and "
                f"this object is missing {', '.join(gone)}. The panel cannot be "
                f"reproduced from it. Do not drop the column - use the object "
                f"the panel was drawn from.")
        _MATRIX["ad"] = ad
        _MATRIX["mono"] = list(PUBLISHED_MONOCYTE)
        _MATRIX["mac"] = list(PUBLISHED_MACROPHAGE)
    return _MATRIX

def draw_scatter(fr, save=True):
    """S9 A: state means on the monocyte/macrophage score plane."""
    base.apply_style()
    summary = fr["summary"]
    fig, ax = style.subplots_mm(A_W, A_H)
    for state, r in summary.iterrows():
        highlight = state == A.TARGET
        ax.scatter(r["monocyte_score"], r["macrophage_score"],
                   s=(70 if highlight else 45) * A.SCALE * A_AREA,
                   c="#B2182B" if highlight else "#4d4d4d",
                   edgecolors="white", linewidths=style.EDGE_PT,
                   zorder=3)
        # Labels start right of the point; _place_labels moves any that collide.
        label = SHOWN[str(state)]
        ax.annotate(label, (r["monocyte_score"], r["macrophage_score"]),
                    textcoords="offset points",
                    xytext=(6 * A.SCALE * A_MARK, 3 * A.SCALE * A_MARK),
                    ha="left", fontsize=style.tick_pt(),
                    color="#B2182B" if highlight else "#333333")
    lim = [min(summary["monocyte_score"].min(),
               summary["macrophage_score"].min()),
           max(summary["monocyte_score"].max(),
               summary["macrophage_score"].max())]
    pad = 0.12 * (lim[1] - lim[0])
    # Extra x room for the labels.
    ax.set_xlim(lim[0] - 3.5 * pad, lim[1] + 5.5 * pad)
    ax.set_ylim(lim[0] - pad, lim[1] + pad)
    ax.set_xlabel("Monocyte signature score")
    ax.set_ylabel("Macrophage signature score")
    base.finish(ax)
    base.fit(fig)
    _place_labels(fig, ax)
    if save:
        base.save(fig, FIG, "S9_A", "S9_A_momac_lineage_scores")
    return fig

#: Candidate label positions, in order tried: (dx, dy) in units of the 6-point
#: offset, and the alignment.
_LABEL_SLOTS = ((1, 0.5, "left", "bottom"), (-1, 0.5, "right", "bottom"),
                (1, -0.5, "left", "top"), (-1, -0.5, "right", "top"),
                (0, 1.2, "center", "bottom"), (0, -1.2, "center", "top"))

def _place_labels(fig, ax, gap_mm=0.5):
    """Give each state label the first slot in `_LABEL_SLOTS` that is `gap_mm`
    clear of every marker, tick label and earlier label; raise if none is."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pad = gap_mm / 25.4 * fig.dpi
    marks = []
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        if lab.get_text():
            bb = lab.get_window_extent(renderer=r)
            marks.append((bb.x0, bb.y0, bb.x1, bb.y1))
    for col in ax.collections:
        off = ax.transData.transform(col.get_offsets())
        for (x, y), s in zip(off, col.get_sizes()):
            rad = (s ** 0.5) / 72 * fig.dpi / 2
            marks.append((x - rad, y - rad, x + rad, y + rad))
    placed = []

    frame = ax.get_window_extent(renderer=r)

    def clear(bb):
        box = (bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad)
        # The label must also stay inside the plotting box.
        if box[0] < frame.x0 or box[3] > frame.y1 + 2 * pad:
            return False
        return all(box[2] < o[0] or box[0] > o[2] or box[3] < o[1] or box[1] > o[3]
                   for o in marks + placed)

    unit = 6 * A.SCALE * A_MARK
    for t in ax.texts:
        for dx, dy, ha, va in _LABEL_SLOTS:
            t.set_position((dx * unit, dy * unit))
            t.set_ha(ha); t.set_va(va)
            bb = t.get_window_extent(renderer=r)
            if clear(bb):
                placed.append((bb.x0, bb.y0, bb.x1, bb.y1))
                break
        else:
            raise RuntimeError(f"S9 A: no clear slot for the label {t.get_text()!r}")
    over = style.overflow_mm(fig)
    if max(over) > 0.0:
        raise RuntimeError(f"S9 A: a moved label left the canvas: {over}")

def _raise_to_floor(fig, pt):
    """Raise any text smaller than `pt` to `pt`."""
    for ax in fig.axes:
        ticks = ax.get_xticklabels() + ax.get_yticklabels()
        for t in [ax.title, ax.xaxis.label, ax.yaxis.label] + list(ax.texts) \
                + ticks:
            if t.get_text() and t.get_fontsize() < pt:
                t.set_fontsize(pt)
        if any(t.get_text() and t.get_fontsize() < pt for t in ticks):
            ax.tick_params(axis="both", labelsize=pt)

def draw_dotplot_live(save=True):
    """S9 B: marker dot plot (standard_scale="var") through _dotplot, with
    Monocyte/Macrophage gene-group brackets."""
    import _dotplot
    m = matrix()
    base.apply_style()
    style.setup_scanpy()

    mono, mac = m["mono"], m["mac"]
    genes = mono + mac
    ad = m["ad"].copy()
    states = sorted(set(ad.obs["minor_cell_state"].astype(str)))
    order = [SHOWN[s] for s in states]          # the cluster order, C0 first
    ad.obs["cell state"] = pd.Categorical(
        [SHOWN[s] for s in ad.obs["minor_cell_state"].astype(str)],
        categories=order)
    fr = _dotplot.frames(ad, genes, "cell state", order, path=A.MOMAC_H5AD)
    fig, _ = _dotplot.draw(
        fr, order, genes, w_mm=B_W, h_mm=B_H, left_mm=17.0, bottom_mm=12.0,
        top_mm=8.0, largest_mm=2.4, panel="S9 B",
        var_groups=[((0, len(mono) - 1), "Monocyte"),
                    ((len(mono), len(genes) - 1), "Macrophage")])
    _key_below(fig)
    _level_group_labels(fig, ["Monocyte", "Macrophage"])
    if save:
        base.save(fig, FIG, "S9_B", "S9_B_momac_lineage_dotplot")
    return fig

def _key_below(fig, key_h_mm=14.0):
    """Move B's key into a strip under the gene names and widen the matrix.
    Expects _dotplot's axes in order: matrix, bracket strip, spacer, then the
    four key axes."""
    W, H = B_W, B_H
    ax = fig.axes
    if len(ax) != 7:
        raise RuntimeError(f"S9 B: expected 7 axes from _dotplot, got {len(ax)}")
    main, group, spacer, keys = ax[0], ax[1], ax[2], ax[3:]
    spacer.set_visible(False)
    left = main.get_position().x0 * W
    right = W - 1.5
    top = main.get_position().y1 * H
    main.set_position([left / W, (12.0 + key_h_mm + 1.0) / H,
                       (right - left) / W, (top - 13.0 - key_h_mm) / H])
    gp = group.get_position()
    group.set_position([left / W, gp.y0, (right - left) / W, gp.height])
    x0 = keys[0].get_position().x0 * W
    size_top = max(k.get_position().y1 for k in keys[:2]) * H
    col_top = max(k.get_position().y1 for k in keys[2:]) * H

    def shift(a, dx, dy):
        p = a.get_position()
        a.set_position([p.x0 + dx / W, p.y0 + dy / H, p.width, p.height])
    for k in keys[:2]:
        shift(k, left - x0, (key_h_mm + 0.5) - size_top)
    for k in keys[2:]:
        shift(k, left + 23.0 - x0, (key_h_mm - 3.0) - col_top)


def _level_group_labels(fig, labels):
    """Set the var-group labels horizontal and the brackets at RULE_PT (scanpy
    exposes neither as a keyword); raise if a label then overlaps the other
    label or a bracket."""
    from matplotlib.patches import PathPatch
    texts, patches = [], []
    for ax in fig.axes:
        for t in ax.texts:
            if t.get_text() in labels:
                t.set_rotation(0)
                texts.append(t)
        for pt in ax.patches:
            if isinstance(pt, PathPatch):
                pt.set_linewidth(style.RULE_PT)
                patches.append(pt)
    if len(texts) != len(labels) or not patches:
        raise RuntimeError(f"S9 B: expected {len(labels)} var-group labels and "
                           f"a bracket patch, found {len(texts)} and {len(patches)}")
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    boxes_ = [t.get_window_extent(renderer=r) for t in texts]
    if boxes_[0].overlaps(boxes_[1]):
        raise RuntimeError("S9 B: the two group labels overlap when set upright")
    for t, bb in zip(texts, boxes_):
        for pt in patches:
            if bb.overlaps(pt.get_window_extent(renderer=r)):
                raise RuntimeError(f"S9 B: the label {t.get_text()!r} sits on the bracket")

def main():
    fr = frames()
    check = "--check" in sys.argv
    return base.run({
        # _panel_scatter reads only `summary`, so None is passed for `obs`.
        "S9_A": (lambda: draw_scatter(fr),
                 lambda: A._panel_scatter(fr["summary"], None)),
        "S9_B": (lambda: draw_dotplot_live(),
                 lambda: A._panel_dotplot(matrix()["ad"], matrix()["mono"],
                                          matrix()["mac"])),
    })

if __name__ == "__main__":
    sys.exit(main())
