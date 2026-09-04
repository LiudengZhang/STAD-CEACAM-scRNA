#!/usr/bin/env python3
"""
Figure 2, printed panel A, RESTYLED (Version B) - square UMAP of epithelial
cells coloured by minor cell state, with on-plot centroid labels.

Version A is
`03_Revised_Panels/Main_Figures/02_Figure_2/02_A/create_epithelial_umap.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md
allows: the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every cell, every category, every colour and every label string is Version A's,
and since 2026-09-01 every label *position* is Version A's too - see LABELS.

  printed panel  Figure 2 A     (PROVENANCE.csv; the directory letter happens
                                 to agree here - it was still looked up)
  printed rect   44.5 x 46.5 mm    (panel_rects.csv)
  Version B box  66.0 x 68.0 mm

MARK
    Version A drew at SCALE = 4 (7.0 x 7.0 cm x 4 = 280 x 280 mm). The only
    body type it sets is the on-plot centroid label at `2.5 * SCALE`; the
    `use_panel_style(font_pt=5)` rcParam reaches nothing, because the panel has
    no ticks, no title, no axis labels and no legend. So

        SCALE = 4, SMALL_PT = 2.5
        MARK  = style.tick_pt() / (SMALL_PT * SCALE) = 7 / 10 = 0.70
        AREA  = MARK ** 2 = 0.49

    MARK is high here for the same reason it is low elsewhere: this panel's
    nominal type was unusually small, so holding the marks' size relative to
    the type moves them by less than in a panel set at 5 or 6 pt.

    Two non-type sizes: the UMAP dot area `size = 4` takes AREA, and the
    leader line `lw = 0.4 * SCALE` takes MARK.

    JUDGEMENT CALL, recorded because it is visible: at MARK the dots go from
    about 0.32 pt to 1.4 pt across while the panel grows only 44.5 -> 66 mm, so
    149,373 cells cover the square more densely than they do in the published
    panel. That is what "hold the marks' size relative to the type" produces
    when the type has to grow 4.4x; matching the published dot *density*
    instead would need a 196 mm square, which is wider than the page. No cell,
    coordinate, colour or label changed.

LABELS - why they are pinned, and where the numbers come from
    Version A places the nine labels with `adjust_text`, which repels a label
    by its *rendered* bounding box. Run unchanged on this canvas the algorithm
    returns a different answer: a 7 pt label on a 66 mm square is a far larger
    fraction of the axes than a 10 pt label on Version A's 280 mm square, the
    repulsion is stronger, and the labels settle elsewhere. Measured on
    2026-09-01, that moved 55 plotted values and put a label over a cluster it
    does not name. A redraw is a visualisation-only change, so that outcome is
    not available.

    WHICH label, corrected. The stage-3 note recorded it as "Stem_SPINK4",
    landing over the KRT19/CEACAM5-6 region. Measurement says otherwise, and
    the measurement is in `check_label_layout.py --control`: mapping that
    render's nine text anchors back into data coordinates - its canvas is
    identical to this one, so the pinned panel calibrates the map to 7e-8 data
    units - it is CEACAM5/6 that crossed. It moved 4.37 data units into the gap
    beside the tan arm, with no dots at all beneath it, and its nearest
    centroid became Stem_SPINK4's (2.90) rather than its own (3.58).
    Stem_SPINK4 itself moved 0.93 units and stayed nearest its own cluster;
    MT1E (2.97 units) came to sit over the orange KRT19 body and MUC5AC (1.80)
    over the pink Chief_Like body. The decision to stop was right; the label it
    was attributed to was not.

    The fix is to stop running a canvas-dependent algorithm and to state the
    answer Version A reached. `PINNED_LABELS` and `PINNED_LEADERS` below are
    not chosen; they are *measured*, by
    `measure_version_a_labels.py`, which runs Version A's own script under
    `compare_panel_content.capture(writes="block")` - so nothing in
    `03_Revised_Panels/` is written - and reads the nine final label positions
    and the six leader-line paths straight out of the figure Version A builds,
    in data coordinates. The raw measurement is kept beside this file in
    `label_pins_version_a.json`. The literals below were generated from it, not
    typed.

    Six leaders, not nine: adjustText draws one only where the target ends up
    at least `min_arrow_len = 5` display pixels outside the label's box, which
    on Version A's canvas is true for PTMA, KRT19, MUC5AC, Stem_SPINK4, MT1E
    and CD74 and false for CEACAM5/6, Chief_Like and Stem_TPX2. Each leader is
    recorded by its two endpoints as Version A *drew* them - already clipped
    against the label box and shrunk - because that clipping is itself a
    function of Version A's canvas and cannot be recomputed here. They are
    stubs: the longest is 0.05 data units, about 0.15 mm on this canvas.
    Reproducing them exactly is what keeps the content gate honest; visually
    they are almost nothing, and at print size the (relatively larger) label
    box covers them.

    The FancyArrowPatch is rebuilt with `shrinkA = shrinkB = 0` and no
    `patchA`, so its path in data coordinates is exactly the two measured
    endpoints and their midpoint - which is what a rad=0 arc3 connection is,
    and what Version A's harvested path already was.

    Two guards run at draw time, against the data actually loaded, so a stale
    pin cannot survive a change in the h5ad:
      1. every pinned label is nearer its own cluster's median than any other
         cluster's - the "same label over the same cluster" property;
      2. every leader's far endpoint is within 0.15 data units of the median of
         the cluster it belongs to.
    Both raise rather than warn.

    No label needed a nudge: at 7 pt on this box no two label boxes overlap and
    none leaves the axes. Measured by `check_label_layout.py`.

OVERFLOW
    `style.overflow_mm` reports 2.54 mm on the left and bottom. It is a false
    positive: the artists are scanpy's "UMAP1" and "UMAP2" axis labels, which
    `frameon=False` leaves in place but never draws - neither `UMAP1` nor
    `UMAP2` appears anywhere in the saved SVG, in this version or in Version A.
    No drawn ink leaves the canvas. The margins below were set so that the label
    boxes, which are drawn, all sit inside it.
"""

import numpy as np
import scanpy as sc
import pandas as pd
from pathlib import Path
from matplotlib.patches import FancyArrowPatch
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402

BASE_DIR = Path(__file__).parent

SCALE = 4                           # Version A's canvas multiplier
SMALL_PT = 2.5                      # Version A's smallest body type

PRINTED_MM = (44.5, 46.5)           # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 66.0, 68.0
MARGIN = dict(left=7.0, right=7.0, top=2.0, bottom=3.5)

# Set2 palette (ColorBrewer) - ordered by category
COLORS = {
    'C0_Epi_PTMA': '#66C2A5',
    'C1_Epi_KRT19': '#FC8D62',
    'C2_Epi_CEACAM6': '#8DA0CB',
    'C3_Epi_Chief_Like_PGC': '#E78AC3',
    'C4_Epi_MUC5AC': '#A6D854',
    'C5_Epi_Stem_Like_TPX2': '#FFD92F',
    'C6_Epi_Stem_Like_SPINK4': '#E5C494',
    'C7_Epi_MT1E': '#B3B3B3',
    'C8_Epi_CD74': '#8DD3C7',
}

SHORT_LABELS = {
    'C0_Epi_PTMA': 'PTMA',
    'C1_Epi_KRT19': 'KRT19',
    'C2_Epi_CEACAM6': 'CEACAM5/6',
    'C3_Epi_Chief_Like_PGC': 'Chief_Like',
    'C4_Epi_MUC5AC': 'MUC5AC',
    'C5_Epi_Stem_Like_TPX2': 'Stem_TPX2',
    'C6_Epi_Stem_Like_SPINK4': 'Stem_SPINK4',
    'C7_Epi_MT1E': 'MT1E',
    'C8_Epi_CD74': 'CD74',
}

# ---------------------------------------------------------------------------
# Measured off Version A by measure_version_a_labels.py on 2026-09-01.
# label -> final position in UMAP data coordinates. Do not hand-edit.
# ---------------------------------------------------------------------------
PINNED_LABELS = {
    'PTMA':        (-1.3673487686183599, 7.648700458808822),
    'KRT19':       (-0.6530729121408712, 2.522255165382308),
    'CEACAM5/6':   (1.5870579527498947, 1.1080370972594444),
    'Chief_Like':  (-4.084364528786972, 0.3954561183890526),
    'MUC5AC':      (-0.6054549738708204, -3.412629859642106),
    'Stem_TPX2':   (-2.0646697236417086, 11.30056164960472),
    'Stem_SPINK4': (4.581559543478653, 3.3731412718733953),
    'MT1E':        (0.4446967487840592, 4.668013317390365),
    'CD74':        (2.140957321673275, 8.887738926216045),
}

# label -> (start, end) of the leader line, in UMAP data coordinates, exactly
# as Version A drew it. A label absent here had no leader in Version A.
PINNED_LEADERS = {
    'PTMA':        ((-0.8578974478527783, 7.901319242287869),
                    (-0.851193573202325, 7.904643455200397)),
    'KRT19':       ((-0.0970986603142503, 2.7789300599185704),
                    (-0.09501596318800054, 2.779891572189104)),
    'MUC5AC':      ((0.11218166749815417, -3.152955187494662),
                    (0.12065659524927241, -3.1498885602373825)),
    'Stem_SPINK4': ((5.3880739796948145, 3.6500102191982258),
                    (5.406764984130859, 3.6564266681671125)),
    'MT1E':        ((0.9270839981257755, 4.92634788467403),
                    (0.9736745953559858, 4.951298713684082)),
    'CD74':        ((2.608670181094002, 9.14865136248212),
                    (2.648776054382324, 9.171024322509766)),
}

LEADER_TOL = 0.15                   # data units, guard 2


def _guard_pins(centroids):
    """Raise unless the pins still describe the data that was just loaded."""
    names = list(centroids)
    for name, (lx, ly) in PINNED_LABELS.items():
        d = {n: float(np.hypot(lx - centroids[n][0], ly - centroids[n][1]))
             for n in names}
        nearest = min(d, key=d.get)
        if nearest != name:
            raise SystemExit(
                f"pinned label {name!r} at ({lx:.4f}, {ly:.4f}) is nearer the "
                f"{nearest!r} cluster ({d[nearest]:.4f}) than its own "
                f"({d[name]:.4f}). The pins no longer describe this data.")
    for name, (_, end) in PINNED_LEADERS.items():
        cx, cy = centroids[name]
        d = float(np.hypot(end[0] - cx, end[1] - cy))
        if d > LEADER_TOL:
            raise SystemExit(
                f"leader for {name!r} ends {d:.4f} data units from that "
                f"cluster's median, over the {LEADER_TOL} tolerance.")
    print(f"  pins guard: 9/9 labels nearest their own cluster; "
          f"6/6 leaders within {LEADER_TOL} of their cluster median")


def main():
    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.5f}")

    print("Loading data...")
    adata = sc.read_h5ad(EPITHELIAL_H5AD)                 # noqa: F405
    print(f"Loaded {adata.n_obs} cells")

    # Rename categories to short labels for legend
    orig_cats = list(adata.obs['minor_cell_state'].cat.categories)
    color_list = [COLORS[c] for c in orig_cats]
    adata.obs['minor_cell_state'] = adata.obs['minor_cell_state'].cat.rename_categories(SHORT_LABELS)
    adata.uns['minor_cell_state_colors'] = color_list

    print("=== Color mapping ===")
    for cat, color in zip(adata.obs['minor_cell_state'].cat.categories, color_list):
        print(f"  {cat}: {color}")

    fig, ax = style.subplots_mm(PANEL_W_MM, PANEL_H_MM)
    style.margins_mm(fig, **MARGIN)

    sc.pl.umap(
        adata,
        color='minor_cell_state',
        ax=ax,
        show=False,
        legend_loc='none',
        title='',
        frameon=False,
        size=4 * AREA,
        alpha=0.6,
    )

    # Remove legend if scanpy created one anyway
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # On-plot centroid labels. Version A ran adjust_text here; this draws the
    # positions adjust_text reached on Version A's canvas, measured off it.
    coords = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    coords['cluster'] = adata.obs['minor_cell_state'].values

    centroids = {}
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        mask = coords['cluster'] == cluster_name
        centroids[cluster_name] = (float(coords.loc[mask, 'UMAP1'].median()),
                                   float(coords.loc[mask, 'UMAP2'].median()))
    _guard_pins(centroids)

    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        cx, cy = PINNED_LABELS[cluster_name]
        ax.text(cx, cy, cluster_name, fontsize=style.tick_pt(),
                ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                          alpha=0.8, edgecolor='none', linewidth=0))

    # The leader lines, in the same order the labels are drawn in.
    for cluster_name in adata.obs['minor_cell_state'].cat.categories:
        if cluster_name not in PINNED_LEADERS:
            continue
        posA, posB = PINNED_LEADERS[cluster_name]
        ax.add_patch(FancyArrowPatch(
            posA=posA, posB=posB,
            arrowstyle='-', connectionstyle='arc3,rad=0.0',
            shrinkA=0, shrinkB=0,
            color='0.4', lw=0.4 * SCALE * MARK,
            transform=ax.transData))

    # Rasterize scatter dots (keeps axes/legend as vectors, dots as embedded raster)
    for coll in ax.collections:
        coll.set_rasterized(True)

    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")

    style.save_panel(fig, BASE_DIR / 'epithelial_umap_minor_states')
    print(f"Saved: {BASE_DIR / 'epithelial_umap_minor_states'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
