#!/usr/bin/env python3
"""
Figure 2 panel A - epithelial UMAP coloured by minor cell state, with the nine
cluster labels drawn on the plot.

The panel is drawn at the millimetre rectangle it prints in, read from
03_Final_Panels/panel_rects.csv through 00_Config/slots.py, so the type size
set here is the type size printed. The margins are millimetres of paper.

  printed panel  Figure 2 A       (PROVENANCE.csv; the directory name agrees,
                                   but the letter was looked up, not inferred)

MARK
    The earlier drawing used a canvas four times the printed size, and the only
    body type it sets is the on-plot cluster label at 2.5 * SCALE - the rcParam
    it installs reaches nothing, because the panel has no ticks, no title, no
    axis labels and no legend. MARK carries the non-type point sizes across to
    the 1:1 canvas:

        MARK = style.tick_pt() / (SMALL_PT * SCALE)
        AREA = MARK ** 2

    The UMAP point area takes AREA and the leader line takes MARK. MARK is high
    here for the same reason it is low elsewhere: this panel's nominal type was
    unusually small, so holding the marks' size relative to the type moves them
    by less than in a panel set at 5 or 6 pt.

    The consequence, recorded because it is visible: the points go from about
    0.32 pt to 1.2 pt across while the panel keeps its printed width, so
    149,373 cells cover the square more densely than they do in the published
    panel. That is what holding the marks' size relative to the type produces
    when the type has to grow by a factor of four on an unchanged footprint. No
    cell, coordinate, colour or label changed.

THE LABEL POSITIONS ARE PINNED, AND THE PINS ARE MEASUREMENTS
    The earlier drawing places the nine labels with `adjust_text`, which repels
    a label by its *rendered* bounding box. Run unchanged on this canvas the
    algorithm returns a different answer: a 6 pt label on a 44.5 mm square is a
    far larger fraction of the axes than a 10 pt label on a 280 mm one, the
    repulsion is stronger, and the labels settle elsewhere - moving 55 plotted
    values and putting a label over a cluster it does not name. A redraw is a
    visualisation-only change, so that outcome is not available.

    `PINNED_LABELS` and `PINNED_LEADERS` are therefore not chosen; they are
    measured, by `measure_version_a_labels.py`, which runs the earlier script
    under `compare_panel_content.capture(writes="block")` and reads the nine
    final label positions and the six leader-line paths straight out of the
    figure it builds, in data coordinates. The raw measurement is kept beside
    this file in `label_pins_version_a.json`; the literals below were generated
    from it, not typed.

    Six leaders, not nine: adjustText draws one only where the target ends up
    at least `min_arrow_len = 5` display pixels outside the label's box, which
    on the earlier canvas is true for PTMA, KRT19, MUC5AC, Stem_SPINK4, MT1E
    and CD74 and false for CEACAM5/6, Chief_Like and Stem_TPX2. Each leader is
    recorded by its two endpoints as they were drawn - already clipped against
    the label box and shrunk - because that clipping is itself a function of
    the earlier canvas and cannot be recomputed here. They are stubs: the
    longest is 0.05 data units.

    The FancyArrowPatch is rebuilt with `shrinkA = shrinkB = 0` and no
    `patchA`, so its path in data coordinates is exactly the two measured
    endpoints and their midpoint - which is what a rad=0 arc3 connection is,
    and what the harvested path already was.

    Two guards run at draw time, against the data actually loaded, so that a
    stale pin cannot survive a change in the h5ad: every pinned label is nearer
    its own cluster's median than any other cluster's, and every leader's far
    endpoint is within 0.15 data units of the median of the cluster it belongs
    to. Both raise rather than warn.

THE NINE LABELS ARE THE SCRIPT'S OWN, NOT THE PAGE'S, AND THAT IS MEASURED
    The shipped page prints these labels with an "Epi_" prefix - "Epi_PTMA",
    "Epi_CEACAM5/6", "Epi_Chief_Like" - set by an assembler that is not on
    disk, at 5.47 pt. They cannot be carried at this figure's 6 pt floor inside
    this panel's published footprint, and the reason is arithmetic rather than
    layout.

    `sc.pl.umap` locks the axes to equal aspect, so the millimetres between two
    labels are set by the panel's height, which is its printed height, and no
    margin recovers any. Measured on the rendered panel: "Epi_CEACAM5/6" sets
    16.23 mm and "Epi_Chief_Like" 14.23 mm, their centres stand 12.87 mm apart,
    and half their widths sum to 15.23 mm - so the pair overlaps by 2.37 mm.
    Widening the box does not move them apart (tried: the horizontal fill moves
    by 0.3 mm and the overlap does not change) and shortening the box makes a
    second pair collide. In the script's own shorter form the same two labels
    sum to 10.99 mm of half-width against the same 12.87 mm, which clears by
    1.88 mm, and no pair on the panel touches.

    So the panel keeps the strings the script sets. The difference from the
    page is reported rather than absorbed, and nothing is shortened to hide it.

THE AXIS LABELS
    scanpy writes "UMAP1" and "UMAP2" onto the axes and `frameon=False` leaves
    them in place without ever drawing them - neither string reaches the saved
    SVG. They are hidden here rather than deleted, so the strings the panel
    declares are unchanged and the overflow check measures only ink that is
    actually drawn.
"""

import numpy as np
import scanpy as sc
import pandas as pd
from pathlib import Path
from matplotlib.patches import FancyArrowPatch
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import *                                       # noqa: E402,F403
import panel_style_cns as style                           # noqa: E402
import slots                                              # noqa: E402

BASE_DIR = Path(__file__).parent

SCALE = 4                           # the earlier canvas multiplier
SMALL_PT = 2.5                      # the earlier smallest body type

PANEL_LETTER = "A"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(2, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(2, PANEL_LETTER)
MARGIN = dict(left=2.9, right=0.6, top=0.6, bottom=0.6)

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
# Measured off the earlier drawing by measure_version_a_labels.py.
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
# as the earlier drawing drew it. A label absent here had no leader there.
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
    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
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

    # scanpy's axis labels are never drawn under frameon=False; hide them so
    # they are not measured either. The strings themselves are unchanged.
    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)

    # On-plot centroid labels: the positions adjust_text reached on the
    # earlier canvas, measured off it rather than recomputed here.
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
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    style.save_panel(fig, BASE_DIR / 'epithelial_umap_minor_states')
    print(f"Saved: {BASE_DIR / 'epithelial_umap_minor_states'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")


if __name__ == '__main__':
    main()
